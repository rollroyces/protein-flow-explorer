"""Tests for :mod:`pfe_design.metrics`.

The metrics are geometric, so every test here is built around a physical
property that has a known answer: a rigid motion must be undone exactly, a
translated interface must move the iRMSD by exactly the translation length,
an isolated sphere must have the analytic surface area, and two atoms that
touch must bury some of each other's surface.
"""

from __future__ import annotations

import math

import pytest
import torch

from pfe_design.config import MetricConfig
from pfe_design.metrics import (
    apply_transform,
    buried_surface_area,
    clash_penalty,
    clash_report,
    interface_rmsd,
    ligand_rmsd,
    sasa,
    superpose,
    vdw_radii,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _rotation(seed: int = 0) -> torch.Tensor:
    """A proper rotation matrix (det = +1) drawn deterministically."""
    gen = torch.Generator().manual_seed(seed)
    a = torch.randn(3, 3, generator=gen, dtype=torch.float64)
    u, _s, vh = torch.linalg.svd(a)
    r = u @ vh
    if float(torch.det(r)) < 0.0:
        r = -r
    return r


def _points(n: int, seed: int = 0, scale: float = 1.0) -> torch.Tensor:
    """``n`` deterministic pseudo-random points."""
    gen = torch.Generator().manual_seed(seed)
    return torch.randn(n, 3, generator=gen, dtype=torch.float64) * scale


def _complex(
    n_receptor: int = 12,
    n_ligand: int = 6,
    seed: int = 0,
) -> tuple[torch.Tensor, torch.Tensor]:
    """A synthetic complex and the mask of its interface residues."""
    coords = _points(n_receptor + n_ligand, seed=seed, scale=8.0)
    iface = torch.zeros(coords.shape[0], dtype=torch.bool)
    iface[n_receptor:] = True
    return coords, iface


# ---------------------------------------------------------------------------
# superpose
# ---------------------------------------------------------------------------


def test_superpose_recovers_a_pure_rotation_exactly():
    """A rotated copy of a point set superposes back with ~zero residual."""
    torch.manual_seed(0)
    mobile = _points(12, seed=1)
    r_true = _rotation(seed=2)
    target = mobile @ r_true.T

    aligned, rotation, _translation = superpose(mobile, target)

    assert aligned.shape == mobile.shape
    assert torch.max((aligned - target).abs()).item() < 1e-10
    # The returned rotation must undo the applied one.
    assert torch.max((rotation - r_true).abs()).item() < 1e-10


def test_superpose_removes_a_pure_translation():
    """A rigid shift is absorbed by the translation and the rotation stays identity."""
    torch.manual_seed(0)
    mobile = _points(10, seed=3)
    shift = torch.tensor([1.5, -2.0, 0.75], dtype=torch.float64)
    target = mobile + shift

    aligned, rotation, translation = superpose(mobile, target)

    assert torch.max((aligned - target).abs()).item() < 1e-10
    assert torch.max((rotation - torch.eye(3, dtype=torch.float64)).abs()).item() < 1e-10
    assert torch.allclose(translation, shift, atol=1e-10)


def test_superpose_returns_an_orthonormal_rotation():
    """The Kabsch rotation is orthonormal: R R^T = I."""
    torch.manual_seed(0)
    mobile = _points(15, seed=4)
    target = _points(15, seed=5)

    _aligned, rotation, _translation = superpose(mobile, target)

    identity = torch.eye(3, dtype=torch.float64)
    assert torch.allclose(rotation @ rotation.T, identity, atol=1e-10)
    assert torch.allclose(rotation.T @ rotation, identity, atol=1e-10)


def test_superpose_rotation_has_determinant_plus_one():
    """A superposition must never introduce a reflection (det = -1)."""
    torch.manual_seed(0)
    mobile = _points(15, seed=6)
    target = _points(15, seed=7)

    _aligned, rotation, _translation = superpose(mobile, target)

    assert float(torch.det(rotation)) == pytest.approx(1.0, abs=1e-10)


def test_superpose_supports_batched_input():
    """Batched (B, N, 3) input gives per-structure rotations and translations."""
    torch.manual_seed(0)
    mobile = _points(9, seed=8)
    r_true = _rotation(seed=9)
    shifted = mobile + torch.tensor([2.0, 0.0, -1.0], dtype=torch.float64)
    mobile_batch = torch.stack([mobile, shifted])
    target_batch = torch.stack([mobile @ r_true.T, shifted + 1.0])

    aligned, rotation, translation = superpose(mobile_batch, target_batch)

    assert aligned.shape == (2, 9, 3)
    assert rotation.shape == (2, 3, 3)
    assert translation.shape == (2, 3)
    for b in range(2):
        assert torch.max((aligned[b] - target_batch[b]).abs()).item() < 1e-10
    assert torch.allclose(torch.det(rotation), torch.ones(2, dtype=torch.float64), atol=1e-10)


def test_superpose_batched_matches_unbatched():
    """Batching must not change the answer for either structure."""
    torch.manual_seed(0)
    mobile = _points(11, seed=10)
    target = _points(11, seed=11)

    aligned_b, rotation_b, translation_b = superpose(
        torch.stack([mobile, mobile]), torch.stack([target, target])
    )
    aligned_u, rotation_u, translation_u = superpose(mobile, target)

    assert torch.allclose(aligned_b[0], aligned_u, atol=1e-10)
    assert torch.allclose(aligned_b[1], aligned_u, atol=1e-10)
    assert torch.allclose(rotation_b[0], rotation_u, atol=1e-10)
    assert torch.allclose(translation_b[1], translation_u, atol=1e-10)


def test_superpose_rejects_fewer_than_three_points():
    """Two points do not define a rotation, so the fit must be refused."""
    torch.manual_seed(0)
    pts = _points(2, seed=12)
    with pytest.raises(ValueError):
        superpose(pts, pts)


def test_superpose_rejects_mismatched_shapes():
    """Mobile and target must describe the same set of atoms."""
    torch.manual_seed(0)
    with pytest.raises(ValueError):
        superpose(_points(5, seed=13), _points(4, seed=14))


# ---------------------------------------------------------------------------
# apply_transform
# ---------------------------------------------------------------------------


def test_apply_transform_reproduces_the_aligned_coordinates():
    """apply_transform(coords, R, t) must equal the superposed output."""
    torch.manual_seed(0)
    mobile = _points(12, seed=15)
    target = mobile @ _rotation(seed=16).T + torch.tensor([1.0, 2.0, 3.0], dtype=torch.float64)

    aligned, rotation, translation = superpose(mobile, target)
    manual = apply_transform(mobile, rotation, translation)

    assert torch.allclose(manual, aligned, atol=1e-12)
    assert torch.allclose(manual, target, atol=1e-10)


def test_apply_transform_preserves_interatomic_distances():
    """A rigid transform cannot change any inter-atomic distance."""
    torch.manual_seed(0)
    coords = _points(10, seed=17)
    rotation = _rotation(seed=18)
    translation = torch.tensor([-3.0, 0.5, 2.0], dtype=torch.float64)

    moved = apply_transform(coords, rotation, translation)
    before = torch.cdist(coords, coords)
    after = torch.cdist(moved, moved)

    assert torch.allclose(before, after, atol=1e-10)


def test_apply_transform_supports_batched_input():
    """Batched coordinates are transformed by their own rotation and shift."""
    torch.manual_seed(0)
    coords = _points(8, seed=19)
    rotation = _rotation(seed=20)
    batch = torch.stack([coords, coords + 4.0])
    rotations = torch.stack([rotation, torch.eye(3, dtype=torch.float64)])
    translations = torch.stack(
        [torch.zeros(3, dtype=torch.float64), torch.tensor([1.0, 1.0, 1.0], dtype=torch.float64)]
    )

    moved = apply_transform(batch, rotations, translations)

    assert moved.shape == (2, 8, 3)
    assert torch.allclose(moved[0], coords @ rotation.T, atol=1e-12)
    assert torch.allclose(moved[1], coords + 4.0 + 1.0, atol=1e-12)


# ---------------------------------------------------------------------------
# interface_rmsd
# ---------------------------------------------------------------------------


def test_interface_rmsd_equals_the_known_interface_displacement():
    """Fitting the receptor leaves the interface shift untouched, so iRMSD = |shift|."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=21)
    shift = torch.tensor([0.3, -0.4, 0.5], dtype=torch.float64)
    mobile = coords.clone()
    mobile[iface] += shift

    value = interface_rmsd(mobile, coords, iface, superposition="receptor")

    assert float(value) == pytest.approx(float(shift.norm()), rel=1e-9)


def test_interface_rmsd_with_interface_superposition_removes_the_shift():
    """Fitting on the interface itself absorbs a pure translation -> iRMSD ~ 0."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=22)
    shift = torch.tensor([0.3, -0.4, 0.5], dtype=torch.float64)
    mobile = coords.clone()
    mobile[iface] += shift

    on_receptor = float(interface_rmsd(mobile, coords, iface, superposition="receptor"))
    on_interface = float(interface_rmsd(mobile, coords, iface, superposition="interface"))

    assert on_interface == pytest.approx(0.0, abs=1e-9)
    assert on_interface < 1e-6 * max(on_receptor, 1.0)


def test_interface_rmsd_is_zero_for_identical_structures():
    """A perfect prediction scores exactly zero."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=23)

    value = interface_rmsd(coords, coords, iface)

    assert float(value) == pytest.approx(0.0, abs=1e-12)


def test_interface_rmsd_returns_zero_for_an_empty_interface():
    """No interface atoms means there is nothing to measure."""
    torch.manual_seed(0)
    coords, _iface = _complex(seed=24)
    empty = torch.zeros(coords.shape[0], dtype=torch.bool)

    value = interface_rmsd(coords, coords, empty)

    assert float(value) == pytest.approx(0.0, abs=1e-12)


def test_interface_rmsd_accepts_unbatched_input():
    """Unbatched input returns a scalar tensor, not a one-element vector."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=25)
    shift = torch.tensor([0.0, 0.0, 1.0], dtype=torch.float64)
    mobile = coords.clone()
    mobile[iface] += shift

    value = interface_rmsd(mobile, coords, iface)

    assert value.dim() == 0
    assert float(value) == pytest.approx(1.0, rel=1e-9)


def test_interface_rmsd_batched_returns_one_value_per_structure():
    """A batch of two structures yields two iRMSD values."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=26)
    shift = torch.tensor([0.0, 2.0, 0.0], dtype=torch.float64)
    mobile = coords.clone()
    mobile[iface] += shift

    values = interface_rmsd(torch.stack([mobile, coords]), torch.stack([coords, coords]), iface)

    assert values.shape == (2,)
    assert float(values[0]) == pytest.approx(2.0, rel=1e-9)
    assert float(values[1]) == pytest.approx(0.0, abs=1e-9)


def test_interface_rmsd_ignores_residues_outside_the_valid_mask():
    """Unmodelled residues must not dilute the interface RMSD."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=27)
    shift = torch.tensor([0.3, -0.4, 0.5], dtype=torch.float64)
    mobile = coords.clone()
    mobile[iface] += shift
    mobile[0] = coords[0]  # this interface atom did not move

    diluted = float(interface_rmsd(mobile, coords, iface))
    valid = torch.ones(coords.shape[0], dtype=torch.bool)
    valid[0] = False
    masked = float(interface_rmsd(mobile, coords, iface, valid_mask=valid))

    assert masked == pytest.approx(float(shift.norm()), rel=1e-9)
    assert diluted < masked


def test_interface_rmsd_rejects_an_unknown_superposition_mode():
    """Only 'receptor' and 'interface' are meaningful superposition modes."""
    torch.manual_seed(0)
    coords, iface = _complex(seed=28)
    with pytest.raises(ValueError):
        interface_rmsd(coords, coords, iface, superposition="sidechain")


# ---------------------------------------------------------------------------
# ligand_rmsd
# ---------------------------------------------------------------------------


def test_ligand_rmsd_equals_the_known_ligand_displacement():
    """The ligand is carried by the receptor fit, so its displacement survives."""
    torch.manual_seed(0)
    coords = _points(18, seed=29, scale=8.0)
    receptor = torch.zeros(18, dtype=torch.bool)
    receptor[:12] = True
    ligand = ~receptor
    shift = torch.tensor([0.0, 0.0, 1.5], dtype=torch.float64)
    mobile = coords.clone()
    mobile[ligand] += shift

    value = ligand_rmsd(mobile, coords, receptor, ligand)

    assert float(value) == pytest.approx(1.5, rel=1e-9)


def test_ligand_rmsd_is_zero_for_a_rigid_motion_of_the_whole_complex():
    """If everything moves together, the receptor fit removes the motion."""
    torch.manual_seed(0)
    coords = _points(18, seed=30, scale=8.0)
    receptor = torch.zeros(18, dtype=torch.bool)
    receptor[:12] = True
    ligand = ~receptor
    moved = coords @ _rotation(seed=31).T + torch.tensor([2.0, -1.0, 0.5], dtype=torch.float64)

    value = ligand_rmsd(moved, coords, receptor, ligand)

    assert float(value) == pytest.approx(0.0, abs=1e-9)


def test_ligand_rmsd_returns_zero_for_an_empty_ligand_mask():
    """Without a ligand there is nothing to measure."""
    torch.manual_seed(0)
    coords = _points(12, seed=32, scale=8.0)
    receptor = torch.ones(12, dtype=torch.bool)
    empty = torch.zeros(12, dtype=torch.bool)

    value = ligand_rmsd(coords, coords, receptor, empty)

    assert float(value) == pytest.approx(0.0, abs=1e-12)


def test_ligand_rmsd_batched_matches_unbatched():
    """Batching must not change the ligand RMSD."""
    torch.manual_seed(0)
    coords = _points(18, seed=33, scale=8.0)
    receptor = torch.zeros(18, dtype=torch.bool)
    receptor[:12] = True
    ligand = ~receptor
    shift = torch.tensor([1.0, 1.0, 0.0], dtype=torch.float64)
    mobile = coords.clone()
    mobile[ligand] += shift

    batched = ligand_rmsd(
        torch.stack([mobile, coords]), torch.stack([coords, coords]), receptor, ligand
    )
    unbatched = ligand_rmsd(mobile, coords, receptor, ligand)

    assert batched.shape == (2,)
    assert float(batched[0]) == pytest.approx(float(unbatched), rel=1e-9)
    assert float(batched[1]) == pytest.approx(0.0, abs=1e-9)


# ---------------------------------------------------------------------------
# sasa
# ---------------------------------------------------------------------------


def test_sasa_of_a_single_atom_matches_the_analytic_sphere_area():
    """An isolated atom exposes exactly 4*pi*(r + probe)^2."""
    torch.manual_seed(0)
    radius = 1.7
    probe = 1.4
    coords = torch.zeros(1, 3, dtype=torch.float64)
    radii = torch.tensor([radius], dtype=torch.float64)

    areas = sasa(coords, radii, probe=probe, n_points=100)

    expected = 4.0 * math.pi * (radius + probe) ** 2
    assert float(areas[0]) == pytest.approx(expected, rel=0.02)


def test_sasa_of_two_distant_atoms_is_the_sum_of_the_two():
    """Far apart atoms cannot occlude each other, so the areas add up."""
    torch.manual_seed(0)
    radius = 1.7
    coords = torch.tensor([[0.0, 0.0, 0.0], [30.0, 0.0, 0.0]], dtype=torch.float64)
    radii = torch.full((2,), radius, dtype=torch.float64)

    together = sasa(coords, radii, n_points=100)
    alone_a = sasa(coords[:1], radii[:1], n_points=100)
    alone_b = sasa(coords[1:], radii[1:], n_points=100)

    assert float(together.sum()) == pytest.approx(
        float(alone_a.sum()) + float(alone_b.sum()), rel=1e-6
    )


def test_sasa_of_two_contacting_atoms_is_smaller_than_the_sum():
    """Atoms in contact bury part of each other's surface."""
    torch.manual_seed(0)
    radius = 1.7
    coords = torch.tensor([[0.0, 0.0, 0.0], [2 * radius, 0.0, 0.0]], dtype=torch.float64)
    radii = torch.full((2,), radius, dtype=torch.float64)

    together = sasa(coords, radii, n_points=100)
    alone_a = sasa(coords[:1], radii[:1], n_points=100)
    alone_b = sasa(coords[1:], radii[1:], n_points=100)

    assert float(together.sum()) < 0.9 * (float(alone_a.sum()) + float(alone_b.sum()))
    assert float(together.sum()) > 0.0


def test_sasa_is_invariant_under_a_rigid_translation():
    """Moving the whole structure cannot change how much surface is exposed."""
    torch.manual_seed(0)
    coords = _points(6, seed=34, scale=4.0)
    radii = torch.full((6,), 1.7, dtype=torch.float64)

    before = sasa(coords, radii, n_points=100)
    after = sasa(coords + torch.tensor([13.0, -7.0, 2.0], dtype=torch.float64), radii, n_points=100)

    assert torch.allclose(before, after, rtol=1e-8, atol=1e-8)


def test_sasa_rejects_radii_that_do_not_match_the_coordinates():
    """The radius vector must have one entry per atom."""
    torch.manual_seed(0)
    with pytest.raises(ValueError):
        sasa(_points(4, seed=35), torch.full((3,), 1.7, dtype=torch.float64))


# ---------------------------------------------------------------------------
# buried_surface_area
# ---------------------------------------------------------------------------


def test_buried_surface_area_is_zero_for_distant_groups():
    """Groups that never touch bury no surface at all."""
    torch.manual_seed(0)
    radius = 1.7
    coords = torch.tensor([[0.0, 0.0, 0.0], [40.0, 0.0, 0.0]], dtype=torch.float64)
    radii = torch.full((2,), radius, dtype=torch.float64)

    result = buried_surface_area(coords, radii, [0], [1])

    assert result.bsa == pytest.approx(0.0, abs=1e-6)


def test_buried_surface_area_is_positive_for_groups_in_contact():
    """Two atoms in contact bury a clearly positive amount of surface."""
    torch.manual_seed(0)
    radius = 1.7
    coords = torch.tensor([[0.0, 0.0, 0.0], [2 * radius, 0.0, 0.0]], dtype=torch.float64)
    radii = torch.full((2,), radius, dtype=torch.float64)

    result = buried_surface_area(coords, radii, [0], [1])

    assert result.bsa > 10.0
    assert result.bsa < result.sasa_free_a + result.sasa_free_b


def test_buried_surface_area_satisfies_its_definition():
    """BSA = SASA(free A) + SASA(free B) - SASA(complex)."""
    torch.manual_seed(0)
    group_a = _points(8, seed=36, scale=4.0)
    group_b = group_a + torch.tensor([8.0, 0.0, 0.0], dtype=torch.float64)
    coords = torch.cat([group_a, group_b])
    radii = torch.full((16,), 1.7, dtype=torch.float64)

    result = buried_surface_area(coords, radii, list(range(8)), list(range(8, 16)))

    assert result.bsa == pytest.approx(
        result.sasa_free_a + result.sasa_free_b - result.sasa_complex, rel=1e-9
    )
    assert result.bsa > 0.0


def test_buried_surface_area_reports_one_value_per_atom():
    """The per-atom tensor covers every atom of the complex."""
    torch.manual_seed(0)
    coords = _points(11, seed=37, scale=4.0)
    radii = torch.full((11,), 1.7, dtype=torch.float64)

    result = buried_surface_area(coords, radii, [0, 1, 2], [3, 4, 5, 6, 7, 8, 9, 10])

    assert result.per_atom_complex.shape == (11,)
    assert float(result.per_atom_complex.sum()) == pytest.approx(result.sasa_complex, rel=1e-9)


# ---------------------------------------------------------------------------
# clash_report
# ---------------------------------------------------------------------------


def test_clash_report_detects_an_overlapping_pair():
    """Two atoms closer than the threshold are a clash."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=torch.float64)

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    assert report.n_clashes == 1
    assert report.n_atoms == 2
    assert report.max_overlap == pytest.approx(1.0, rel=1e-9)
    assert report.penalty == pytest.approx(1.0, rel=1e-9)


def test_clash_report_ignores_a_separated_pair():
    """Atoms further apart than the threshold are not a clash."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=torch.float64)

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    assert report.n_clashes == 0
    assert report.penalty == 0.0
    assert report.max_overlap == 0.0
    assert report.pairs.shape == (0, 2)


def test_clash_report_penalty_is_the_sum_of_squared_overlaps():
    """The penalty is the sum of the squared overlaps of every clashing pair."""
    torch.manual_seed(0)
    # Two well-separated clashes: atoms 0/1 overlap by 1.0 A and 2/3 by 0.5 A.
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.5, 0.0, 0.0], [5.0, 0.0, 0.0]], dtype=torch.float64
    )

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    expected_overlaps = [1.0, 0.5]
    assert report.n_clashes == len(expected_overlaps)
    assert report.penalty == pytest.approx(sum(o ** 2 for o in expected_overlaps), rel=1e-9)


def test_clash_report_max_overlap_is_the_largest_overlap():
    """max_overlap reports the single worst pair, not the average."""
    torch.manual_seed(0)
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [0.0, 1.5, 0.0]], dtype=torch.float64
    )

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    assert report.max_overlap == pytest.approx(1.5, rel=1e-9)


def test_clash_report_pairs_reports_the_offending_indices():
    """The reported index pairs point at the atoms that clash."""
    torch.manual_seed(0)
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [8.0, 0.0, 0.0], [8.5, 0.0, 0.0]], dtype=torch.float64
    )

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    assert report.pairs.shape == (1, 2)
    assert report.pairs[0].tolist() == [1, 2]


def test_clash_report_ignores_intra_residue_pairs_when_residue_indices_are_given():
    """A covalent bond is short, so bonded atoms must not be reported as clashes."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=torch.float64)
    residue_index = torch.tensor([7, 7], dtype=torch.long)

    report = clash_report(coords, MetricConfig(clash_distance=2.0), residue_index=residue_index)

    assert report.n_clashes == 0
    assert report.penalty == 0.0


def test_clash_report_counts_intra_residue_pairs_without_residue_indices():
    """Without residue indices nothing is excluded, so the short pair clashes."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=torch.float64)

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    assert report.n_clashes == 1
    assert report.penalty == pytest.approx(1.0, rel=1e-9)


def test_clash_report_vdw_mode_uses_the_van_der_waals_radii():
    """In vdw mode the threshold is the summed radii minus the tolerance."""
    torch.manual_seed(0)
    cfg = MetricConfig(clash_mode="vdw", vdw_tolerance=0.4)
    radii = vdw_radii(["C", "C"])
    limit = 2 * 1.70 - 0.4

    clashing = clash_report(
        torch.tensor([[0.0, 0.0, 0.0], [limit - 0.5, 0.0, 0.0]], dtype=torch.float64), cfg, radii
    )
    clear = clash_report(
        torch.tensor([[0.0, 0.0, 0.0], [limit + 0.5, 0.0, 0.0]], dtype=torch.float64), cfg, radii
    )

    assert clashing.n_clashes == 1
    assert clashing.max_overlap == pytest.approx(0.5, rel=1e-6)
    assert clashing.penalty == pytest.approx(0.25, rel=1e-6)
    assert clear.n_clashes == 0


def test_clash_report_vdw_mode_requires_radii():
    """The vdw threshold cannot be computed without radii."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=torch.float64)
    with pytest.raises(ValueError):
        clash_report(coords, MetricConfig(clash_mode="vdw"))


def test_clash_report_counts_the_number_of_atoms_analysed():
    """The report records how many atoms were considered."""
    torch.manual_seed(0)
    coords = _points(9, seed=38, scale=5.0)

    report = clash_report(coords, MetricConfig(clash_distance=2.0))

    assert report.n_atoms == 9


# ---------------------------------------------------------------------------
# clash_penalty
# ---------------------------------------------------------------------------


def test_clash_penalty_returns_one_value_per_structure_in_a_batch():
    """A batch of structures yields one penalty each."""
    torch.manual_seed(0)
    clashing = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=torch.float64)
    clean = torch.tensor([[0.0, 0.0, 0.0], [6.0, 0.0, 0.0]], dtype=torch.float64)

    values = clash_penalty(torch.stack([clashing, clean, clashing]))

    assert values.shape == (3,)
    assert float(values[0]) == pytest.approx(1.0, rel=1e-9)
    assert float(values[1]) == pytest.approx(0.0, abs=1e-12)
    assert float(values[2]) == pytest.approx(1.0, rel=1e-9)


def test_clash_penalty_of_a_clean_structure_is_zero():
    """A structure without any overlap scores exactly zero."""
    torch.manual_seed(0)
    coords = _points(10, seed=39, scale=20.0)

    value = clash_penalty(coords, MetricConfig(clash_distance=2.0))

    assert float(value) == pytest.approx(0.0, abs=1e-12)


def test_clash_penalty_agrees_with_clash_report():
    """The scalar penalty is the same number as the detailed report."""
    torch.manual_seed(0)
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.2, 0.0], [9.0, 0.0, 0.0]], dtype=torch.float64
    )
    cfg = MetricConfig(clash_distance=2.0)

    report = clash_report(coords, cfg)
    value = clash_penalty(coords, cfg)

    assert float(value) == pytest.approx(report.penalty, rel=1e-12)


def test_clash_penalty_batch_honours_residue_indices():
    """A 1-D residue index applies to every structure of the batch."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]], dtype=torch.float64)
    residue_index = torch.tensor([0, 0], dtype=torch.long)

    values = clash_penalty(
        torch.stack([coords, coords]), MetricConfig(clash_distance=2.0), residue_index=residue_index
    )

    assert values.shape == (2,)
    assert float(values[0]) == pytest.approx(0.0, abs=1e-12)
    assert float(values[1]) == pytest.approx(0.0, abs=1e-12)


def _three_part_complex():
    """A receptor, an interface patch and a ligand, as separate blocks."""
    generator = torch.Generator().manual_seed(3)
    receptor = torch.randn((20, 3), generator=generator, dtype=torch.float64) * 6.0
    interface = torch.randn((8, 3), generator=generator, dtype=torch.float64) * 1.2 + torch.tensor(
        [6.0, 0.0, 0.0], dtype=torch.float64
    )
    ligand = torch.randn((6, 3), generator=generator, dtype=torch.float64) * 0.8 + torch.tensor(
        [8.0, 0.0, 0.0], dtype=torch.float64
    )
    target = torch.cat([receptor, interface, ligand], dim=0)
    interface_mask = torch.zeros(34, dtype=torch.bool)
    interface_mask[20:28] = True
    receptor_mask = torch.zeros(34, dtype=torch.bool)
    receptor_mask[:20] = True
    return target, interface_mask, receptor_mask


def test_irmsd_of_a_rigid_motion_is_zero():
    """A rigid motion of the whole complex is undone exactly."""
    target, interface_mask, receptor_mask = _three_part_complex()
    angle = 0.7
    rotation = torch.tensor(
        [
            [math.cos(angle), -math.sin(angle), 0.0],
            [math.sin(angle), math.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ],
        dtype=torch.float64,
    )
    mobile = target @ rotation.T + torch.tensor([3.0, -2.0, 7.0], dtype=torch.float64)

    value = interface_rmsd(mobile, target, interface_mask, receptor_mask=receptor_mask)

    assert float(value) == pytest.approx(0.0, abs=1e-9)


def test_irmsd_equals_the_interface_displacement():
    """Displacing the interface by d gives an iRMSD of exactly d."""
    target, interface_mask, receptor_mask = _three_part_complex()
    displacement = 2.5
    mobile = target.clone()
    mobile[20:28] += torch.tensor([displacement, 0.0, 0.0], dtype=torch.float64)

    value = interface_rmsd(mobile, target, interface_mask, receptor_mask=receptor_mask)

    assert float(value) == pytest.approx(displacement, rel=1e-9)


def test_irmsd_is_not_contaminated_by_the_ligand():
    """The ligand's own atoms are kept out of the fit.

    Without an explicit receptor mask the fit uses everything that is not
    interface, which for a complex with a ligand silently includes the
    ligand - and a displaced ligand then tilts the fit.
    """
    target, interface_mask, receptor_mask = _three_part_complex()
    displacement = 2.5
    mobile = target.clone()
    mobile[28:] += torch.tensor([0.0, 0.0, displacement], dtype=torch.float64)

    with_ligand = interface_rmsd(mobile, target, interface_mask)
    without_ligand = interface_rmsd(
        mobile, target, interface_mask, receptor_mask=receptor_mask
    )

    assert float(without_ligand) == pytest.approx(0.0, abs=1e-9)
    assert float(with_ligand) > 0.5
