"""Cross-check the surface-area calculation against Biopython.

`pfe_design.metrics.sasa` is a vectorised Shrake-Rupley implementation in
PyTorch. Biopython ships an independent implementation of the same method
(`Bio.PDB.SASA.ShrakeRupley`), so the two can be compared directly: if they
agree on the same coordinates, the PyTorch version is not merely plausible, it
agrees with the reference implementation.

Two caveats are worth stating, because they shape the assertions:

* The two implementations distribute their sphere points differently
  (Biopython includes the poles, this one uses midpoint sampling), so an
  individual atom can differ by more than the total does. For a heavily
  buried atom only a handful of points are exposed, so the relative error of
  that atom's area is large while the error of the total stays small.
* The error therefore shrinks as the sampling gets finer. Asserting that it
  shrinks is a stronger statement than asserting a fixed tolerance, because
  it distinguishes sampling noise from a systematic error.

The whole module is skipped when Biopython is not installed, which keeps the
package's test-suite runnable with only its declared dependencies.
"""

from __future__ import annotations

import math

import pytest

torch = pytest.importorskip("torch")

sasa_module = pytest.importorskip("Bio.PDB.SASA")

from Bio.PDB.Atom import Atom  # noqa: E402
from Bio.PDB.Chain import Chain  # noqa: E402
from Bio.PDB.Model import Model  # noqa: E402
from Bio.PDB.Residue import Residue  # noqa: E402
from Bio.PDB.Structure import Structure  # noqa: E402

from pfe_design.config import MetricConfig  # noqa: E402
from pfe_design.metrics import (  # noqa: E402
    buried_surface_area,
    sasa,
    sasa_biopython,
)


def _biopython_sasa(coords, elements, n_points=100):
    """Per-atom SASA from Biopython's Shrake-Rupley implementation."""
    structure = Structure("test")
    model = Model(0)
    structure.add(model)
    chain = Chain("A")
    model.add(chain)
    for index, (xyz, element) in enumerate(zip(coords, elements, strict=False)):
        residue = Residue((" ", index + 1, " "), element, "")
        atom = Atom(
            name=element,
            coord=xyz.tolist(),
            bfactor=0.0,
            occupancy=1.0,
            altloc=" ",
            fullname=element,
            serial_number=index + 1,
            element=element,
        )
        residue.add(atom)
        chain.add(residue)
    calculator = sasa_module.ShrakeRupley(n_points=n_points)
    calculator.compute(structure, level="R")
    values = []
    for residue in structure.get_residues():
        for atom in residue:
            values.append(float(atom.sasa))
    return values


def _cluster(seed=1, n_atoms=12, spread=4.0):
    generator = torch.Generator().manual_seed(seed)
    coords = torch.randn((n_atoms, 3), generator=generator, dtype=torch.float64) * spread
    elements = ["C", "N", "O", "S"] * (n_atoms // 4 + 1)
    return coords, elements[:n_atoms]


def _radii(elements):
    return torch.tensor(
        [sasa_module.ATOMIC_RADII.get(e, 1.70) for e in elements], dtype=torch.float64
    )


def test_sasa_of_an_isolated_atom_matches_the_analytic_area():
    """A single atom has the analytic surface area in both implementations."""
    coords = torch.tensor([[0.0, 0.0, 0.0]], dtype=torch.float64)
    radii = torch.tensor([1.70], dtype=torch.float64)
    analytic = 4.0 * math.pi * (1.70 + 1.4) ** 2

    mine = float(sasa(coords, radii, probe=1.4, n_points=200)[0])
    theirs = _biopython_sasa(coords, ["C"])[0]

    assert mine == pytest.approx(analytic, rel=0.01)
    assert theirs == pytest.approx(analytic, rel=0.01)
    assert mine == pytest.approx(theirs, rel=0.02)


def test_total_sasa_of_a_cluster_matches_biopython():
    """The total accessible area of a cluster agrees to within a per cent."""
    coords, elements = _cluster()
    radii = _radii(elements)

    mine = float(sasa(coords, radii, probe=1.4, n_points=100).sum())
    theirs = sum(_biopython_sasa(coords, elements, n_points=100))

    assert mine == pytest.approx(theirs, rel=0.01)


def test_disagreement_shrinks_as_the_sampling_gets_finer():
    """The residual difference is sampling noise, not a systematic error.

    If the two implementations were doing something different - a different
    probe radius, say - the error would not shrink with more points. The
    error is measured as the RMS per-atom difference relative to the total
    area, because the relative error of a single heavily buried atom is
    dominated by its small denominator and stays noisy at any resolution.
    """
    coords, elements = _cluster()
    radii = _radii(elements)

    errors = []
    for n_points in (100, 4000):
        mine = sasa(coords, radii, probe=1.4, n_points=n_points).tolist()
        theirs = _biopython_sasa(coords, elements, n_points=n_points)
        total = sum(theirs)
        rms = math.sqrt(sum((a - b) ** 2 for a, b in zip(mine, theirs, strict=False)) / len(mine))
        errors.append(rms / total)

    assert errors[1] < errors[0] / 4, f"the error did not shrink: {errors}"
    assert errors[1] < 0.001, f"the residual error is still large: {errors[1]}"


def test_buried_surface_area_matches_the_reference_definition():
    """BSA computed with the PyTorch implementation matches Biopython."""
    coords, elements = _cluster(seed=5, n_atoms=10)
    radii = _radii(elements)
    group_a = torch.arange(5)
    group_b = torch.arange(5, 10)

    result = buried_surface_area(
        coords, radii, group_a, group_b, MetricConfig(bsa_n_points=100)
    )

    free_a = sum(_biopython_sasa(coords[group_a], [elements[i] for i in group_a.tolist()]))
    free_b = sum(_biopython_sasa(coords[group_b], [elements[i] for i in group_b.tolist()]))
    complex_sasa = sum(_biopython_sasa(coords, elements))

    # BSA is a small difference between large numbers, so the sampling error
    # of the three areas is amplified by roughly an order of magnitude in the
    # difference.
    assert result.bsa == pytest.approx(free_a + free_b - complex_sasa, rel=0.05)


def test_buried_surface_area_is_zero_for_two_far_apart_groups():
    """Groups that never touch bury nothing."""
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [0.0, 0.0, 1.5], [50.0, 0.0, 0.0], [50.0, 0.0, 1.5]],
        dtype=torch.float64,
    )
    radii = torch.tensor([1.70, 1.70, 1.70, 1.70], dtype=torch.float64)
    result = buried_surface_area(
        coords, radii, torch.tensor([0, 1]), torch.tensor([2, 3]), MetricConfig()
    )

    assert result.bsa == pytest.approx(0.0, abs=0.5)


def test_sasa_agrees_with_biopython_for_a_system_larger_than_one_chunk():
    """The cross-check at a size that exercises the chunking path.

    The calculation is processed in chunks of atoms, and every earlier
    cross-check used a system small enough to fit in a single chunk. This
    one does not, so it covers the code path where the self-occlusion mask
    has to be built on the atom axis rather than the chunk axis.
    """
    coords, elements = _cluster(seed=9, n_atoms=60, spread=5.0)
    radii = _radii(elements)

    mine = float(sasa(coords, radii, probe=1.4, n_points=100).sum())
    theirs = sum(_biopython_sasa(coords, elements, n_points=100))

    assert mine == pytest.approx(theirs, rel=0.02)


def test_the_module_itself_can_use_biopython():
    """The reference implementation is available from the metrics module.

    The PyTorch calculation is the fast one; this is the reference, and a
    user who wants the reference number rather than the fast one should be
    able to have it without leaving the package.
    """
    coords, elements = _cluster(seed=4, n_atoms=24)
    radii = _radii(elements)

    mine = sasa_biopython(coords, elements, probe=1.4, n_points=100)
    theirs = sasa(coords, radii, probe=1.4, n_points=100)

    assert mine.shape == (24,)
    assert float(mine.sum()) == pytest.approx(float(theirs.sum()), rel=0.02)


def test_buried_surface_area_can_be_computed_by_biopython():
    """Both implementations of the burial agree to within the sampling noise.

    Over a sample of 48 systems the two differ by 0.1 % to 13 % at 100
    points per atom, with a median near 1 %: a buried area is a difference
    between two large surface areas, so the sampling error is amplified.
    The bound below is therefore generous - what matters is the next test.
    """
    coords, elements = _cluster(seed=8, n_atoms=20)
    radii = _radii(elements)
    group_a = torch.arange(10)
    group_b = torch.arange(10, 20)

    fast = buried_surface_area(coords, radii, group_a, group_b, MetricConfig())
    reference = buried_surface_area(
        coords,
        radii,
        group_a,
        group_b,
        MetricConfig(),
        implementation="biopython",
        elements=elements,
    )

    assert reference.bsa == pytest.approx(fast.bsa, rel=0.2)


def test_the_burial_agreement_improves_as_the_sampling_refines():
    """A systematic difference would not shrink; sampling noise does.

    This is the assertion that tells the two implementations apart from a
    pair that merely look similar: the difference between them falls as the
    number of sphere points rises.
    """
    coords, elements = _cluster(seed=8, n_atoms=20)
    radii = _radii(elements)
    group_a = torch.arange(10)
    group_b = torch.arange(10, 20)

    differences = []
    for n_points in (100, 400):
        fast = buried_surface_area(
            coords, radii, group_a, group_b, MetricConfig(bsa_n_points=n_points)
        )
        reference = buried_surface_area(
            coords,
            radii,
            group_a,
            group_b,
            MetricConfig(bsa_n_points=n_points),
            implementation="biopython",
            elements=elements,
        )
        differences.append(abs(fast.bsa - reference.bsa) / abs(reference.bsa))

    assert differences[1] < differences[0], f"the difference did not shrink: {differences}"


def test_the_biopython_path_says_what_it_needs():
    """Without the elements there is nothing to resolve a radius from."""
    coords, elements = _cluster(seed=4, n_atoms=12)
    radii = _radii(elements)

    with pytest.raises(ValueError, match="elements"):
        buried_surface_area(
            coords,
            radii,
            torch.arange(6),
            torch.arange(6, 12),
            MetricConfig(),
            implementation="biopython",
        )


def test_the_two_implementations_agree_on_an_empty_system():
    """An empty structure has no area, in either implementation."""
    empty = torch.zeros((0, 3), dtype=torch.float64)

    assert sasa_biopython(empty, [], 1.4, 50).shape == (0,)
    assert sasa(empty, torch.zeros(0, dtype=torch.float64), 1.4, 50).shape == (0,)
