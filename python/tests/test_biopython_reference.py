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
from pfe_design.metrics import buried_surface_area, sasa  # noqa: E402


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
