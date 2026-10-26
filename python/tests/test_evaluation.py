"""Tests for the sequence helpers and the one-call complex evaluation.

Two things are covered here:

* the residue-type <-> sequence mapping, which is what lets a structure
  predictor be handed a sequence instead of a list of integers, and
* :func:`pfe_design.metrics.evaluate_complex`, the façade that computes the
  four geometry measurements together.
"""

from __future__ import annotations

import pytest
import torch

from pfe_design import (
    AA1_ORDER,
    MetricConfig,
    canonical_sequence,
    residue_types_from_sequence,
    sequence_from_residue_types,
)
from pfe_design.metrics import evaluate_complex
from pfe_design.packing import ESMFoldBackend, PackingConfig, SideChainPackingHook

# ---------------------------------------------------------------------------
# sequences
# ---------------------------------------------------------------------------


def test_the_twenty_amino_acids_round_trip():
    """Every standard residue survives the trip through the index."""
    indices = residue_types_from_sequence(AA1_ORDER)
    assert len(indices) == 20
    assert sequence_from_residue_types(indices) == AA1_ORDER


def test_an_unknown_residue_is_written_as_unknown():
    """Index 0 is unknown, and is not guessed at."""
    assert sequence_from_residue_types([0, 1, 0]) == "XAX"
    assert residue_types_from_sequence("XAX") == [0, 1, 0]


def test_an_index_outside_the_known_range_is_written_as_unknown():
    """An index nobody defined must not crash or be interpreted."""
    assert sequence_from_residue_types([1, 99, 2]) == "AXR"


def test_three_letter_codes_are_accepted():
    """A structure file gives three-letter codes."""
    assert canonical_sequence("ALA GLY SER") == "AGS"


def test_modified_residues_become_the_residue_they_derive_from():
    """A predictor cannot be asked about MSE; it can be asked about MET."""
    assert canonical_sequence("MSE", format="three-letter") == "M"
    assert canonical_sequence("SEC", format="three-letter") == "C"
    assert canonical_sequence("PYL", format="three-letter") == "K"
    # Separated codes are unambiguous, so the default reads them the same way.
    assert canonical_sequence("MSE ALA") == "MA"


def test_ambiguous_residues_become_x():
    """ASX and GLX have no single-letter equivalent."""
    assert canonical_sequence("UNK ASX GLX") == "XXX"
    assert canonical_sequence("UNK", format="three-letter") == "X"


def test_whitespace_in_a_sequence_is_ignored():
    """A sequence written out in columns still reads correctly."""
    assert canonical_sequence("arn dC q") == "ARNDCQ"


def test_a_lower_case_sequence_is_read():
    """Case is not information."""
    assert canonical_sequence("arndcq") == "ARNDCQ"


# ---------------------------------------------------------------------------
# the packing request carries a sequence
# ---------------------------------------------------------------------------


def test_the_packing_request_carries_a_sequence():
    """The request is what a structure predictor is handed."""
    hook = SideChainPackingHook(PackingConfig())
    coords = torch.randn(5, 3)

    request = hook.prepare(coords, residue_types_from_sequence("ACDEF"))

    assert request.sequence == "ACDEF"


def test_an_explicit_sequence_is_kept():
    """A caller may know a better sequence than the residue types imply."""
    hook = SideChainPackingHook(PackingConfig())
    coords = torch.randn(5, 3)

    request = hook.prepare(coords, residue_types_from_sequence("ACDEF"), sequence="ACDEW")

    assert request.sequence == "ACDEW"


def test_the_esmfold_backend_receives_a_sequence():
    """The seam can actually be driven: the predictor gets the sequence."""
    seen = {}

    def predictor(request):
        seen["sequence"] = request.sequence
        return request.backbone_xyz

    backend = ESMFoldBackend(predictor)
    hook = SideChainPackingHook(PackingConfig(backend="esmfold"), backend=backend)
    coords = torch.randn(5, 3)
    request = hook.prepare(coords, residue_types_from_sequence("ACDEF"))

    backend.pack(request)

    assert seen["sequence"] == "ACDEF"


# ---------------------------------------------------------------------------
# evaluate_complex
# ---------------------------------------------------------------------------


def _complex():
    """A receptor, an interface patch and a ligand, built on a lattice.

    The points are on a 5 Angstrom lattice so that no two of them are within
    a clash threshold: a random cloud would have accidental clashes, and the
    test would then be measuring the wrong thing.
    """
    receptor = torch.tensor(
        [[0.0, y, z] for y in (0.0, 5.0, 10.0, 15.0) for z in (0.0, 5.0, 10.0, 15.0, 20.0)],
        dtype=torch.float64,
    )
    interface = torch.tensor(
        [[-4.0, y, z] for y in (2.5, 7.5, 12.5) for z in (2.5, 7.5)],
        dtype=torch.float64,
    )[:6]
    ligand = torch.tensor(
        [[-9.5, y, z] for y in (5.0, 10.0) for z in (5.0, 10.0, 15.0)], dtype=torch.float64
    )
    target = torch.cat([receptor, interface, ligand], dim=0)

    interface_mask = torch.zeros(32, dtype=torch.bool)
    interface_mask[20:26] = True
    receptor_mask = torch.zeros(32, dtype=torch.bool)
    receptor_mask[:20] = True
    ligand_mask = torch.zeros(32, dtype=torch.bool)
    ligand_mask[26:] = True
    group_a = torch.arange(20)
    group_b = torch.arange(20, 32)
    return target, interface_mask, receptor_mask, ligand_mask, group_a, group_b


def test_evaluate_complex_computes_all_four_measurements():
    """One call gives the interface RMSD, the ligand RMSD, the area and the clashes."""
    target, interface, receptor, ligand, group_a, group_b = _complex()
    radii = torch.full((32,), 1.70, dtype=torch.float64)

    result = evaluate_complex(
        target,
        target,
        interface,
        receptor_mask=receptor,
        ligand_mask=ligand,
        group_a=group_a,
        group_b=group_b,
        radii=radii,
        cfg=MetricConfig(),
    )

    assert result.irmsd == pytest.approx(0.0, abs=1e-9)
    assert result.ligand_rmsd == pytest.approx(0.0, abs=1e-9)
    assert result.bsa is not None and result.bsa > 0.0
    assert result.n_clashes == 0
    assert result.is_clean


def test_evaluate_complex_finds_the_clashes():
    """A pair that is too close is reported."""
    target, interface, receptor, ligand, group_a, group_b = _complex()
    moved = target.clone()
    moved[1] = moved[0] + torch.tensor([0.2, 0.0, 0.0], dtype=torch.float64)

    result = evaluate_complex(moved, target, interface, receptor_mask=receptor)

    assert result.n_clashes == 1
    assert result.clash_penalty > 0.0
    assert result.max_overlap > 0.0
    assert not result.is_clean


def test_evaluate_complex_reports_what_it_was_not_given():
    """A measurement that needs more information is None, not a guess."""
    target, interface, receptor, ligand, group_a, group_b = _complex()

    result = evaluate_complex(target, target, interface)

    assert result.irmsd is not None
    assert result.ligand_rmsd is None
    assert result.bsa is None


def test_evaluate_complex_refuses_a_batch():
    """The batched metrics have their own entry points."""
    target, interface, receptor, ligand, group_a, group_b = _complex()

    with pytest.raises(ValueError, match="one complex"):
        evaluate_complex(target.unsqueeze(0), target.unsqueeze(0), interface.unsqueeze(0))


def test_a_one_letter_sequence_is_not_read_as_three_letter_codes():
    """Every three-letter code is also a run of three one-letter codes.

    'MET' is one residue written in three-letter form and three residues
    written in one-letter form, and reading it the wrong way changes the
    length of the sequence - which is the sequence a predictor is given.
    """
    assert canonical_sequence("MET") == "MET"
    assert len(residue_types_from_sequence("MET")) == 3
    assert canonical_sequence("METAL") == "METAL"
    assert canonical_sequence("SERINE") == "SERINE"
    assert canonical_sequence("GLYCINE") == "GLYCINE"


def test_a_three_letter_sequence_is_read_as_three_letter_codes():
    """Separated codes are unambiguous, so they are read that way."""
    assert canonical_sequence("MSE ALA SEC GLY") == "MACG"
    assert canonical_sequence("ALA GLY SER") == "AGS"
    assert canonical_sequence("MSE,ALA") == "MA"


def test_the_format_can_be_stated_explicitly():
    """A caller who knows the format does not have to rely on the guess."""
    assert canonical_sequence("ALAGLYSER", format="three-letter") == "AGS"
    assert canonical_sequence("ALAGLYSER", format="one-letter") == "ALAGLYSER"
    assert canonical_sequence("MSE ALA", format="one-letter") == "MSEALA"


def test_an_unknown_format_is_refused():
    """A typo in the format must not be silently ignored."""
    with pytest.raises(ValueError, match="unknown sequence format"):
        residue_types_from_sequence("ACDEF", format="guess")


def test_two_overlapping_partners_are_refused():
    """A group cannot bury itself.

    With overlapping groups the union collapses to one of them, and the
    burial comes out as that group's own surface - a large positive number
    for a pair of molecules that are really one.
    """
    generator = torch.Generator().manual_seed(2)
    coords = torch.randn((8, 3), generator=generator, dtype=torch.float64) * 3.0
    radii = torch.full((8,), 1.70, dtype=torch.float64)

    with pytest.raises(ValueError, match="cannot bury itself"):
        from pfe_design.metrics import buried_surface_area

        buried_surface_area(coords, radii, torch.arange(5), torch.arange(3, 8))


def test_an_empty_interface_is_reported_as_not_measurable():
    """The facade reports None rather than a perfect zero."""
    target, interface, receptor, ligand, group_a, group_b = _complex()
    empty = torch.zeros_like(interface)

    result = evaluate_complex(target, target, empty)

    assert result.irmsd is None


def test_the_ligand_is_never_part_of_its_own_fit():
    """Without an explicit receptor, the ligand is excluded from the fit.

    The ligand is measured by the receptor's fit, so letting the ligand
    into the fit would let it follow itself and report a perfect score.
    """
    target, interface, receptor, ligand, group_a, group_b = _complex()
    moved = target.clone()
    moved[26:] += torch.tensor([0.0, 3.0, 0.0], dtype=torch.float64)   # the ligand moves

    implicit = evaluate_complex(moved, target, interface, ligand_mask=ligand)
    explicit = evaluate_complex(
        moved, target, interface, receptor_mask=receptor, ligand_mask=ligand
    )

    assert implicit.ligand_rmsd == pytest.approx(3.0, rel=1e-6)
    assert explicit.ligand_rmsd == pytest.approx(3.0, rel=1e-6)


def test_boolean_group_masks_are_read_as_atom_lists():
    """A boolean mask is a list of atoms, not a list of indices."""
    generator = torch.Generator().manual_seed(3)
    coords = torch.randn((6, 3), generator=generator, dtype=torch.float64) * 3.0
    radii = torch.full((6,), 1.70, dtype=torch.float64)
    group_a = torch.tensor([True, True, False, False, False, False])
    group_b = torch.tensor([False, False, True, False, False, False])

    from pfe_design.metrics import buried_surface_area

    result = buried_surface_area(coords, radii, group_a, group_b)

    assert result.bsa >= 0.0


def test_overlapping_boolean_masks_are_still_refused():
    """The overlap check works on masks as well as on index lists."""
    generator = torch.Generator().manual_seed(3)
    coords = torch.randn((6, 3), generator=generator, dtype=torch.float64) * 3.0
    radii = torch.full((6,), 1.70, dtype=torch.float64)
    group_a = torch.tensor([True, True, False, False, False, False])
    overlapping = torch.tensor([True, False, True, False, False, False])

    from pfe_design.metrics import buried_surface_area

    with pytest.raises(ValueError, match="cannot bury itself"):
        buried_surface_area(coords, radii, group_a, overlapping)


def test_an_unknown_implementation_is_refused():
    """A typo must not fall back silently to the other implementation."""
    generator = torch.Generator().manual_seed(3)
    coords = torch.randn((6, 3), generator=generator, dtype=torch.float64) * 3.0
    radii = torch.full((6,), 1.70, dtype=torch.float64)

    from pfe_design.metrics import buried_surface_area

    with pytest.raises(ValueError, match="unknown implementation"):
        buried_surface_area(
            coords,
            radii,
            torch.arange(3),
            torch.arange(3, 6),
            implementation="biopythn",
        )


def test_the_ligand_does_not_tilt_the_interface_fit_either():
    """The interface RMSD uses the same ligand-free fit as the ligand RMSD."""
    generator = torch.Generator().manual_seed(9)
    target = torch.randn((34, 3), generator=generator, dtype=torch.float64) * 5.0
    interface = torch.zeros(34, dtype=torch.bool)
    interface[20:28] = True
    ligand = torch.zeros(34, dtype=torch.bool)
    ligand[28:] = True
    moved = target.clone()
    moved[28:] += torch.tensor([0.0, 0.0, 5.0], dtype=torch.float64)

    result = evaluate_complex(moved, target, interface, ligand_mask=ligand)

    assert result.irmsd == pytest.approx(0.0, abs=1e-6)
    assert result.ligand_rmsd == pytest.approx(5.0, rel=1e-6)
