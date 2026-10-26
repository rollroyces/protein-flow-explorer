"""Tests for :mod:`pfe_design.conditioning`.

The conditioner turns a target structure plus a binary complex mask into the
per-residue context consumed by the vector field.  Two contracts matter most:

* the **mask convention** - the first ``n_target`` columns of the mask are the
  fixed target, the trailing columns are the designed (dynamic) residues, and
* the **pocket resolution order** - an explicit centre wins, otherwise the
  design region defines the site, otherwise the whole target is the pocket.
"""

from __future__ import annotations

import pytest
import torch

from pfe_design.conditioning import (
    ConditioningOutput,
    PocketConditioner,
    build_complex_mask,
    dynamic_mask_from_layout,
)
from pfe_design.config import ModelConfig, TargetConditioningConfig

N_TARGET = 4
N_DESIGN = 3
FEATURE_DIM = 16

# A target whose first two residues form a tight cluster and whose last two
# residues are far away from it.
TARGET_COORDS = torch.tensor(
    [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [50.0, 0.0, 0.0],
        [100.0, 0.0, 0.0],
    ],
    dtype=torch.float32,
).unsqueeze(0)

# The designed region sits on top of the third target residue (index 2).
DESIGN_NEAR_RESIDUE_2 = torch.tensor(
    [
        [50.0, 0.0, 0.0],
        [50.5, 0.0, 0.0],
        [51.0, 0.0, 0.0],
    ],
    dtype=torch.float32,
).unsqueeze(0)


def _conditioner(**cfg_kwargs) -> PocketConditioner:
    """A small, deterministic conditioner."""
    torch.manual_seed(0)
    cfg = TargetConditioningConfig(feature_dim=FEATURE_DIM, **cfg_kwargs)
    model_cfg = ModelConfig(hidden_dim=32, n_heads=4, n_layers=2)
    return PocketConditioner(cfg, model_cfg)


def _mask(n_target: int = N_TARGET, n_design: int = N_DESIGN) -> torch.Tensor:
    return build_complex_mask(n_target, n_design, batch_size=1)


def _valid(n_target: int = N_TARGET) -> torch.Tensor:
    return torch.ones(1, n_target, dtype=torch.float32)


def _complex(target: torch.Tensor, design: torch.Tensor) -> torch.Tensor:
    return torch.cat([target, design], dim=1)


# ------------------------------------------------------------- build_complex_mask


def test_build_complex_mask_has_expected_shape():
    """The mask has one column per residue of the complex."""
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=1)
    assert mask.shape == (1, N_TARGET + N_DESIGN)


def test_build_complex_mask_marks_design_columns_dynamic():
    """The trailing ``n_design`` columns are the dynamic (designed) residues."""
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=1)
    assert torch.equal(mask[0, N_TARGET:], torch.ones(N_DESIGN))
    assert torch.equal(mask[0, :N_TARGET], torch.zeros(N_TARGET))


def test_build_complex_mask_target_columns_are_zero():
    """The target columns are exactly zero - the target is fixed."""
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=1)
    assert torch.count_nonzero(mask[0, :N_TARGET]).item() == 0


def test_build_complex_mask_values_are_binary():
    """The mask only ever contains 0.0 and 1.0."""
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=2)
    unique = torch.unique(mask)
    assert set(unique.tolist()) == {0.0, 1.0}


def test_build_complex_mask_batch_rows_are_identical():
    """Every complex in the batch shares the same layout."""
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=4)
    assert mask.shape[0] == 4
    for row in range(1, 4):
        assert torch.equal(mask[0], mask[row])


def test_build_complex_mask_without_design_region_is_all_zero():
    """A target-only complex has no dynamic residues at all."""
    mask = build_complex_mask(N_TARGET, 0, batch_size=1)
    assert mask.shape == (1, N_TARGET)
    assert torch.count_nonzero(mask).item() == 0


def test_build_complex_mask_rejects_negative_counts():
    """Negative residue counts are a programming error."""
    with pytest.raises(ValueError):
        build_complex_mask(-1, N_DESIGN)


def test_dynamic_mask_from_layout_matches_build_complex_mask():
    """The readability alias returns exactly the same mask."""
    a = dynamic_mask_from_layout(N_TARGET, N_DESIGN, batch_size=2)
    b = build_complex_mask(N_TARGET, N_DESIGN, batch_size=2)
    assert torch.equal(a, b)


# ----------------------------------------------------------------- input validation


def test_conditioner_rejects_dynamic_target_columns():
    """The target must be fixed, so its mask columns have to be zero."""
    conditioner = _conditioner()
    mask = _mask()
    mask[0, 0] = 1.0
    with pytest.raises(ValueError):
        conditioner(TARGET_COORDS, _valid(), mask)


def test_conditioner_rejects_dynamic_target_column_in_any_batch_row():
    """The rule holds for every row of the batch, not only the first."""
    conditioner = _conditioner()
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=2)
    mask[1, 2] = 1.0
    with pytest.raises(ValueError):
        conditioner(
            TARGET_COORDS.expand(2, -1, -1).contiguous(),
            _valid().expand(2, -1).contiguous(),
            mask,
        )


def test_conditioner_rejects_mask_narrower_than_target():
    """The complex mask must have at least as many columns as the target."""
    conditioner = _conditioner()
    narrow = build_complex_mask(N_TARGET - 1, N_DESIGN, batch_size=1)
    with pytest.raises(ValueError):
        conditioner(TARGET_COORDS, _valid(), narrow)


def test_conditioner_rejects_wrong_target_valid_shape():
    """``target_valid`` must have exactly one entry per target residue."""
    conditioner = _conditioner()
    with pytest.raises(ValueError):
        conditioner(TARGET_COORDS, torch.ones(1, N_TARGET + 1), _mask())


def test_conditioner_rejects_wrong_target_coordinate_shape():
    """The target coordinates must be ``(batch, n_target, 3)``."""
    conditioner = _conditioner()
    with pytest.raises(ValueError):
        conditioner(TARGET_COORDS[0], _valid(), _mask())


def test_conditioner_rejects_wrong_batch_size_in_mask():
    """The batch dimension of the mask must match the target."""
    conditioner = _conditioner()
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=2)
    with pytest.raises(ValueError):
        conditioner(TARGET_COORDS, _valid(), mask)


# ------------------------------------------------------------------ pocket resolution


def test_pocket_uses_explicit_centre_as_given():
    """An explicit pocket centre is used verbatim."""
    conditioner = _conditioner()
    centre = torch.tensor([[0.0, 0.0, 0.0]])
    cond = conditioner(TARGET_COORDS, _valid(), _mask(), pocket_center=centre)
    assert torch.equal(cond.pocket_center, centre)


def test_pocket_membership_with_explicit_centre_is_nearby_target_only():
    """Only the target residues within the cutoff line the explicit pocket."""
    conditioner = _conditioner()
    centre = torch.tensor([[0.0, 0.0, 0.0]])
    cond = conditioner(TARGET_COORDS, _valid(), _mask(), pocket_center=centre)
    expected = torch.tensor([[1.0, 1.0, 0.0, 0.0]])
    assert torch.equal(cond.pocket_mask, expected)


def test_far_target_residue_is_not_a_pocket_member():
    """A target residue far from the explicit centre does not line the pocket."""
    conditioner = _conditioner()
    centre = torch.tensor([[0.0, 0.0, 0.0]])
    cond = conditioner(TARGET_COORDS, _valid(), _mask(), pocket_center=centre)
    assert cond.pocket_mask[0, 2].item() == 0.0
    assert cond.pocket_mask[0, 3].item() == 0.0
    assert cond.pocket_mask[0, 0].item() == 1.0


def test_pocket_is_derived_from_design_region_proximity():
    """Without an explicit centre the design region defines the binding site."""
    conditioner = _conditioner()
    complex_xyz = _complex(TARGET_COORDS, DESIGN_NEAR_RESIDUE_2)
    cond = conditioner(TARGET_COORDS, _valid(), _mask(), complex_xyz=complex_xyz)
    expected = torch.tensor([[0.0, 0.0, 1.0, 0.0]])
    assert torch.equal(cond.pocket_mask, expected)


def test_pocket_derived_from_design_excludes_distant_target_residues():
    """Target residues away from the designed region do not line the pocket."""
    conditioner = _conditioner()
    complex_xyz = _complex(TARGET_COORDS, DESIGN_NEAR_RESIDUE_2)
    cond = conditioner(TARGET_COORDS, _valid(), _mask(), complex_xyz=complex_xyz)
    assert cond.pocket_mask[0, 0].item() == 0.0
    assert cond.pocket_mask[0, 1].item() == 0.0
    assert cond.pocket_mask[0, 3].item() == 0.0


def test_pocket_is_whole_target_when_neither_centre_nor_design_given():
    """With no centre and no design coordinates the pocket is the whole target."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert torch.equal(cond.pocket_mask, _valid())


def test_pocket_centre_falls_back_to_target_centroid():
    """Without any hint the pocket centre is the centroid of the target."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    expected = TARGET_COORDS[0].mean(dim=0)
    assert torch.allclose(cond.pocket_center[0], expected, atol=1e-5)


def test_pocket_excludes_unmodelled_target_residues():
    """An unmodelled target residue never lines the pocket."""
    conditioner = _conditioner()
    valid = _valid()
    valid[0, 1] = 0.0
    cond = conditioner(TARGET_COORDS, valid, _mask(), pocket_center=torch.zeros(1, 3))
    assert cond.pocket_mask[0, 1].item() == 0.0
    assert cond.pocket_mask[0, 0].item() == 1.0


def test_pocket_is_capped_to_the_configured_number_of_residues():
    """The pocket keeps the configured number of closest residues."""
    conditioner = _conditioner(max_pocket_residues=2)
    cond = conditioner(
        TARGET_COORDS, _valid(), _mask(), pocket_center=torch.zeros(1, 3)
    )
    # All four target residues line the pocket, but only the two closest are
    # kept (distances 0 and 2 from the centre).
    assert cond.pocket_mask[0, 0].item() == 1.0
    assert cond.pocket_mask[0, 1].item() == 1.0
    assert cond.pocket_mask[0, 2].item() == 0.0
    assert cond.pocket_mask[0, 3].item() == 0.0
    assert cond.pocket_mask.sum().item() == 2.0


def test_pocket_membership_flag_can_be_switched_off():
    """Without the membership flag the context no longer knows the pocket."""
    with_flag = _conditioner(use_pocket_membership=True)(
        TARGET_COORDS, _valid(), _mask(), pocket_center=torch.zeros(1, 3)
    ).context
    without_flag = _conditioner(use_pocket_membership=False)(
        TARGET_COORDS, _valid(), _mask(), pocket_center=torch.zeros(1, 3)
    ).context
    assert not torch.allclose(with_flag, without_flag, atol=1e-6)


# ------------------------------------------------------------- output shape invariants


def test_conditioning_output_context_has_expected_shape():
    """The context is ``(batch, n_complex, feature_dim)``."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert cond.context.shape == (1, N_TARGET + N_DESIGN, FEATURE_DIM)


def test_conditioning_output_is_a_conditioning_output_instance():
    """The conditioner returns a ``ConditioningOutput``."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert isinstance(cond, ConditioningOutput)


def test_dynamic_and_fixed_masks_partition_the_complex():
    """``dynamic_mask + fixed_mask`` is one everywhere, with no overlap."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    total = cond.dynamic_mask + cond.fixed_mask
    assert torch.equal(total, torch.ones_like(total))


def test_fixed_mask_covers_exactly_the_target_columns():
    """The fixed residues are the target ones."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert torch.equal(cond.fixed_mask[0, :N_TARGET], torch.ones(N_TARGET))
    assert torch.equal(cond.fixed_mask[0, N_TARGET:], torch.zeros(N_DESIGN))


def test_fixed_xyz_holds_target_coordinates_in_first_columns():
    """The fixed coordinates are the target ones, in the leading columns."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert torch.equal(cond.fixed_xyz[0, :N_TARGET], TARGET_COORDS[0])


def test_fixed_xyz_is_zero_on_the_designed_residues():
    """The designed rows of ``fixed_xyz`` are zeroed out."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert torch.count_nonzero(cond.fixed_xyz[0, N_TARGET:]).item() == 0


def test_pocket_radius_is_positive():
    """The pocket radius is a strictly positive length."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert cond.pocket_radius.shape == (1,)
    assert bool((cond.pocket_radius > 0).all())


def test_pocket_radius_spans_the_pocket_members():
    """The radius reaches the furthest residue that lines the pocket."""
    conditioner = _conditioner()
    cond = conditioner(
        TARGET_COORDS, _valid(), _mask(), pocket_center=torch.zeros(1, 3)
    )
    assert cond.pocket_radius[0].item() == pytest.approx(2.0, abs=1e-5)


def test_valid_mask_marks_unmodelled_target_residues_invalid():
    """An unmodelled target residue is flagged invalid in the output."""
    conditioner = _conditioner()
    valid = _valid()
    valid[0, 2] = 0.0
    cond = conditioner(TARGET_COORDS, valid, _mask())
    assert cond.valid_mask[0, 2].item() == 0.0
    assert cond.valid_mask[0, 0].item() == 1.0


def test_valid_mask_marks_the_designed_residues_valid():
    """Every designed residue is a modelled residue."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert torch.equal(cond.valid_mask[0, N_TARGET:], torch.ones(N_DESIGN))


def test_residue_counts_are_consistent():
    """The reported counts match the complex layout."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert cond.n_target == N_TARGET
    assert cond.n_design == N_DESIGN
    assert cond.n_residues == N_TARGET + N_DESIGN


# ----------------------------------------------------------------- context content


def test_context_depends_on_the_target_structure():
    """Two different target structures give two different contexts."""
    conditioner = _conditioner()
    clustered = torch.tensor(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0], [6.0, 0.0, 0.0]]
    ).unsqueeze(0)
    spread = torch.tensor(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0], [60.0, 0.0, 0.0]]
    ).unsqueeze(0)
    centre = torch.zeros(1, 3)
    ctx_a = conditioner(clustered, _valid(), _mask(), pocket_center=centre).context
    ctx_b = conditioner(spread, _valid(), _mask(), pocket_center=centre).context
    assert not torch.allclose(ctx_a, ctx_b, atol=1e-5)
    assert (ctx_a - ctx_b).abs().max().item() > 1e-3


def test_context_changes_when_the_pocket_centre_moves():
    """Moving the pocket centre changes the geometric features."""
    conditioner = _conditioner()
    ctx_a = conditioner(
        TARGET_COORDS, _valid(), _mask(), pocket_center=torch.tensor([[0.0, 0.0, 0.0]])
    ).context
    ctx_b = conditioner(
        TARGET_COORDS, _valid(), _mask(), pocket_center=torch.tensor([[3.0, 0.0, 0.0]])
    ).context
    assert not torch.allclose(ctx_a, ctx_b, atol=1e-5)
    assert (ctx_a - ctx_b).abs().max().item() > 1e-3


def test_context_is_reproducible_for_the_same_input():
    """The conditioner is deterministic given the same input."""
    conditioner = _conditioner()
    first = conditioner(TARGET_COORDS, _valid(), _mask()).context
    second = conditioner(TARGET_COORDS, _valid(), _mask()).context
    assert torch.equal(first, second)


def test_context_is_finite():
    """No NaN or inf leaks into the context, even for a distant target."""
    conditioner = _conditioner()
    cond = conditioner(TARGET_COORDS, _valid(), _mask())
    assert bool(torch.isfinite(cond.context).all())


def test_conditioning_handles_a_batch_of_two_complexes():
    """Batching does not change the result of either complex."""
    conditioner = _conditioner()
    # The second target is a genuinely different arrangement, not a copy.
    target_b = TARGET_COORDS.clone()
    target_b[0, 1, 0] = 20.0
    target = torch.cat([TARGET_COORDS, target_b], dim=0)
    valid = torch.ones(2, N_TARGET)
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=2)

    batched = conditioner(target, valid, mask)
    assert batched.context.shape == (2, N_TARGET + N_DESIGN, FEATURE_DIM)
    assert batched.pocket_mask.shape == (2, N_TARGET)

    # Each row of the batch matches the result of conditioning that row alone:
    # no information leaks between the complexes of a batch.
    solo_a = conditioner(TARGET_COORDS, _valid(), _mask())
    solo_b = conditioner(target_b, _valid(), _mask())
    assert torch.allclose(batched.context[0], solo_a.context[0], atol=1e-6)
    assert torch.allclose(batched.context[1], solo_b.context[0], atol=1e-6)
    # The two complexes really are different from each other.
    assert not torch.allclose(batched.context[0], batched.context[1], atol=1e-6)


# ---------------------------------------------------------------------- gradients


def test_detach_target_blocks_gradient_into_target_coordinates():
    """With ``detach_target`` the target receives no gradient."""
    conditioner = _conditioner(detach_target=True)
    target = TARGET_COORDS.clone().requires_grad_(True)
    cond = conditioner(target, _valid(), _mask())
    cond.context.sum().backward()
    assert target.grad is None


def test_detach_target_also_detaches_the_fixed_coordinates():
    """The carried fixed coordinates do not leak a gradient either."""
    conditioner = _conditioner(detach_target=True)
    target = TARGET_COORDS.clone().requires_grad_(True)
    cond = conditioner(target, _valid(), _mask())
    assert cond.fixed_xyz.grad_fn is None


def test_without_detach_target_the_gradient_reaches_the_target():
    """With ``detach_target=False`` the target is trainable through the context."""
    conditioner = _conditioner(detach_target=False)
    target = TARGET_COORDS.clone().requires_grad_(True)
    cond = conditioner(target, _valid(), _mask())
    cond.context.sum().backward()
    assert target.grad is not None
    assert torch.count_nonzero(target.grad).item() > 0


def test_detach_target_still_produces_a_usable_context():
    """Detaching does not change the values, only the graph."""
    detached = _conditioner(detach_target=True)
    attached = _conditioner(detach_target=False)
    ctx_detached = detached(TARGET_COORDS, _valid(), _mask()).context
    ctx_attached = attached(TARGET_COORDS, _valid(), _mask()).context
    assert torch.allclose(ctx_detached, ctx_attached, atol=1e-6)
