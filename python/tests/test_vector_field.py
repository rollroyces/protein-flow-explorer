"""Tests for :mod:`pfe_design.vector_field`.

The field is the network that predicts the flow-matching velocity of every
residue.  Two properties are load bearing for the rest of the pipeline:

* an untrained field predicts **exactly zero** velocity (the identity flow),
  which makes the masking guarantees of the integrator testable in isolation;
* invalid (unmodelled) residues are masked out of the attention keys, so they
  cannot influence the residues that *are* modelled.
"""

from __future__ import annotations

import pytest
import torch

from pfe_design.conditioning import PocketConditioner, build_complex_mask
from pfe_design.config import ModelConfig, TargetConditioningConfig
from pfe_design.vector_field import BackboneVelocityField, PairwiseBiasAttention, rbf_features

FEATURE_DIM = 16
HIDDEN_DIM = 32


def _configs(feature_dim: int = FEATURE_DIM):
    """A small model config and a matching conditioning config."""
    model_cfg = ModelConfig(hidden_dim=HIDDEN_DIM, n_heads=4, n_layers=2)
    cond_cfg = TargetConditioningConfig(feature_dim=feature_dim)
    return model_cfg, cond_cfg


def _field(feature_dim: int = FEATURE_DIM, seed: int = 0) -> BackboneVelocityField:
    """A freshly initialised, reproducible field."""
    torch.manual_seed(seed)
    model_cfg, cond_cfg = _configs(feature_dim)
    return BackboneVelocityField(model_cfg, cond_cfg)


def _live_field(feature_dim: int = FEATURE_DIM, seed: int = 0) -> BackboneVelocityField:
    """A field whose output layer is *not* zero.

    A fresh field is the identity flow (it predicts exactly zero), so any test
    that compares two velocities of a fresh field would compare two arrays of
    zeros and pass vacuously.  These tests need a field that actually reacts
    to its input.
    """
    field = _field(feature_dim, seed)
    torch.manual_seed(seed + 1000)
    torch.nn.init.normal_(field.to_velocity.weight, std=0.1)
    torch.nn.init.normal_(field.to_velocity.bias, std=0.1)
    return field


def _conditioner(feature_dim: int = FEATURE_DIM) -> PocketConditioner:
    torch.manual_seed(0)
    model_cfg, cond_cfg = _configs(feature_dim)
    return PocketConditioner(cond_cfg, model_cfg)


def _target(n_target: int = 4) -> torch.Tensor:
    """A target with one point per residue, spread along x."""
    coords = torch.arange(n_target * 3, dtype=torch.float32).reshape(1, n_target, 3)
    return coords * 2.0


def _cond(
    n_target: int = 4,
    n_design: int = 3,
    feature_dim: int = FEATURE_DIM,
    target: torch.Tensor | None = None,
    valid: torch.Tensor | None = None,
):
    """A conditioning output for a small target + design complex."""
    if target is None:
        target = _target(n_target)
    if valid is None:
        valid = torch.ones(1, n_target, dtype=torch.float32)
    mask = build_complex_mask(n_target, n_design, batch_size=1)
    return _conditioner(feature_dim)(target, valid, mask), target


def _x(n_target: int = 4, n_design: int = 3, seed: int = 0) -> torch.Tensor:
    torch.manual_seed(seed)
    return torch.randn(1, n_target + n_design, 3)


# ------------------------------------------------------------------- output shape


def test_field_output_has_shape_batch_residues_xyz():
    """The velocity has the same shape as the coordinates."""
    field = _field()
    cond, _ = _cond()
    v = field(_x(), None, cond)
    assert v.shape == (1, 4 + 3, 3)


def test_field_output_is_finite():
    """A forward pass produces finite numbers."""
    field = _live_field()
    cond, _ = _cond()
    v = field(_x(), None, cond)
    assert torch.count_nonzero(v).item() > 0
    assert bool(torch.isfinite(v).all())


def test_field_requires_a_conditioning_output():
    """The field is target-conditioned, so a missing context is an error."""
    field = _field()
    with pytest.raises(ValueError):
        field(_x(), None, None)


def test_field_rejects_coordinates_with_wrong_shape():
    """The coordinates must be ``(batch, n_residues, 3)``."""
    field = _field()
    cond, _ = _cond()
    with pytest.raises(ValueError):
        field(torch.randn(7, 3), None, cond)


def test_field_rejects_residue_type_out_of_range():
    """A residue-type index beyond the embedding table is rejected."""
    field = _field()
    cond, _ = _cond()
    n_res = 4 + 3
    residue_types = torch.zeros(1, n_res, dtype=torch.long)
    residue_types[0, 0] = field.cfg.n_residue_types
    with pytest.raises(ValueError):
        field(_x(), residue_types, cond)


def test_field_accepts_the_highest_valid_residue_type():
    """The last embedding slot is usable."""
    field = _field()
    cond, _ = _cond()
    n_res = 4 + 3
    residue_types = torch.zeros(1, n_res, dtype=torch.long)
    residue_types[0, 0] = field.cfg.n_residue_types - 1
    v = field(_x(), residue_types, cond)
    assert v.shape == (1, n_res, 3)


def test_field_rejects_residue_type_with_wrong_shape():
    """The residue-type tensor must have one entry per residue."""
    field = _field()
    cond, _ = _cond()
    with pytest.raises(ValueError):
        field(_x(), torch.zeros(1, 5, dtype=torch.long), cond)


def test_field_rejects_conditioning_with_wrong_context_width():
    """The context width must match the configured feature dimension."""
    model_cfg, cond_cfg = _configs(feature_dim=FEATURE_DIM)
    field = BackboneVelocityField(model_cfg, cond_cfg)
    cond, _ = _cond(feature_dim=FEATURE_DIM + 8)
    assert cond.context.shape[-1] != FEATURE_DIM
    with pytest.raises(ValueError):
        field(_x(), None, cond)


def test_field_rejects_conditioning_with_wrong_residue_count():
    """The context must describe exactly the residues that are being moved."""
    field = _field()
    cond, _ = _cond(n_target=4, n_design=3)
    with pytest.raises(ValueError):
        field(_x(n_target=6, n_design=3), None, cond)


# --------------------------------------------------------------- zero initialisation


def test_fresh_field_predicts_exactly_zero_velocity():
    """The zero-initialised output layer makes the untrained field the identity."""
    field = _field()
    cond, _ = _cond()
    v = field(_x(), None, cond)
    assert torch.count_nonzero(v).item() == 0


def test_fresh_field_predicts_zero_for_any_input():
    """Zero velocity does not depend on the input configuration."""
    field = _field()
    cond, _ = _cond()
    other = torch.randn(1, 7, 3) * 100.0
    assert torch.count_nonzero(field(other, None, cond)).item() == 0


def test_output_layer_weights_and_bias_are_zero():
    """The zero initialisation is visible in the parameters themselves."""
    field = _field()
    assert torch.count_nonzero(field.to_velocity.weight).item() == 0
    assert torch.count_nonzero(field.to_velocity.bias).item() == 0


# ---------------------------------------------------------------------- determinism


def test_field_is_deterministic_for_the_same_input():
    """Two calls with identical inputs give bit-identical velocities."""
    field = _live_field()
    cond, _ = _cond()
    x = _x()
    first = field(x, None, cond)
    second = field(x, None, cond)
    assert torch.count_nonzero(first).item() > 0
    assert torch.equal(first, second)


def test_field_is_deterministic_with_explicit_residue_types():
    """Determinism also holds when residue types are supplied."""
    field = _live_field()
    cond, _ = _cond()
    x = _x()
    residue_types = torch.arange(7).reshape(1, 7) % 20
    first = field(x, residue_types, cond)
    second = field(x, residue_types, cond)
    assert torch.equal(first, second)


def test_field_output_depends_on_the_input():
    """The field is not a constant function - it reacts to the coordinates."""
    field = _live_field()
    cond, _ = _cond()
    a = field(_x(seed=0), None, cond)
    b = field(_x(seed=1), None, cond)
    assert not torch.allclose(a, b, atol=1e-6)


# -------------------------------------------------------------- invalid residues


def test_far_invalid_residue_does_not_change_the_valid_velocities():
    """An unmodelled residue far away cannot influence the modelled ones."""
    field = _live_field()
    conditioner = _conditioner()

    n_target, n_design = 5, 3
    near = torch.tensor(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [2.0, 2.0, 0.0]]
    )
    far_a = torch.tensor([[500.0, 0.0, 0.0]])
    far_b = torch.tensor([[800.0, 0.0, 0.0]])
    target_a = torch.cat([near, far_a], dim=0).unsqueeze(0)
    target_b = torch.cat([near, far_b], dim=0).unsqueeze(0)
    # The fifth target residue is unmodelled.
    valid = torch.tensor([[1.0, 1.0, 1.0, 1.0, 0.0]])
    mask = build_complex_mask(n_target, n_design, batch_size=1)
    design = torch.zeros(1, n_design, 3)

    x_a = torch.cat([target_a, design], dim=1)
    x_b = torch.cat([target_b, design], dim=1)

    cond_a = conditioner(target_a, valid, mask)
    cond_b = conditioner(target_b, valid, mask)
    v_a = field(x_a, None, cond_a)
    v_b = field(x_b, None, cond_b)

    # The unmodelled residue really is masked out of the attention keys.
    assert cond_a.valid_mask[0, n_target - 1].item() == 0.0
    # The comparison is not vacuous: the field does produce a velocity.
    assert torch.count_nonzero(v_a[:, : n_target - 1]).item() > 0
    # ... and the modelled residues see exactly the same velocity.
    assert torch.equal(v_a[:, : n_target - 1], v_b[:, : n_target - 1])


def test_far_invalid_residue_does_not_change_the_designed_velocities():
    """The same holds for the designed residues."""
    field = _live_field()
    conditioner = _conditioner()

    n_target, n_design = 5, 3
    near = torch.tensor(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [2.0, 2.0, 0.0]]
    )
    target_a = torch.cat([near, torch.tensor([[500.0, 0.0, 0.0]])], dim=0).unsqueeze(0)
    target_b = torch.cat([near, torch.tensor([[800.0, 0.0, 0.0]])], dim=0).unsqueeze(0)
    valid = torch.tensor([[1.0, 1.0, 1.0, 1.0, 0.0]])
    mask = build_complex_mask(n_target, n_design, batch_size=1)
    design = torch.zeros(1, n_design, 3)

    cond_a = conditioner(target_a, valid, mask)
    cond_b = conditioner(target_b, valid, mask)
    v_a = field(torch.cat([target_a, design], dim=1), None, cond_a)
    v_b = field(torch.cat([target_b, design], dim=1), None, cond_b)
    assert torch.count_nonzero(v_a[:, n_target:]).item() > 0
    assert torch.equal(v_a[:, n_target:], v_b[:, n_target:])


def test_modelled_neighbours_do_change_the_velocity():
    """The control: a modelled neighbour *does* influence the velocity."""
    field = _live_field()
    conditioner = _conditioner()

    n_target, n_design = 4, 3
    target_a = torch.tensor(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [2.0, 2.0, 0.0]]
    ).unsqueeze(0)
    target_b = target_a + 6.0
    valid = torch.ones(1, n_target)
    mask = build_complex_mask(n_target, n_design, batch_size=1)
    design = torch.zeros(1, n_design, 3)

    v_a = field(
        torch.cat([target_a, design], dim=1),
        None,
        conditioner(target_a, valid, mask),
    )
    v_b = field(
        torch.cat([target_b, design], dim=1),
        None,
        conditioner(target_b, valid, mask),
    )
    assert torch.count_nonzero(v_a[:, :n_target]).item() > 0
    assert not torch.allclose(v_a[:, :n_target], v_b[:, :n_target], atol=1e-6)


# ------------------------------------------------------------------- rbf features


def test_rbf_features_have_one_column_per_bin():
    """The radial basis encoding is ``(..., n_bins)``."""
    distances = torch.rand(2, 5) * 20.0
    encoded = rbf_features(distances, n_bins=16, max_distance=20.0)
    assert encoded.shape == (2, 5, 16)


def test_rbf_features_are_maximal_at_zero_distance():
    """The first basis function peaks at a distance of zero."""
    distances = torch.zeros(1, 1)
    encoded = rbf_features(distances, n_bins=16, max_distance=20.0)
    assert encoded[0, 0, 0].item() == pytest.approx(1.0, abs=1e-6)


def test_rbf_features_vanish_beyond_the_maximum_distance():
    """Distances far beyond the cutoff contribute nothing."""
    distances = torch.full((1, 1), 500.0)
    encoded = rbf_features(distances, n_bins=16, max_distance=20.0)
    assert torch.count_nonzero(encoded).item() == 0


def test_the_attention_layer_works_without_a_pair_mask():
    """The pair mask is optional in the signature, so it has to be optional."""
    attention = PairwiseBiasAttention(hidden_dim=8, n_heads=2, n_rbf_bins=4)
    h = torch.randn(1, 3, 8)
    features = rbf_features(torch.rand(1, 3, 3), 4, 20.0)

    out = attention(h, features, None, None)

    assert out.shape == (1, 3, 8)


def test_a_float64_structure_goes_through_the_field():
    """A float64 structure must not fail inside the network.

    The conditioner was fixed to accept float64 coordinates; the network
    takes the same coordinates, and the pair features are built from them,
    so the fix has to reach all the way through.
    """
    torch.manual_seed(0)
    model_cfg, cond_cfg = _configs(FEATURE_DIM)
    field = BackboneVelocityField(model_cfg, cond_cfg)

    target = (_target() * 2.0).to(torch.float64)
    design = torch.randn(1, 3, 3, dtype=torch.float64)
    mask = build_complex_mask(4, 3, batch_size=1)
    valid = torch.ones(1, 4, dtype=torch.float64)
    cond = PocketConditioner(cond_cfg, model_cfg)(target, valid, mask)
    coords = torch.cat([target, design], dim=1)

    velocity = field(coords, None, cond)

    assert velocity.shape == (1, 7, 3)
    assert bool(torch.isfinite(velocity).all())


def test_the_pair_features_keep_the_dtype_that_was_asked_for():
    """A float64 distance tensor must not drag the basis to float64."""
    distances = torch.rand(2, 3, dtype=torch.float64) * 10.0

    features = rbf_features(distances, 4, 20.0, dtype=torch.float32)

    assert features.dtype == torch.float32
