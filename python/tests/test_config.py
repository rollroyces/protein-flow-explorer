"""Tests for the pydantic configuration schemas in :mod:`pfe_design.config`.

The schemas are the first line of defence of the pipeline: an invalid shape
must be rejected before a single tensor is allocated.  These tests therefore
pin down the validation behaviour (ranges, cross-field rules) and the
immutability guarantees (frozen + ``extra="forbid"``).
"""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from pfe_design.config import (
    FlowIntegrationConfig,
    MetricConfig,
    ModelConfig,
    PackingBackend,
    PackingConfig,
    Solver,
    TargetConditioningConfig,
)

# --------------------------------------------------------------------- ModelConfig


def test_model_config_rejects_hidden_dim_not_divisible_by_n_heads():
    """A head count that does not divide the hidden dim is rejected."""
    with pytest.raises(ValidationError):
        ModelConfig(hidden_dim=128, n_heads=3)


def test_model_config_accepts_hidden_dim_divisible_by_n_heads():
    """A divisible configuration is accepted and keeps the requested values."""
    cfg = ModelConfig(hidden_dim=128, n_heads=8)
    assert cfg.hidden_dim == 128
    assert cfg.n_heads == 8


def test_model_config_head_dim_reports_per_head_width():
    """``head_dim`` is the hidden dim divided by the head count."""
    cfg = ModelConfig(hidden_dim=96, n_heads=8)
    assert cfg.head_dim == 12
    assert cfg.head_dim * cfg.n_heads == cfg.hidden_dim


def test_model_config_head_dim_scales_with_head_count():
    """Halving the head count doubles the width of a single head."""
    few_heads = ModelConfig(hidden_dim=64, n_heads=2)
    many_heads = ModelConfig(hidden_dim=64, n_heads=8)
    assert few_heads.head_dim == 32
    assert many_heads.head_dim == 8


def test_model_config_rejects_non_positive_hidden_dim():
    """The hidden dimension must be strictly positive."""
    with pytest.raises(ValidationError):
        ModelConfig(hidden_dim=0)


def test_model_config_rejects_dropout_of_one():
    """Dropout is bounded to ``[0, 1)`` - a rate of 1.0 is rejected."""
    with pytest.raises(ValidationError):
        ModelConfig(dropout=1.0)


def test_model_config_rejects_single_residue_type():
    """Index 0 is reserved for "unknown", so at least two types are required."""
    with pytest.raises(ValidationError):
        ModelConfig(n_residue_types=1)


def test_model_config_defaults_are_usable():
    """The defaults describe a small, valid network."""
    cfg = ModelConfig()
    assert cfg.hidden_dim % cfg.n_heads == 0
    assert cfg.head_dim == cfg.hidden_dim // cfg.n_heads
    assert cfg.n_layers >= 1


# ------------------------------------------------------------- FlowIntegrationConfig


def test_flow_integration_rejects_stochastic_with_zero_noise_scale():
    """A stochastic run without a noise scale is a configuration error."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(stochastic=True, noise_scale=0.0)


def test_flow_integration_accepts_stochastic_with_positive_noise_scale():
    """A stochastic run with a positive noise scale is accepted."""
    cfg = FlowIntegrationConfig(stochastic=True, noise_scale=0.5)
    assert cfg.stochastic is True
    assert cfg.noise_scale == pytest.approx(0.5)


def test_flow_integration_allows_zero_noise_when_deterministic():
    """``noise_scale`` is ignored by the deterministic solvers, so 0.0 is fine."""
    cfg = FlowIntegrationConfig(stochastic=False, noise_scale=0.0)
    assert cfg.stochastic is False


def test_flow_integration_rejects_negative_step_size():
    """The step size must be strictly positive."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(step_size=-0.01)


def test_flow_integration_rejects_zero_steps():
    """At least one integration step is required."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(n_steps=0)


def test_flow_integration_rejects_negative_seed():
    """A seed is a non-negative integer."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(seed=-1)


def test_flow_integration_rejects_unknown_dtype():
    """Only the two float dtypes the integrator implements are allowed."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(dtype="float16")


def test_flow_integration_rejects_non_positive_max_displacement():
    """The optional displacement cap must be strictly positive."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(max_displacement=0.0)


def test_flow_integration_parses_solver_from_string():
    """The solver is given as a string and resolved to the enum member."""
    cfg = FlowIntegrationConfig(solver="heun")
    assert cfg.solver is Solver.HEUN


# ------------------------------------------------------------------ other schemas


def test_metric_config_rejects_negative_clash_distance():
    """A clash threshold must be strictly positive."""
    with pytest.raises(ValidationError):
        MetricConfig(clash_distance=-1.0)


def test_metric_config_rejects_zero_clash_distance():
    """A zero clash threshold would flag every atom pair."""
    with pytest.raises(ValidationError):
        MetricConfig(clash_distance=0.0)


def test_metric_config_rejects_unknown_clash_mode():
    """The clash mode is a closed set of choices."""
    with pytest.raises(ValidationError):
        MetricConfig(clash_mode="whatever")


def test_target_conditioning_rejects_non_positive_pocket_cutoff():
    """The pocket cutoff must be strictly positive."""
    with pytest.raises(ValidationError):
        TargetConditioningConfig(pocket_cutoff=0.0)


def test_target_conditioning_rejects_zero_feature_dim():
    """The context width must be at least one."""
    with pytest.raises(ValidationError):
        TargetConditioningConfig(feature_dim=0)


def test_target_conditioning_rejects_zero_pocket_cap():
    """The pocket cap must allow at least one residue."""
    with pytest.raises(ValidationError):
        TargetConditioningConfig(max_pocket_residues=0)


def test_packing_config_rejects_zero_iterations():
    """The refinement hook needs at least one iteration."""
    with pytest.raises(ValidationError):
        PackingConfig(max_iterations=0)


def test_packing_config_defaults_to_mock_backend():
    """The packing hook has a backend that needs no external software."""
    assert PackingConfig().backend is PackingBackend.MOCK


# ------------------------------------------------------------- immutability rules


def test_schema_is_frozen_against_attribute_mutation():
    """Mutating a validated attribute is rejected."""
    cfg = ModelConfig()
    with pytest.raises(ValidationError):
        cfg.hidden_dim = 256


def test_schema_is_frozen_against_nested_mutation():
    """The frozen rule holds for every schema, not just ``ModelConfig``."""
    cfg = FlowIntegrationConfig()
    with pytest.raises(ValidationError):
        cfg.n_steps = 10


def test_schema_forbids_unknown_keyword():
    """A typo in a keyword is an error rather than a silently ignored setting."""
    with pytest.raises(ValidationError):
        ModelConfig(hidden_dim=128, n_heads=4, hidden_dimm=128)


def test_schema_forbids_unknown_keyword_on_integration_config():
    """The closed keyword set holds for the integration schema too."""
    with pytest.raises(ValidationError):
        FlowIntegrationConfig(n_steps=10, noize_scale=0.1)


def test_frozen_schema_is_hashable():
    """A frozen configuration can be attached to a provenance record."""
    cfg = ModelConfig(hidden_dim=64, n_heads=4)
    assert isinstance(hash(cfg), int)


def test_solver_enum_exposes_the_two_deterministic_solvers():
    """The integrator implements the Euler and the Heun solver."""
    assert Solver.EULER.value == "euler"
    assert Solver.HEUN.value == "heun"
    assert {s.value for s in Solver} == {"euler", "heun"}


def test_packing_backend_enum_exposes_known_backends():
    """The packing backends are a closed, explicit set."""
    assert PackingBackend.NONE.value == "none"
    assert PackingBackend.MOCK.value == "mock"
    assert {b.value for b in PackingBackend} == {
        "none",
        "mock",
        "rosetta",
        "openmm",
        "esmfold",
    }
