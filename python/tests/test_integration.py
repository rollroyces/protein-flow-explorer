"""Tests for :mod:`pfe_design.integration`.

This is the core guarantee of the module: during an ODE / SDE trajectory the
fixed target residues must receive **zero velocity** and **zero displacement**.
Both guarantees are checked independently:

* *zero velocity* - the predicted velocity is masked before it is used, so a
  fixed residue never moves even if the model asks it to;
* *zero update* - the fixed coordinates are re-imposed after every step, so
  even floating-point drift cannot move them.

The fixed residues are compared with ``torch.equal`` (bit-identical), not with
``allclose``: the constraint is exact, so the test should be exact too.
"""

from __future__ import annotations

import pytest
import torch

from pfe_design.conditioning import PocketConditioner, build_complex_mask
from pfe_design.config import (
    FlowIntegrationConfig,
    ModelConfig,
    Solver,
    TargetConditioningConfig,
)
from pfe_design.integration import apply_mask_constraints, integrate, sample_binder

N_TARGET = 4
N_DESIGN = 3
N_COMPLEX = N_TARGET + N_DESIGN

TARGET_COORDS = torch.tensor(
    [
        [0.0, 0.0, 0.0],
        [2.0, 0.0, 0.0],
        [4.0, 0.0, 0.0],
        [6.0, 0.0, 0.0],
    ],
    dtype=torch.float32,
).unsqueeze(0)

DESIGN_COORDS = torch.tensor(
    [
        [1.0, 1.0, 1.0],
        [3.0, 1.0, 1.0],
        [5.0, 1.0, 1.0],
    ],
    dtype=torch.float32,
).unsqueeze(0)


def _cond():
    """A conditioning output for a small target + design complex."""
    torch.manual_seed(0)
    cond_cfg = TargetConditioningConfig(feature_dim=16)
    model_cfg = ModelConfig(hidden_dim=32, n_heads=4, n_layers=2)
    conditioner = PocketConditioner(cond_cfg, model_cfg)
    valid = torch.ones(1, N_TARGET, dtype=torch.float32)
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=1)
    cond = conditioner(TARGET_COORDS, valid, mask)
    return cond


def _x0() -> torch.Tensor:
    """The starting coordinates of the complex."""
    return torch.cat([TARGET_COORDS, DESIGN_COORDS], dim=1)


def _constant_model(speed: float = 0.5):
    """A velocity model that always asks for the same constant velocity."""

    def model(x, residue_types, cond):
        return torch.full_like(x, speed)

    return model


def _axis_model(speed: float = 10.0):
    """A velocity model that pushes everything along +x."""

    def model(x, residue_types, cond):
        v = torch.zeros_like(x)
        v[..., 0] = speed
        return v

    return model


def _huge_on_fixed_model(huge: float = 1.0e4, modest: float = 0.25):
    """A model that tries to blow the fixed residues away."""

    def model(x, residue_types, cond):
        fixed = cond.fixed_mask.unsqueeze(-1).to(x.dtype)
        return modest * torch.ones_like(x) + huge * fixed

    return model


def _cfg(**kwargs) -> FlowIntegrationConfig:
    """A short deterministic integration configuration."""
    base = dict(n_steps=10, step_size=0.1, solver=Solver.EULER)
    base.update(kwargs)
    return FlowIntegrationConfig(**base)


def _stochastic_cfg(seed: int = 1234, **kwargs) -> FlowIntegrationConfig:
    base = dict(
        n_steps=10,
        step_size=0.1,
        solver=Solver.EULER,
        stochastic=True,
        noise_scale=0.5,
        seed=seed,
    )
    base.update(kwargs)
    return FlowIntegrationConfig(**base)


# ------------------------------------------------------- fixed residues stay fixed


def test_euler_keeps_fixed_residues_bit_identical():
    """The Euler solver leaves the target coordinates bit-identical."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(), cond, _cfg(solver=Solver.EULER))
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])


def test_heun_keeps_fixed_residues_bit_identical():
    """The Heun solver leaves the target coordinates bit-identical."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(), cond, _cfg(solver=Solver.HEUN))
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])


def test_stochastic_sde_keeps_fixed_residues_bit_identical():
    """The stochastic solver leaves the target coordinates bit-identical."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(), cond, _stochastic_cfg())
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])


def test_stochastic_heun_keeps_fixed_residues_bit_identical():
    """The noisy Heun solver also leaves the target coordinates bit-identical."""
    cond = _cond()
    x0 = _x0()
    result = integrate(
        x0, _constant_model(), cond, _stochastic_cfg(solver=Solver.HEUN)
    )
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])


def test_fixed_residues_stay_fixed_without_the_write_back():
    """The zero-velocity guarantee alone is enough to hold the target in place."""
    cond = _cond()
    x0 = _x0()
    result = integrate(
        x0,
        _constant_model(),
        cond,
        _cfg(enforce_hard_constraints=False),
    )
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])


# ----------------------------------------------------------- reported diagnostics


def test_euler_reports_zero_max_fixed_velocity():
    """The reported maximum velocity on the fixed residues is exactly zero."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg(solver=Solver.EULER))
    assert result.max_fixed_velocity == 0.0


def test_heun_reports_zero_max_fixed_velocity():
    """The reported maximum velocity on the fixed residues is exactly zero."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg(solver=Solver.HEUN))
    assert result.max_fixed_velocity == 0.0


def test_stochastic_reports_zero_max_fixed_velocity():
    """The reported maximum velocity on the fixed residues is exactly zero."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _stochastic_cfg())
    assert result.max_fixed_velocity == 0.0


def test_euler_reports_zero_max_fixed_displacement():
    """The reported maximum displacement of the fixed residues is exactly zero."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg(solver=Solver.EULER))
    assert result.max_fixed_displacement == 0.0


def test_heun_reports_zero_max_fixed_displacement():
    """The reported maximum displacement of the fixed residues is exactly zero."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg(solver=Solver.HEUN))
    assert result.max_fixed_displacement == 0.0


def test_stochastic_reports_zero_max_fixed_displacement():
    """The reported maximum displacement of the fixed residues is exactly zero."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _stochastic_cfg())
    assert result.max_fixed_displacement == 0.0


# ----------------------------------------------------- the designed region moves


def test_euler_moves_the_designed_residues():
    """The designed residues move - the constraint is not trivially satisfied."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(0.5), cond, _cfg(solver=Solver.EULER))
    moved = result.coordinates[:, N_TARGET:] - x0[:, N_TARGET:]
    assert moved.norm().item() > 0.1


def test_heun_moves_the_designed_residues():
    """The Heun solver also moves the designed residues."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(0.5), cond, _cfg(solver=Solver.HEUN))
    moved = result.coordinates[:, N_TARGET:] - x0[:, N_TARGET:]
    assert moved.norm().item() > 0.1


def test_stochastic_moves_the_designed_residues():
    """The stochastic solver also moves the designed residues."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(0.5), cond, _stochastic_cfg())
    moved = result.coordinates[:, N_TARGET:] - x0[:, N_TARGET:]
    assert moved.norm().item() > 0.1


def test_euler_displacement_matches_the_expected_constant_velocity():
    """A constant velocity integrates to ``n_steps * step_size * speed``."""
    cond = _cond()
    x0 = _x0()
    cfg = _cfg(n_steps=10, step_size=0.1, solver=Solver.EULER)
    result = integrate(x0, _constant_model(0.5), cond, cfg)
    expected = 10 * 0.1 * 0.5
    assert result.coordinates[0, N_TARGET, 0].item() == pytest.approx(
        x0[0, N_TARGET, 0].item() + expected, abs=1e-4
    )


def test_a_zero_velocity_field_leaves_the_complex_untouched():
    """The untrained field is the identity flow - nothing moves at all."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _constant_model(0.0), cond, _cfg())
    assert torch.equal(result.coordinates, x0)


# ------------------------------------------------------- apply_mask_constraints


def test_apply_mask_constraints_zeroes_the_fixed_velocity():
    """The velocity of a fixed residue is exactly zero."""
    cond = _cond()
    velocity = torch.randn(1, N_COMPLEX, 3)
    constrained = apply_mask_constraints(velocity, cond)
    assert torch.count_nonzero(constrained[:, :N_TARGET]).item() == 0


def test_apply_mask_constraints_leaves_the_dynamic_velocity_untouched():
    """The velocity of a designed residue is passed through unchanged."""
    cond = _cond()
    velocity = torch.randn(1, N_COMPLEX, 3)
    constrained = apply_mask_constraints(velocity, cond)
    assert torch.equal(constrained[:, N_TARGET:], velocity[:, N_TARGET:])


def test_apply_mask_constraints_does_not_modify_its_input():
    """The masking is not done in place on the caller's tensor."""
    cond = _cond()
    velocity = torch.randn(1, N_COMPLEX, 3)
    original = velocity.clone()
    apply_mask_constraints(velocity, cond)
    assert torch.equal(velocity, original)


# ------------------------------------------------------------ hostile velocity models


def test_huge_velocity_on_fixed_residues_is_fully_constrained():
    """Even a model that asks for an enormous velocity cannot move the target."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _huge_on_fixed_model(), cond, _cfg(solver=Solver.EULER))
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])
    assert result.max_fixed_velocity == 0.0
    assert result.max_fixed_displacement == 0.0


def test_huge_velocity_on_fixed_residues_is_constrained_by_heun():
    """The same holds for the Heun solver."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _huge_on_fixed_model(), cond, _cfg(solver=Solver.HEUN))
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])
    assert result.max_fixed_velocity == 0.0


def test_huge_velocity_on_fixed_residues_still_moves_the_design():
    """The designed residues are still free to move under the same model."""
    cond = _cond()
    x0 = _x0()
    result = integrate(x0, _huge_on_fixed_model(), cond, _cfg())
    moved = result.coordinates[:, N_TARGET:] - x0[:, N_TARGET:]
    assert moved.norm().item() > 0.1


# ------------------------------------------------------------- max_displacement cap


def test_max_displacement_caps_the_per_step_movement():
    """No residue moves further than the configured cap in one step."""
    cond = _cond()
    x0 = _x0()
    cfg = _cfg(n_steps=4, step_size=0.1, max_displacement=0.05)
    result = integrate(x0, _axis_model(10.0), cond, cfg)
    # Without the cap the single step would have been 1.0 Angstrom.
    total = (result.coordinates[:, N_TARGET:] - x0[:, N_TARGET:]).norm(dim=-1)
    assert torch.allclose(total, torch.full_like(total, 0.2), atol=1e-4)


def test_max_displacement_is_never_exceeded():
    """The per-step displacement never exceeds the cap."""
    cond = _cond()
    x0 = _x0()
    cap = 0.05
    cfg = _cfg(n_steps=4, step_size=0.1, max_displacement=cap)
    result = integrate(x0, _axis_model(10.0), cond, cfg)
    step = (result.coordinates[0, N_TARGET, 0] - x0[0, N_TARGET, 0]).item() / 4
    assert step == pytest.approx(cap, abs=1e-5)


def test_max_displacement_does_not_disturb_the_fixed_residues():
    """The cap is applied to the designed residues only."""
    cond = _cond()
    x0 = _x0()
    cfg = _cfg(n_steps=4, step_size=0.1, max_displacement=0.05)
    result = integrate(x0, _axis_model(10.0), cond, cfg)
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])


def test_max_displacement_leaves_small_moves_untouched():
    """A movement below the cap is not scaled down."""
    cond = _cond()
    x0 = _x0()
    cfg = _cfg(n_steps=1, step_size=0.1, max_displacement=5.0)
    result = integrate(x0, _constant_model(0.5), cond, cfg)
    expected = x0[0, N_TARGET, 0].item() + 0.05
    assert result.coordinates[0, N_TARGET, 0].item() == pytest.approx(
        expected, abs=1e-5
    )


# ------------------------------------------------------------- stochastic runs


def test_stochastic_run_is_reproducible_with_the_same_seed():
    """The same seed gives the same trajectory."""
    cond = _cond()
    x0 = _x0()
    first = integrate(x0, _constant_model(0.5), cond, _stochastic_cfg(seed=1234))
    second = integrate(x0, _constant_model(0.5), cond, _stochastic_cfg(seed=1234))
    assert torch.equal(first.coordinates, second.coordinates)


def test_stochastic_run_differs_with_a_different_seed():
    """A different seed explores a different noise realisation."""
    cond = _cond()
    x0 = _x0()
    first = integrate(x0, _constant_model(0.5), cond, _stochastic_cfg(seed=1234))
    second = integrate(x0, _constant_model(0.5), cond, _stochastic_cfg(seed=4321))
    assert not torch.allclose(first.coordinates, second.coordinates, atol=1e-6)


def test_stochastic_run_differs_from_the_deterministic_run():
    """The diffusion term actually contributes to the trajectory."""
    cond = _cond()
    x0 = _x0()
    deterministic = integrate(x0, _constant_model(0.5), cond, _cfg())
    stochastic = integrate(x0, _constant_model(0.5), cond, _stochastic_cfg())
    assert not torch.allclose(
        deterministic.coordinates[:, N_TARGET:],
        stochastic.coordinates[:, N_TARGET:],
        atol=1e-6,
    )


# ------------------------------------------------------------------- result record


def test_integration_result_reports_the_number_of_steps():
    """The result records how many steps were taken."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg(n_steps=7))
    assert result.n_steps == 7


def test_integration_result_reports_the_final_velocities():
    """The result carries the final velocity field."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg())
    assert result.velocities.shape == (1, N_COMPLEX, 3)


def test_integration_result_coordinates_keep_the_complex_shape():
    """The final coordinates have the shape of the complex."""
    cond = _cond()
    result = integrate(_x0(), _constant_model(), cond, _cfg())
    assert result.coordinates.shape == (1, N_COMPLEX, 3)


def test_integration_rejects_coordinates_with_wrong_shape():
    """The starting coordinates must be ``(batch, n_residues, 3)``."""
    cond = _cond()
    with pytest.raises(ValueError):
        integrate(torch.zeros(7, 3), _constant_model(), cond, _cfg())


def test_integration_rejects_a_model_with_the_wrong_output_shape():
    """A model that returns the wrong shape is a programming error."""

    def bad_model(x, residue_types, cond):
        return torch.zeros_like(x)[..., :2]

    cond = _cond()
    with pytest.raises(ValueError):
        integrate(_x0(), bad_model, cond, _cfg())


def test_integration_keeps_every_row_of_a_batch_fixed():
    """The guarantee holds for every complex in the batch."""
    torch.manual_seed(0)
    cond_cfg = TargetConditioningConfig(feature_dim=16)
    model_cfg = ModelConfig(hidden_dim=32, n_heads=4, n_layers=2)
    conditioner = PocketConditioner(cond_cfg, model_cfg)
    target = torch.cat([TARGET_COORDS, TARGET_COORDS + 10.0], dim=0)
    valid = torch.ones(2, N_TARGET)
    mask = build_complex_mask(N_TARGET, N_DESIGN, batch_size=2)
    cond = conditioner(target, valid, mask)
    design = torch.cat([DESIGN_COORDS, DESIGN_COORDS + 10.0], dim=0)
    x0 = torch.cat([target, design], dim=1)

    result = integrate(x0, _constant_model(), cond, _cfg())
    assert result.coordinates.shape == (2, N_COMPLEX, 3)
    assert torch.equal(result.coordinates[:, :N_TARGET], x0[:, :N_TARGET])
    assert result.max_fixed_velocity == 0.0
    assert result.max_fixed_displacement == 0.0


# ------------------------------------------------------------------ sample_binder


def _sample_cfg(**kwargs) -> FlowIntegrationConfig:
    base = dict(n_steps=10, step_size=0.1, seed=7)
    base.update(kwargs)
    return FlowIntegrationConfig(**base)


def test_sample_binder_produces_a_complex_of_the_right_length():
    """The sampled complex is the target followed by the designed region."""
    result, cond = sample_binder(
        target_xyz=TARGET_COORDS,
        target_valid=torch.ones(1, N_TARGET),
        n_design=N_DESIGN,
        model=_constant_model(0.5),
        cfg=_sample_cfg(),
        init_xyz=DESIGN_COORDS,
    )
    assert result.coordinates.shape == (1, N_TARGET + N_DESIGN, 3)
    assert cond.n_target == N_TARGET
    assert cond.n_design == N_DESIGN


def test_sample_binder_keeps_the_target_fixed():
    """The target residues come out bit-identical."""
    result, _ = sample_binder(
        target_xyz=TARGET_COORDS,
        target_valid=torch.ones(1, N_TARGET),
        n_design=N_DESIGN,
        model=_constant_model(0.5),
        cfg=_sample_cfg(),
        init_xyz=DESIGN_COORDS,
    )
    assert torch.equal(result.coordinates[:, :N_TARGET], TARGET_COORDS[:, :N_TARGET])


def test_sample_binder_moves_the_designed_region():
    """The designed residues are actually moved by the sampler."""
    result, _ = sample_binder(
        target_xyz=TARGET_COORDS,
        target_valid=torch.ones(1, N_TARGET),
        n_design=N_DESIGN,
        model=_constant_model(0.5),
        cfg=_sample_cfg(),
        init_xyz=DESIGN_COORDS,
    )
    moved = result.coordinates[:, N_TARGET:] - DESIGN_COORDS
    assert moved.norm().item() > 0.1


def test_sample_binder_reports_zero_constraint_violations():
    """The sampler reports the same guarantees as a plain integration run."""
    result, _ = sample_binder(
        target_xyz=TARGET_COORDS,
        target_valid=torch.ones(1, N_TARGET),
        n_design=N_DESIGN,
        model=_huge_on_fixed_model(),
        cfg=_sample_cfg(),
        init_xyz=DESIGN_COORDS,
    )
    assert result.max_fixed_velocity == 0.0
    assert result.max_fixed_displacement == 0.0
    assert torch.equal(result.coordinates[:, :N_TARGET], TARGET_COORDS[:, :N_TARGET])


def test_sample_binder_is_reproducible_with_a_seed():
    """A seeded run gives the same complex twice."""
    kwargs = dict(
        target_xyz=TARGET_COORDS,
        target_valid=torch.ones(1, N_TARGET),
        n_design=N_DESIGN,
        model=_constant_model(0.5),
        cfg=_sample_cfg(seed=7),
    )
    first, _ = sample_binder(**kwargs)
    second, _ = sample_binder(**kwargs)
    assert torch.equal(first.coordinates, second.coordinates)


def test_sample_binder_rejects_a_non_positive_design_length():
    """At least one residue has to be designed."""
    with pytest.raises(ValueError):
        sample_binder(
            target_xyz=TARGET_COORDS,
            target_valid=torch.ones(1, N_TARGET),
            n_design=0,
            model=_constant_model(),
            cfg=_sample_cfg(),
        )


def test_sample_binder_uses_a_zero_velocity_field_by_default():
    """Without a model the sampler uses the untrained, zero-velocity field."""
    result, _ = sample_binder(
        target_xyz=TARGET_COORDS,
        target_valid=torch.ones(1, N_TARGET),
        n_design=N_DESIGN,
        cfg=_sample_cfg(),
        init_xyz=DESIGN_COORDS,
    )
    assert torch.equal(result.coordinates[:, N_TARGET:], DESIGN_COORDS)
