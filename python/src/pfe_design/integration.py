"""Constrained ODE / SDE integration of the conditioned vector field.

The objective this module implements:

    Enforce zero velocity / zero update constraints on masked target residues
    during ODE/SDE integration trajectories.

Two independent guarantees are provided, because either one alone is
insufficient:

1. **Zero velocity.**  The predicted velocity is multiplied by the dynamic mask
   before it is used, so a fixed residue never receives a velocity, whatever
   the network predicts.
2. **Zero update.**  After every step the fixed residues are written back to
   their original coordinates, so even floating-point drift cannot move them.

Both guarantees hold for the deterministic (Euler / Heun) and the stochastic
(SDE) integrators.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import torch
from torch import Tensor

from .conditioning import ConditioningOutput, PocketConditioner, build_complex_mask
from .config import (
    FlowIntegrationConfig,
    ModelConfig,
    Solver,
    TargetConditioningConfig,
)
from .vector_field import BackboneVelocityField

__all__ = [
    "IntegrationResult",
    "VelocityModel",
    "apply_mask_constraints",
    "integrate",
    "sample_binder",
]

# A velocity model is any callable ``(x, residue_types, cond) -> velocity``.
VelocityModel = Callable[[Tensor, Tensor | None, ConditioningOutput], Tensor]

# How far the fixed residues of the starting state may differ from the fixed
# coordinates carried by the conditioning.  A millisecond of Angstrom is
# float32 rounding at protein coordinates; anything larger is a mismatch.
_FIXED_COORDINATE_TOLERANCE = 1e-3


@dataclass(frozen=True)
class IntegrationResult:
    """The outcome of an integration run.

    Attributes
    ----------
    coordinates:
        Final coordinates of the complex, ``(B, N, 3)``.
    velocities:
        Final velocities as applied, ``(B, N, 3)`` - zero on the fixed
        residues.
    n_steps:
        Number of integration steps that were taken.
    max_fixed_velocity:
        Largest absolute velocity *applied* to a fixed residue during the run.
        It is zero by construction, because the mask zeroes the velocity
        before it is used; it is reported so that a caller can see that the
        applied update really was zero.
    max_raw_fixed_velocity:
        Largest absolute velocity the model *wanted* on a fixed residue, before
        the mask was applied.  This is the interesting number: it says how hard
        the constraint is pulling, and it is what would have moved the target
        if the mask were not there.
    max_fixed_displacement:
        Largest displacement of a fixed residue from its starting position
        during the run.
    """

    coordinates: Tensor
    velocities: Tensor
    n_steps: int
    max_fixed_velocity: float
    max_raw_fixed_velocity: float
    max_fixed_displacement: float


def apply_mask_constraints(
    velocity: Tensor, cond: ConditioningOutput
) -> Tensor:
    """Zero the velocity of every fixed residue.

    Parameters
    ----------
    velocity:
        Predicted velocities, ``(B, N, 3)``.
    cond:
        Conditioning output carrying the fixed mask.

    Returns
    -------
    Tensor
        Velocities with the fixed residues set to exactly zero.
    """
    fixed = cond.fixed_mask.unsqueeze(-1).to(velocity.dtype)
    return velocity * (1.0 - fixed)


def _reimpose_constraints(x: Tensor, cond: ConditioningOutput) -> Tensor:
    """Write the fixed residues back onto their original coordinates."""
    fixed = cond.fixed_mask.unsqueeze(-1).to(torch.bool)
    return torch.where(fixed, cond.fixed_xyz.to(x.dtype), x)


def _evaluate(
    model: VelocityModel,
    x: Tensor,
    residue_types: Tensor | None,
    cond: ConditioningOutput,
) -> Tensor:
    """Evaluate the model.

    The raw velocity is returned unmasked; the caller applies the mask, so
    that the unmasked value can also be reported.
    """
    velocity = model(x, residue_types, cond)
    if velocity.shape != x.shape:
        raise ValueError(
            f"the velocity model returned {tuple(velocity.shape)}, expected {tuple(x.shape)}"
        )
    return velocity


def _validate_mask(cond: ConditioningOutput) -> None:
    """Reject a mask that is not a clean binary dynamic/fixed partition.

    A fractional value would be read as a *weight* by the velocity mask and
    as a *flag* by the coordinate write-back, and the two readings disagree:
    a residue with a mask of 0.5 would receive half a velocity and then be
    snapped back to a coordinate it never had.  That failure is silent, so it
    is rejected here instead.
    """
    dynamic = cond.dynamic_mask
    binary = (dynamic == 0) | (dynamic == 1)
    if not bool(binary.all()):
        raise ValueError(
            "the dynamic mask must be binary (1 = designed and free to move, "
            f"0 = fixed target); it contains values such as "
            f"{sorted(set(dynamic.flatten().tolist()))[:5]}"
        )
    if not torch.equal(cond.fixed_mask, 1.0 - dynamic):
        raise ValueError("the fixed mask must be the complement of the dynamic mask")


def _validate_start(x: Tensor, cond: ConditioningOutput) -> None:
    """Reject a starting state that does not match the conditioning.

    The integrator re-imposes the conditioning's fixed coordinates after every
    step, so a start whose target is somewhere else would have the target
    teleport on the first step.  That is not a subtle degradation, it is a
    discontinuity, and it is cheaper to refuse it.
    """
    if x.shape[1] != cond.n_residues:
        raise ValueError(
            f"x0 has {x.shape[1]} residues but the conditioning describes "
            f"{cond.n_residues}"
        )
    fixed = cond.fixed_mask.to(torch.bool)
    if not bool(fixed.any()):
        return
    difference = (x - cond.fixed_xyz.to(x.dtype)).abs() * fixed.unsqueeze(-1)
    worst = float(difference.max())
    if worst > _FIXED_COORDINATE_TOLERANCE:
        raise ValueError(
            f"the fixed residues of x0 differ from the conditioning's fixed "
            f"coordinates by up to {worst:.4f} Angstrom; the integrator would "
            f"move the target to the conditioning's coordinates, so start from "
            f"the same target"
        )


def integrate(
    x0: Tensor,
    model: VelocityModel,
    cond: ConditioningOutput,
    cfg: FlowIntegrationConfig | None = None,
    residue_types: Tensor | None = None,
) -> IntegrationResult:
    """Integrate the conditioned vector field under the mask constraints.

    Parameters
    ----------
    x0:
        Starting coordinates of the complex, ``(B, N, 3)``.
    model:
        Callable returning the velocity at ``(x, residue_types, cond)``.
    cond:
        Conditioning output; its fixed mask and fixed coordinates drive the
        constraints.
    cfg:
        Integration settings (steps, step size, solver, noise).
    residue_types:
        Optional residue-type indices, ``(B, N)``.

    Returns
    -------
    IntegrationResult
    """
    cfg = cfg or FlowIntegrationConfig()
    device = torch.device(cfg.device)
    dtype = getattr(torch, cfg.dtype)

    x = x0.to(device=device, dtype=dtype)
    if x.dim() != 3 or x.shape[-1] != 3:
        raise ValueError(f"x0 must have shape (batch, n_residues, 3); got {tuple(x0.shape)}")
    cond = cond.to(device)
    _validate_mask(cond)
    _validate_start(x, cond)

    generator = None
    if cfg.stochastic:
        generator = torch.Generator(device=device)
        if cfg.seed is not None:
            generator.manual_seed(int(cfg.seed))
        elif cfg.seed is None:
            generator.seed()

    dt = float(cfg.step_size)
    dynamic = cond.dynamic_mask.unsqueeze(-1).to(dtype)
    fixed = cond.fixed_mask.unsqueeze(-1).to(dtype)
    start = x.clone()

    max_fixed_velocity = 0.0
    max_raw_fixed_velocity = 0.0
    max_fixed_displacement = 0.0
    velocity = torch.zeros_like(x)

    for _ in range(cfg.n_steps):
        noise = None
        if cfg.stochastic:
            sigma = float(cfg.noise_scale)
            noise = (
                torch.randn(x.shape, device=device, dtype=dtype, generator=generator)
                * (dt ** 0.5)
                * sigma
                * dynamic
            )

        if cfg.solver is Solver.EULER:
            raw = _evaluate(model, x, residue_types, cond)
            velocity = apply_mask_constraints(raw, cond)
            x_next = x + dt * velocity
            if noise is not None:
                x_next = x_next + noise
            raw_velocities = (raw,)
        else:  # Heun
            raw1 = _evaluate(model, x, residue_types, cond)
            v1 = apply_mask_constraints(raw1, cond)
            # The predictor and the corrector share one Brownian increment, so
            # the diffusion term enters the update itself and not only the
            # point the second velocity is evaluated at.
            x_pred = x + dt * v1
            if noise is not None:
                x_pred = x_pred + noise
            if cfg.enforce_hard_constraints:
                x_pred = _reimpose_constraints(x_pred, cond)
            raw2 = _evaluate(model, x_pred, residue_types, cond)
            v2 = apply_mask_constraints(raw2, cond)
            velocity = 0.5 * (v1 + v2)
            x_next = x + dt * velocity
            if noise is not None:
                x_next = x_next + noise
            raw_velocities = (raw1, raw2)

        if cfg.max_displacement is not None:
            delta = x_next - x
            max_norm = torch.linalg.vector_norm(delta, dim=-1, keepdim=True)
            factor = (cfg.max_displacement / max_norm.clamp(min=1e-12)).clamp(max=1.0)
            x_next = x + delta * factor

        # Diagnostics: the applied velocity on a fixed residue is zero by
        # construction; the raw velocity is what the model asked for.
        fixed_velocity = velocity * fixed
        max_fixed_velocity = max(
            max_fixed_velocity, float(fixed_velocity.abs().max().item())
        )
        for raw in raw_velocities:
            raw_fixed = raw * fixed
            max_raw_fixed_velocity = max(
                max_raw_fixed_velocity, float(raw_fixed.abs().max().item())
            )

        if cfg.enforce_hard_constraints:
            x_next = _reimpose_constraints(x_next, cond)

        moved = (x_next - start) * fixed
        max_fixed_displacement = max(
            max_fixed_displacement, float(moved.abs().max().item())
        )
        x = x_next

    return IntegrationResult(
        coordinates=x,
        velocities=velocity,
        n_steps=cfg.n_steps,
        max_fixed_velocity=max_fixed_velocity,
        max_raw_fixed_velocity=max_raw_fixed_velocity,
        max_fixed_displacement=max_fixed_displacement,
    )


def sample_binder(
    target_xyz: Tensor,
    target_valid: Tensor,
    n_design: int,
    model: VelocityModel | None = None,
    cfg: FlowIntegrationConfig | None = None,
    cond_cfg: TargetConditioningConfig | None = None,
    model_cfg: ModelConfig | None = None,
    residue_types: Tensor | None = None,
    pocket_center: Tensor | None = None,
    init_xyz: Tensor | None = None,
) -> tuple[IntegrationResult, ConditioningOutput]:
    """Generate a binder against a target pocket.

    Parameters
    ----------
    target_xyz:
        Target backbone coordinates, ``(B, Nt, 3)``.
    target_valid:
        1 for a modelled target residue, ``(B, Nt)``.
    n_design:
        Number of residues to design.
    model:
        Velocity model; a default
        :class:`~pfe_design.vector_field.BackboneVelocityField` is created when
        omitted.
    cfg:
        Integration settings.
    cond_cfg, model_cfg:
        Settings for the conditioner and the network.
    residue_types:
        Residue types for the whole complex, ``(B, Nt + n_design)``.
    pocket_center:
        The site to aim at, ``(B, 3)`` or ``(3,)``.  When omitted the pocket
        is taken as the centroid of the target.
    init_xyz:
        Starting coordinates of the designed region, ``(B, n_design, 3)``.
        Random coordinates are used when omitted.

    Returns
    -------
    (IntegrationResult, ConditioningOutput)
    """
    cfg = cfg or FlowIntegrationConfig()
    cond_cfg = cond_cfg or TargetConditioningConfig()
    model_cfg = model_cfg or ModelConfig()

    device = torch.device(cfg.device)
    dtype = getattr(torch, cfg.dtype)
    target_xyz = target_xyz.to(device=device, dtype=dtype)
    target_valid = target_valid.to(device=device, dtype=dtype)
    batch, n_target, _ = target_xyz.shape

    if n_design <= 0:
        raise ValueError("n_design must be positive")

    mask = build_complex_mask(
        n_target, n_design, batch_size=batch, device=device, dtype=dtype
    )

    if init_xyz is None:
        generator = torch.Generator(device=device)
        if cfg.seed is not None:
            generator.manual_seed(int(cfg.seed))
        init_xyz = torch.randn(
            (batch, n_design, 3), device=device, dtype=dtype, generator=generator
        )
    init_xyz = init_xyz.to(device=device, dtype=dtype)

    if pocket_center is not None:
        pocket_center = pocket_center.to(device=device, dtype=dtype)
        if pocket_center.dim() == 1:
            pocket_center = pocket_center.unsqueeze(0).expand(batch, 3)
        if pocket_center.shape != (batch, 3):
            raise ValueError(
                f"pocket_center must have shape ({batch}, 3); got {tuple(pocket_center.shape)}"
            )

    x0 = torch.cat([target_xyz, init_xyz], dim=1)
    complex_xyz = x0

    conditioner = PocketConditioner(cond_cfg, model_cfg).to(device=device, dtype=dtype)
    cond = conditioner(
        target_xyz=target_xyz,
        target_valid=target_valid,
        dynamic_mask=mask,
        complex_xyz=complex_xyz,
        pocket_center=pocket_center,
    )

    if model is None:
        field = BackboneVelocityField(model_cfg, cond_cfg).to(device=device, dtype=dtype)
        model = field

    if residue_types is not None:
        residue_types = residue_types.to(device)

    result = integrate(x0, model, cond, cfg, residue_types=residue_types)
    return result, cond
