"""Regression tests for the defects found in review.

Every test in this file corresponds to a specific defect that was found after
the first version of the package was working and passing its own tests. They
are collected together so that the reason each one exists is written down
next to it.

The defects, in the order they were found:

1. A fractional mask value was read as a *weight* by the velocity mask and as
   a *flag* by the coordinate write-back. A residue with a mask of 0.5 was
   therefore written back to the origin - silently.
2. The stochastic Heun solver added its noise only to the predictor, so the
   diffusion term never reached the update and the "stochastic" solver was
   bit-identical to the deterministic one.
3. ``max_fixed_velocity`` was computed from the already-masked velocity, so it
   was zero by construction and could never report a violation.
4. The pocket fallback was a batch-global test, so one member of a batch could
   change the pocket of another.
5. The buried surface area summed the complex over *all* atoms, including
   atoms outside the two partners, which could make the burial negative.
6. A batch of structures with per-structure radii broadcast the radii against
   the wrong axis and counted some pairs twice.
7. An unmodelled residue next door changed the velocities of the modelled
   residues, because the pair features were built from the unmasked
   coordinates.
8. The integrator accepted a starting state whose residue count did not match
   the conditioning, and accepted a target that was somewhere else entirely.
9. The surface-area calculation only worked for systems of thirty-two atoms
   or fewer: the self-occlusion mask was built on the chunk axis instead of
   the atom axis, so the moment a system needed more than one chunk the
   shapes stopped matching and it raised. Every earlier test used a system
   small enough to fit in one chunk.
"""

from __future__ import annotations

import pytest
import torch

from pfe_design import (
    FlowIntegrationConfig,
    MetricConfig,
    ModelConfig,
    Solver,
    TargetConditioningConfig,
    build_complex_mask,
    buried_surface_area,
    clash_penalty,
    integrate,
)
from pfe_design.conditioning import PocketConditioner
from pfe_design.metrics import sasa
from pfe_design.vector_field import BackboneVelocityField


def _target(n_target: int = 4) -> torch.Tensor:
    return torch.arange(n_target * 3, dtype=torch.float32).reshape(1, n_target, 3) * 2.0


def _conditioning(n_target: int = 4, n_design: int = 3):
    target = _target(n_target)
    valid = torch.ones(1, n_target)
    mask = build_complex_mask(n_target, n_design, batch_size=1)
    conditioner = PocketConditioner(TargetConditioningConfig(feature_dim=16))
    return conditioner(target, valid, mask), target, mask


# 1 - a fractional mask value must be refused, not silently misread


def test_a_fractional_mask_is_rejected():
    """0.5 would be a weight to one constraint and a flag to the other."""
    conditioner = PocketConditioner()
    target = _target()
    mask = build_complex_mask(4, 3, batch_size=1)
    mask[0, 4] = 0.5

    with pytest.raises(ValueError, match="binary"):
        conditioner(target, torch.ones(1, 4), mask)


def test_the_integrator_rejects_a_fractional_mask():
    """A mask that did not survive the conditioner is refused there too."""
    cond, target, mask = _conditioning()
    broken = type(cond)(
        context=cond.context,
        pocket_center=cond.pocket_center,
        pocket_radius=cond.pocket_radius,
        dynamic_mask=cond.dynamic_mask.clone(),
        fixed_mask=cond.fixed_mask.clone(),
        valid_mask=cond.valid_mask,
        fixed_xyz=cond.fixed_xyz,
        pocket_mask=cond.pocket_mask,
    )
    broken.dynamic_mask[0, 4] = 0.5
    broken.fixed_mask[0, 4] = 0.5

    x0 = torch.cat([target, torch.randn(1, 3, 3)], dim=1)
    with pytest.raises(ValueError, match="binary"):
        integrate(x0, lambda x, rt, c: torch.zeros_like(x), broken, FlowIntegrationConfig(n_steps=2))


# 2 - the stochastic Heun solver must actually diffuse


def test_stochastic_heun_differs_from_deterministic_heun():
    """The diffusion term has to reach the update, not only the predictor."""
    cond, target, _ = _conditioning()
    x0 = torch.cat([target, torch.randn(1, 3, 3)], dim=1)
    def model(x, rt, c):
        return x * 0.1 + 1.0

    deterministic = integrate(
        x0, model, cond, FlowIntegrationConfig(n_steps=5, step_size=0.1, solver=Solver.HEUN)
    )
    stochastic = integrate(
        x0,
        model,
        cond,
        FlowIntegrationConfig(
            n_steps=5, step_size=0.1, solver=Solver.HEUN, stochastic=True, noise_scale=0.5, seed=7
        ),
    )

    assert not torch.allclose(deterministic.coordinates, stochastic.coordinates)


def test_stochastic_heun_depends_on_the_seed():
    """Two seeds have to give two trajectories."""
    cond, target, _ = _conditioning()
    x0 = torch.cat([target, torch.randn(1, 3, 3)], dim=1)
    def model(x, rt, c):
        return x * 0.1 + 1.0

    first = integrate(
        x0,
        model,
        cond,
        FlowIntegrationConfig(
            n_steps=5, step_size=0.1, solver=Solver.HEUN, stochastic=True, noise_scale=0.5, seed=7
        ),
    )
    second = integrate(
        x0,
        model,
        cond,
        FlowIntegrationConfig(
            n_steps=5, step_size=0.1, solver=Solver.HEUN, stochastic=True, noise_scale=0.5, seed=99
        ),
    )

    assert not torch.allclose(first.coordinates, second.coordinates)


def test_the_noise_never_moves_a_fixed_residue():
    """The diffusion term is masked like the velocity is."""
    cond, target, _ = _conditioning()
    x0 = torch.cat([target, torch.randn(1, 3, 3)], dim=1)
    def model(x, rt, c):
        return torch.ones_like(x)

    result = integrate(
        x0,
        model,
        cond,
        FlowIntegrationConfig(
            n_steps=20, step_size=0.05, stochastic=True, noise_scale=1.0, seed=3
        ),
    )

    assert torch.equal(result.coordinates[:, :4], target)


# 3 - the diagnostic has to be able to report something


def test_the_raw_velocity_diagnostic_reports_what_the_model_wanted():
    """A hostile model has to show up in the diagnostic."""
    cond, target, _ = _conditioning()
    x0 = torch.cat([target, torch.randn(1, 3, 3)], dim=1)

    def hostile(x, rt, c):
        out = torch.ones_like(x)
        out[:, :4] = 1e4
        return out

    result = integrate(x0, hostile, cond, FlowIntegrationConfig(n_steps=3, step_size=0.01))

    assert result.max_raw_fixed_velocity > 1e3
    assert result.max_fixed_velocity == 0.0
    assert torch.equal(result.coordinates[:, :4], target)


# 4 - a batch member must not change another member's pocket


def test_a_batch_members_pocket_is_independent_of_its_neighbours():
    """The fallback for an empty pocket is applied per structure."""
    near = torch.tensor([[[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [4.0, 0.0, 0.0], [6.0, 0.0, 0.0]]])
    far = near + 100.0
    design_far = far + 0.5
    design_near = near + 0.5

    batched_target = torch.cat([near, near], dim=0)
    batched_complex = torch.cat(
        [batched_target, torch.cat([design_far, design_near], dim=0)[:, :3]], dim=1
    )
    mask = build_complex_mask(4, 3, batch_size=2)
    conditioner = PocketConditioner()

    batched = conditioner(
        batched_target, torch.ones(2, 4), mask, complex_xyz=batched_complex
    )
    solo = conditioner(
        near, torch.ones(1, 4), mask[:1], complex_xyz=batched_complex[:1]
    )

    assert torch.allclose(batched.pocket_center[0], solo.pocket_center[0])
    assert int(batched.pocket_mask[0].sum()) > 0


# 5 - the buried surface area has to stay physical


def test_buried_surface_area_ignores_atoms_outside_the_two_partners():
    """A third group must not make the burial negative."""
    generator = torch.Generator().manual_seed(11)
    coords = torch.randn((30, 3), generator=generator, dtype=torch.float64) * 3.0
    radii = torch.full((30,), 1.7, dtype=torch.float64)

    result = buried_surface_area(
        coords, radii, torch.arange(10), torch.arange(10, 20), MetricConfig()
    )

    assert result.bsa > 0.0
    assert result.bsa == pytest.approx(
        result.sasa_free_a + result.sasa_free_b - result.sasa_complex, rel=1e-9
    )


def test_buried_surface_area_is_unchanged_by_an_unrelated_group():
    """The same two partners give the same answer whatever else is present."""
    generator = torch.Generator().manual_seed(12)
    partners = torch.randn((20, 3), generator=generator, dtype=torch.float64) * 3.0
    extra = torch.randn((10, 3), generator=generator, dtype=torch.float64) * 3.0
    radii = torch.full((30,), 1.7, dtype=torch.float64)

    alone = buried_surface_area(
        partners, radii[:20], torch.arange(10), torch.arange(10, 20), MetricConfig()
    )
    with_extra = buried_surface_area(
        torch.cat([partners, extra]),
        radii,
        torch.arange(10),
        torch.arange(10, 20),
        MetricConfig(),
    )

    assert with_extra.bsa == pytest.approx(alone.bsa, rel=1e-9)


# 6 - per-structure radii in a batch


def test_per_structure_radii_are_not_broadcast_against_the_wrong_axis():
    """Two structures with their own radii get their own penalty."""
    coords = torch.stack(
        [torch.tensor([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]], dtype=torch.float64)] * 2
    )
    radii = torch.full((2, 2), 1.7, dtype=torch.float64)

    values = clash_penalty(
        coords, MetricConfig(clash_distance=2.0, clash_mode="vdw"), radii=radii
    )

    assert values.shape == (2,)
    assert float(values[0]) == pytest.approx(4.0, rel=1e-9)
    assert float(values[1]) == pytest.approx(4.0, rel=1e-9)


# 7 - an unmodelled residue must not reach the modelled ones


def test_an_unmodelled_residue_cannot_change_the_modelled_velocities():
    """A residue that is not in the structure has no coordinates to speak of."""
    torch.manual_seed(0)
    model_cfg = ModelConfig(hidden_dim=32, n_layers=2, n_heads=4)
    field = BackboneVelocityField(model_cfg, TargetConditioningConfig(feature_dim=16))

    target = torch.randn(1, 5, 3) * 5.0
    design = torch.randn(1, 3, 3)
    mask = build_complex_mask(5, 3, batch_size=1)
    valid = torch.tensor([[1.0, 1.0, 1.0, 1.0, 0.0]])  # the last target residue is not modelled
    cond = PocketConditioner(TargetConditioningConfig(feature_dim=16))(target, valid, mask)

    near = torch.cat([target, design], dim=1)
    far = near.clone()
    far[0, 4] = torch.tensor([100.0, 100.0, 100.0])

    velocity_near = field(near, None, cond)
    velocity_far = field(far, None, cond)

    modelled = [0, 1, 2, 3, 5, 6, 7]
    difference = (velocity_near - velocity_far).abs()[:, modelled].max()
    assert float(difference.detach()) < 1e-6


# 8 - the integrator refuses an inconsistent start


def test_the_integrator_refuses_a_start_of_the_wrong_length():
    """A conditioning for seven residues does not describe a one-residue state."""
    cond, _target, _ = _conditioning()

    with pytest.raises(ValueError, match="residues"):
        integrate(
            torch.randn(1, 1, 3),
            lambda x, rt, c: torch.zeros_like(x),
            cond,
            FlowIntegrationConfig(n_steps=2),
        )


def test_the_integrator_refuses_a_target_that_is_somewhere_else():
    """The fixed residues of the start must be the ones the conditioning names."""
    cond, target, _ = _conditioning()
    x0 = torch.cat([target + 50.0, torch.randn(1, 3, 3)], dim=1)

    with pytest.raises(ValueError, match="differ"):
        integrate(
            x0,
            lambda x, rt, c: torch.zeros_like(x),
            cond,
            FlowIntegrationConfig(n_steps=2),
        )


# 9 - the surface-area calculation must not depend on the chunk size


def test_sasa_is_the_same_however_the_atoms_are_chunked():
    """A system larger than one chunk must give the same answer.

    The self-occlusion mask used to be built on the chunk axis rather than
    the atom axis, so the calculation only worked while the whole system
    fitted in a single chunk - which every small test did.
    """
    generator = torch.Generator().manual_seed(5)
    coords = torch.randn((80, 3), generator=generator, dtype=torch.float64) * 4.0
    radii = torch.full((80,), 1.70, dtype=torch.float64)

    one_chunk = sasa(coords, radii, 1.4, 60, chunk=256)
    many_chunks = sasa(coords, radii, 1.4, 60, chunk=8)

    assert one_chunk.shape == (80,)
    assert torch.allclose(one_chunk, many_chunks, atol=1e-9)


def test_sasa_handles_a_system_larger_than_one_chunk():
    """A hundred atoms is a small protein, and it has to work."""
    generator = torch.Generator().manual_seed(6)
    coords = torch.randn((100, 3), generator=generator, dtype=torch.float64) * 5.0
    radii = torch.full((100,), 1.70, dtype=torch.float64)

    values = sasa(coords, radii, 1.4, 50)

    assert values.shape == (100,)
    assert bool((values >= 0).all())
    assert bool(torch.isfinite(values).all())
    # A densely packed random cloud has atoms that are completely buried, so
    # a zero area is legitimate; the total is not.
    assert float(values.sum()) > 0.0
