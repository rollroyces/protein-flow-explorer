"""Tests for :mod:`pfe_design.packing`.

The packing module is the seam between a generated backbone and a real
structure.  The properties that matter are physical: overlapping atoms must
end up further apart, residues that are declared fixed must come back
bit-identical, and the hook must never hand back a structure that is worse
than the one it was given.  The external backends are stubs, so what is
tested there is the seam: that they refuse to run without an installed
packer, and that they really do call a callable that was injected.
"""

from __future__ import annotations

import pytest
import torch
from torch import Tensor

from pfe_design.config import MetricConfig, PackingBackend, PackingConfig
from pfe_design.metrics import clash_report
from pfe_design.packing import (
    BackendUnavailableError,
    ESMFoldBackend,
    MockPackingBackend,
    OpenMMBackend,
    PackingBackendBase,
    PackingRequest,
    PackingResult,
    RosettaBackend,
    SideChainPackingHook,
    build_backend,
    relax_clashes,
)

# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------


def _clashing_structure(n: int = 8) -> Tensor:
    """Well separated atoms, two of which deliberately overlap."""
    coords = torch.zeros(n, 3, dtype=torch.float64)
    for i in range(n):
        coords[i, 0] = 15.0 * i
    coords[0] = torch.tensor([0.0, 0.0, 0.0], dtype=torch.float64)
    coords[1] = torch.tensor([0.3, 0.0, 0.0], dtype=torch.float64)
    return coords


def _penalty(coords: Tensor) -> float:
    return clash_report(coords, MetricConfig(clash_distance=2.0)).penalty


class _CollapsingBackend(PackingBackendBase):
    """A backend that answers with deliberately terrible coordinates."""

    name = "collapsing"

    @classmethod
    def available(cls) -> bool:
        return True

    def pack(self, request: PackingRequest) -> Tensor:
        n = request.backbone_xyz.shape[0]
        offsets = 0.01 * torch.arange(n, dtype=torch.float64).reshape(-1, 1)
        return offsets * torch.tensor([[1.0, 0.0, 0.0]], dtype=torch.float64)


# ---------------------------------------------------------------------------
# relax_clashes
# ---------------------------------------------------------------------------


def test_relax_clashes_pushes_overlapping_atoms_to_the_threshold():
    """Two atoms 0.5 A apart end up at (almost) the requested distance."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=torch.float64)

    relaxed = relax_clashes(coords, clash_distance=2.0)

    distance = float(torch.linalg.vector_norm(relaxed[0] - relaxed[1]))
    assert distance >= 2.0 - 0.05
    assert distance > 1.5  # it really moved


def test_relax_clashes_separates_a_cluster_of_overlapping_atoms():
    """A whole cluster is pushed apart, not just the worst pair."""
    torch.manual_seed(0)
    coords = torch.randn(6, 3, dtype=torch.float64) * 0.5

    relaxed = relax_clashes(coords, clash_distance=2.0, max_iterations=300)

    distances = torch.cdist(relaxed, relaxed)
    off_diagonal = distances[~torch.eye(6, dtype=torch.bool)]
    assert float(off_diagonal.min()) >= 2.0 - 0.01


def test_relax_clashes_leaves_non_overlapping_atoms_unchanged():
    """Atoms that do not clash receive no displacement at all."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [12.0, 0.0, 0.0]], dtype=torch.float64)

    relaxed = relax_clashes(coords, clash_distance=2.0)

    assert torch.equal(relaxed, coords)


def test_relax_clashes_keeps_fixed_residues_bit_identical():
    """A fixed residue must come back bit-identical, not merely close."""
    torch.manual_seed(0)
    coords = torch.tensor(
        [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [0.0, 5.0, 0.0]], dtype=torch.float64
    )
    fixed = torch.tensor([False, True, False])

    relaxed = relax_clashes(coords, fixed_mask=fixed, clash_distance=2.0)

    assert torch.equal(relaxed[1], coords[1])
    assert torch.equal(relaxed[2], coords[2])
    assert not torch.equal(relaxed[0], coords[0])


def test_relax_clashes_ignores_intra_residue_pairs():
    """Atoms of one residue are bonded, so they must not be pushed apart."""
    torch.manual_seed(0)
    coords = torch.tensor([[0.0, 0.0, 0.0], [0.5, 0.0, 0.0]], dtype=torch.float64)

    same_residue = relax_clashes(coords, clash_distance=2.0, residue_index=torch.tensor([3, 3]))
    different_residue = relax_clashes(
        coords, clash_distance=2.0, residue_index=torch.tensor([3, 4])
    )

    assert torch.equal(same_residue, coords)
    assert not torch.equal(different_residue, coords)


def test_relax_clashes_removes_the_clash_penalty():
    """After relaxation no pair is closer than the threshold."""
    torch.manual_seed(0)
    coords = _clashing_structure()

    relaxed = relax_clashes(coords, clash_distance=2.0, max_iterations=300)

    assert _penalty(coords) > 1.0
    assert _penalty(relaxed) == pytest.approx(0.0, abs=1e-6)


def test_relax_clashes_is_deterministic():
    """The relaxation is a fixed rule, so the same input gives the same output."""
    torch.manual_seed(0)
    coords = _clashing_structure()

    first = relax_clashes(coords, clash_distance=2.0)
    second = relax_clashes(coords, clash_distance=2.0)

    assert torch.equal(first, second)


def test_relax_clashes_preserves_shape_and_dtype():
    """The relaxed coordinates have the same layout as the input."""
    torch.manual_seed(0)
    coords = _clashing_structure()

    relaxed = relax_clashes(coords, clash_distance=2.0)

    assert relaxed.shape == coords.shape
    assert relaxed.dtype == coords.dtype


# ---------------------------------------------------------------------------
# PackingRequest / PackingResult
# ---------------------------------------------------------------------------


def test_packing_request_reports_its_residue_count():
    """The request knows how many residues it carries."""
    torch.manual_seed(0)
    coords = _clashing_structure(6)
    request = PackingRequest(
        backbone_xyz=coords,
        residue_types=torch.zeros(6, dtype=torch.long),
        fixed_mask=torch.zeros(6, dtype=torch.bool),
        config=PackingConfig(),
    )

    assert request.n_residues == 6


def test_packing_request_carries_its_configuration_and_metadata():
    """The settings and the free-form context travel with the request."""
    torch.manual_seed(0)
    coords = _clashing_structure(4)
    cfg = PackingConfig(clash_distance=2.5)
    request = PackingRequest(
        backbone_xyz=coords,
        residue_types=torch.zeros(4, dtype=torch.long),
        fixed_mask=torch.zeros(4, dtype=torch.bool),
        config=cfg,
        metadata={"seed": 3},
    )

    assert request.config is cfg
    assert request.metadata["seed"] == 3


def test_packing_result_improved_is_true_when_the_penalty_falls():
    """A lower clash penalty counts as an improvement."""
    result = PackingResult(
        coordinates=torch.zeros(2, 3),
        initial_clash_penalty=2.0,
        clash_penalty=1.0,
        iterations=10,
    )

    assert result.improved is True


def test_packing_result_improved_is_false_when_the_penalty_is_unchanged():
    """No change is not an improvement."""
    result = PackingResult(
        coordinates=torch.zeros(2, 3),
        initial_clash_penalty=1.0,
        clash_penalty=1.0,
        iterations=0,
    )

    assert result.improved is False


def test_packing_result_improved_is_false_when_the_penalty_rises():
    """A worse structure is never reported as an improvement."""
    result = PackingResult(
        coordinates=torch.zeros(2, 3),
        initial_clash_penalty=1.0,
        clash_penalty=3.0,
        iterations=0,
    )

    assert result.improved is False


# ---------------------------------------------------------------------------
# SideChainPackingHook.prepare
# ---------------------------------------------------------------------------


def test_prepare_rejects_coordinates_that_are_not_n_by_three():
    """The coordinate tensor must be (n_residues, 3)."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    with pytest.raises(ValueError):
        hook.prepare(torch.zeros(5, 2, dtype=torch.float64))


def test_prepare_rejects_residue_types_that_do_not_match_the_residue_count():
    """One residue type per residue, no more and no fewer."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    coords = torch.zeros(5, 3, dtype=torch.float64)
    with pytest.raises(ValueError):
        hook.prepare(coords, torch.zeros(4, dtype=torch.long))


def test_prepare_rejects_a_fixed_mask_that_does_not_match_the_residue_count():
    """The fixed mask must have one entry per residue."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    coords = torch.zeros(5, 3, dtype=torch.float64)
    with pytest.raises(ValueError):
        hook.prepare(coords, torch.zeros(5, dtype=torch.long), torch.zeros(3, dtype=torch.bool))


def test_prepare_returns_a_request_that_carries_the_metadata():
    """The context handed to prepare reaches the backend unchanged."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    coords = _clashing_structure(5)

    request = hook.prepare(
        coords,
        torch.zeros(5, dtype=torch.long),
        torch.zeros(5, dtype=torch.bool),
        metadata={"target": "1ABC", "cycle": 7},
    )

    assert isinstance(request, PackingRequest)
    assert request.n_residues == 5
    assert request.metadata == {"target": "1ABC", "cycle": 7}
    assert request.config is hook.cfg
    assert torch.allclose(request.backbone_xyz, coords)


def test_prepare_defaults_the_residue_types_to_zeros():
    """Unknown residue types are encoded as zero."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    coords = _clashing_structure(4)

    request = hook.prepare(coords)

    assert request.residue_types.tolist() == [0, 0, 0, 0]
    assert request.fixed_mask.tolist() == [False] * 4


# ---------------------------------------------------------------------------
# SideChainPackingHook.refine / forward
# ---------------------------------------------------------------------------


def test_refine_removes_the_clash_penalty():
    """The built-in relaxation resolves the overlaps of the input."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    coords = _clashing_structure()
    request = hook.prepare(coords)

    refined, iterations = hook.refine(request)

    assert refined.shape == coords.shape
    assert iterations >= 0
    assert _penalty(refined) < _penalty(coords)
    # No pair may stay meaningfully closer than the clash threshold.
    distances = torch.cdist(refined, refined)
    off_diagonal = distances[~torch.eye(refined.shape[0], dtype=torch.bool)]
    assert float(off_diagonal.min()) >= 2.0 - 0.05
    assert _penalty(refined) < 1e-3 * _penalty(coords)


def test_forward_never_worsens_the_clash_penalty():
    """The hook must not hand back a structure that clashes more."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.MOCK))
    request = hook.prepare(_clashing_structure())

    result = hook(request)

    assert result.clash_penalty <= result.initial_clash_penalty + 1e-12
    assert result.improved is True


def test_forward_keeps_fixed_residues_unchanged():
    """A fixed residue that sits in a clash must still not move."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.MOCK))
    coords = _clashing_structure()
    fixed = torch.zeros(coords.shape[0], dtype=torch.bool)
    fixed[1] = True
    request = hook.prepare(coords, fixed_mask=fixed)

    result = hook(request)

    assert torch.equal(result.coordinates[1], coords[1])
    assert torch.equal(result.coordinates[5], coords[5])


def test_forward_reports_the_mock_backend_and_accepts_its_result():
    """The mock backend is available, so it runs and its result is kept."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.MOCK))
    request = hook.prepare(_clashing_structure())

    result = hook(request)

    assert result.backend_used == "mock"
    assert result.accepted_backend is True


def test_forward_without_a_backend_does_not_use_one():
    """With no backend configured the hook stops after the refinement."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure())

    result = hook(request)

    assert result.backend_used is None
    assert result.accepted_backend is False
    assert result.clash_penalty <= result.initial_clash_penalty


def test_forward_rejects_a_worse_backend_result():
    """A backend that makes things worse is discarded."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(
        PackingConfig(backend=PackingBackend.MOCK), backend=_CollapsingBackend()
    )
    request = hook.prepare(_clashing_structure())
    refined, _ = hook.refine(request)

    result = hook(request)

    assert result.backend_used == "collapsing"
    assert result.accepted_backend is False
    assert result.clash_penalty <= result.initial_clash_penalty
    assert torch.allclose(result.coordinates, refined, atol=1e-12)


def test_forward_accepts_a_worse_backend_result_when_reject_worse_is_off():
    """With reject_worse disabled the backend answer is taken as it is."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(
        PackingConfig(backend=PackingBackend.MOCK, reject_worse=False),
        backend=_CollapsingBackend(),
    )
    request = hook.prepare(_clashing_structure())

    result = hook(request)

    assert result.accepted_backend is True
    assert result.backend_used == "collapsing"
    assert result.clash_penalty > result.initial_clash_penalty


def test_refine_moves_fixed_residues_when_keep_fixed_is_disabled():
    """The fixed-residue guarantee is a configuration choice, not a hard rule."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(
        PackingConfig(backend=PackingBackend.NONE, keep_fixed_residues=False)
    )
    coords = _clashing_structure()
    fixed = torch.zeros(coords.shape[0], dtype=torch.bool)
    fixed[1] = True
    request = hook.prepare(coords, fixed_mask=fixed)

    refined, _ = hook.refine(request)

    assert not torch.equal(refined[1], coords[1])


# ---------------------------------------------------------------------------
# backends
# ---------------------------------------------------------------------------


def test_rosetta_backend_refuses_to_run_without_pyosetta():
    """The PyRosetta stub has no packer installed, so it must refuse."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    backend = RosettaBackend()
    with pytest.raises(BackendUnavailableError):
        backend.pack(request)


def test_openmm_backend_refuses_to_run_without_openmm():
    """The OpenMM stub has no minimiser installed, so it must refuse."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    backend = OpenMMBackend()
    with pytest.raises(BackendUnavailableError):
        backend.pack(request)


def test_esmfold_backend_refuses_to_run_without_a_predictor():
    """The ESMFold stub has no predictor installed, so it must refuse."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    backend = ESMFoldBackend()
    with pytest.raises(BackendUnavailableError):
        backend.pack(request)


def test_backend_stubs_report_themselves_as_unavailable():
    """None of the external backers are installed in this environment."""
    assert RosettaBackend.available() is False
    assert OpenMMBackend.available() is False
    assert ESMFoldBackend.available() is False


def test_base_backend_refuses_to_run():
    """The base class has no packer at all."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    with pytest.raises(BackendUnavailableError):
        PackingBackendBase().pack(request)


def test_build_backend_returns_none_for_the_none_backend():
    """'none' means no backend at all."""
    assert build_backend("none") is None
    assert build_backend(PackingBackend.NONE) is None


def test_build_backend_returns_a_mock_backend_for_the_mock_name():
    """'mock' maps to the dependency-free stand-in."""
    backend = build_backend("mock")
    assert isinstance(backend, MockPackingBackend)
    assert backend.name == "mock"
    assert MockPackingBackend.available() is True


def test_build_backend_rejects_an_unknown_backend_name():
    """An unknown name is a configuration error."""
    with pytest.raises(ValueError):
        build_backend("not-a-backend")


def test_mock_backend_relaxes_the_clashes_of_a_request():
    """The mock backend is not a no-op: it resolves the overlaps."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure())

    packed = MockPackingBackend().pack(request)

    assert packed.shape == request.backbone_xyz.shape
    assert _penalty(packed) < _penalty(request.backbone_xyz)


def test_rosetta_backend_uses_an_injected_packer():
    """The injected callable is really the thing that produces the coordinates."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    seen: list[PackingRequest] = []

    def packer(req: PackingRequest) -> Tensor:
        seen.append(req)
        return req.backbone_xyz + 1.0

    packed = RosettaBackend(packer=packer).pack(request)

    assert len(seen) == 1
    assert seen[0] is request
    assert torch.allclose(packed, request.backbone_xyz + 1.0)


def test_openmm_backend_uses_an_injected_minimiser():
    """The injected callable is really the thing that produces the coordinates."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    seen: list[PackingRequest] = []

    def minimiser(req: PackingRequest) -> Tensor:
        seen.append(req)
        return req.backbone_xyz * 2.0

    packed = OpenMMBackend(minimiser=minimiser).pack(request)

    assert len(seen) == 1
    assert seen[0] is request
    assert torch.allclose(packed, request.backbone_xyz * 2.0)


def test_esmfold_backend_uses_an_injected_predictor():
    """The injected callable is really the thing that produces the coordinates."""
    torch.manual_seed(0)
    hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.NONE))
    request = hook.prepare(_clashing_structure(4))
    seen: list[PackingRequest] = []

    def predictor(req: PackingRequest) -> Tensor:
        seen.append(req)
        return req.backbone_xyz - 2.0

    packed = ESMFoldBackend(predictor=predictor).pack(request)

    assert len(seen) == 1
    assert seen[0] is request
    assert torch.allclose(packed, request.backbone_xyz - 2.0)


def test_a_backend_that_returns_nan_is_discarded():
    """A diverged minimisation must not be scored as a perfect zero.

    A NaN comparison is false, so a clash score computed on NaN coordinates
    comes out as 0.0 - "perfectly clash-free" - and a naive hook would
    accept a structure that has fallen apart.
    """
    import torch

    from pfe_design.packing import PackingBackendBase

    class NanBackend(PackingBackendBase):
        name = "nan"

        @classmethod
        def available(cls) -> bool:
            return True

        def pack(self, request):
            return torch.full_like(request.backbone_xyz, float("nan"))

    coords = torch.tensor([[0.0, 0.0, 0.0], [0.0, 0.0, 1.9]], dtype=torch.float64)
    hook = SideChainPackingHook(PackingConfig(), backend=NanBackend())
    request = hook.prepare(coords, torch.zeros(2, dtype=torch.long))

    result = hook(request)

    assert result.accepted_backend is False
    assert result.rejection is not None and "non-finite" in result.rejection
    assert bool(torch.isfinite(result.coordinates).all())
    # The coordinates are the refined ones, not the backend's NaNs.
    assert bool(torch.isfinite(result.coordinates).all())


def test_a_backend_that_returns_the_wrong_shape_is_discarded():
    """A backend that returns something else is not scored at all."""
    import torch

    from pfe_design.packing import PackingBackendBase

    class WrongShapeBackend(PackingBackendBase):
        name = "wrong-shape"

        @classmethod
        def available(cls) -> bool:
            return True

        def pack(self, request):
            return request.backbone_xyz[:1]

    coords = torch.tensor([[0.0, 0.0, 0.0], [0.0, 0.0, 1.9]], dtype=torch.float64)
    hook = SideChainPackingHook(PackingConfig(), backend=WrongShapeBackend())
    request = hook.prepare(coords, torch.zeros(2, dtype=torch.long))

    result = hook(request)

    assert result.accepted_backend is False
    assert result.rejection is not None and "expected" in result.rejection


def test_a_backend_cannot_move_the_fixed_residues():
    """The target is rigid: what comes out of the hook keeps it in place."""
    import torch

    from pfe_design.packing import PackingBackendBase

    class DriftingBackend(PackingBackendBase):
        name = "drifting"

        @classmethod
        def available(cls) -> bool:
            return True

        def pack(self, request):
            moved = request.backbone_xyz + 5.0
            return moved

    coords = torch.tensor([[0.0, 0.0, 0.0], [0.0, 0.0, 1.9], [0.0, 0.0, 4.0]],
                           dtype=torch.float64)
    fixed = torch.tensor([True, False, False])
    hook = SideChainPackingHook(PackingConfig(), backend=DriftingBackend())
    request = hook.prepare(coords, torch.zeros(3, dtype=torch.long), fixed)

    result = hook(request)

    assert result.accepted_backend is True
    assert torch.equal(result.coordinates[0], coords[0])
