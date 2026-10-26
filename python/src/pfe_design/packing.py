"""Refinement and side-chain packing hook.

The hook is the bridge between a generated backbone and a real structure: it
takes the coordinates the model produced, prepares them for packing, relaxes
the remaining steric clashes while the fixed target residues stay exactly where
they were, and can then hand the result to an external packer.

External packers are integrated as *stubs with a real seam*: each backend
declares whether it is installed, and a backend can be replaced by a plain
callable so that a site-specific installation can be dropped in without
touching this module.

* ``MockPackingBackend`` - deterministic, dependency-free relaxation.
* ``RosettaBackend`` - PyRosetta fast-relax / fixbb interface.
* ``OpenMMBackend`` - energy minimisation with an OpenMM system.
* ``ESMFoldBackend`` - side-chain reconstruction from a predicted structure.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field, replace

import torch
from torch import Tensor, nn

from .config import PackingBackend, PackingConfig
from .metrics import clash_report

__all__ = [
    "BackendUnavailableError",
    "PackingRequest",
    "PackingResult",
    "PackingBackendBase",
    "MockPackingBackend",
    "RosettaBackend",
    "OpenMMBackend",
    "ESMFoldBackend",
    "SideChainPackingHook",
    "relax_clashes",
    "build_backend",
]


class BackendUnavailableError(RuntimeError):
    """Raised when a packing backend is asked to run but is not installed."""


@dataclass(frozen=True)
class PackingRequest:
    """Everything a backend needs to pack a structure.

    Attributes
    ----------
    backbone_xyz:
        Backbone coordinates of the complex, ``(N, 3)``.
    residue_types:
        Residue-type index per residue, ``(N,)``.
    fixed_mask:
        1 for the residues that must not move (the target), ``(N,)``.
    config:
        The packing settings.
    """

    backbone_xyz: Tensor
    residue_types: Tensor
    fixed_mask: Tensor
    config: PackingConfig
    metadata: dict = field(default_factory=dict)

    @property
    def n_residues(self) -> int:
        return int(self.backbone_xyz.shape[0])


@dataclass(frozen=True)
class PackingResult:
    """The outcome of the refinement and packing hook."""

    coordinates: Tensor
    """Final coordinates, ``(N, 3)``."""

    initial_clash_penalty: float
    """Clash penalty before refinement."""

    clash_penalty: float
    """Clash penalty after refinement (and after the backend, if accepted)."""

    iterations: int
    """Number of relaxation iterations that were run."""

    backend_used: str | None = None
    """Name of the backend that produced the final coordinates, if any."""

    accepted_backend: bool = False
    """Whether the backend result was kept (it may have been rejected)."""

    @property
    def improved(self) -> bool:
        """True when the hook reduced the clash penalty."""
        return self.clash_penalty < self.initial_clash_penalty


_FALLBACK_DIRECTIONS = torch.tensor(
    [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0]]
)


def _unit_directions(delta: Tensor, distance: Tensor) -> Tensor:
    """Unit vectors pointing from atom j to atom i.

    Two atoms at exactly the same position have no direction between them, so
    the repulsive displacement would be zero and they could never be
    separated. For those pairs a deterministic direction derived from the atom
    indices is used instead; it is antisymmetric, so the two atoms are pushed
    apart rather than both staying put.
    """
    safe = distance.clamp(min=1e-9)
    direction = delta / safe[..., None]
    coincident = distance < 1e-8
    if bool(coincident.any()):
        n = delta.shape[0]
        idx = torch.arange(n, device=delta.device)
        lo = torch.minimum(idx[:, None], idx[None, :])
        hi = torch.maximum(idx[:, None], idx[None, :])
        choice = ((lo * 73856093) ^ (hi * 19349663)) % 6
        fallback = _FALLBACK_DIRECTIONS.to(device=delta.device, dtype=delta.dtype)[choice]
        sign = torch.where(idx[:, None] < idx[None, :], 1.0, -1.0).to(delta.dtype)
        direction = torch.where(
            coincident[..., None], fallback * sign[..., None], direction
        )
    return direction


def _relax_clashes(
    coords: Tensor,
    fixed_mask: Tensor | None = None,
    clash_distance: float = 2.0,
    max_iterations: int = 100,
    step_size: float = 0.05,
    residue_index: Tensor | None = None,
) -> tuple[Tensor, int]:
    """The relaxation loop, also reporting how many sweeps it actually ran."""
    x = coords.to(torch.float64).clone()
    n = x.shape[0]
    if n < 2:
        return x.to(coords.dtype), 0

    if fixed_mask is None:
        fixed = torch.zeros(n, dtype=torch.bool, device=x.device)
    else:
        fixed = torch.as_tensor(fixed_mask, device=x.device).reshape(-1).to(torch.bool)

    delta = x[:, None, :] - x[None, :, :]
    distance = torch.linalg.vector_norm(delta, dim=-1)
    direction = _unit_directions(delta, distance)

    overlap = (float(clash_distance) - distance).clamp(min=0.0)
    eye = torch.eye(n, dtype=torch.bool, device=x.device)
    overlap = overlap.masked_fill(eye, 0.0)
    same = None
    if residue_index is not None:
        res = torch.as_tensor(residue_index, device=x.device).reshape(-1).to(torch.long)
        same = res[:, None] == res[None, :]
        overlap = overlap.masked_fill(same, 0.0)

    iterations = 0
    for iteration in range(int(max_iterations)):
        if float(overlap.max()) <= 0.0:
            break
        iterations = iteration + 1
        push = 0.5 * float(step_size) * overlap[..., None] * direction
        movement = push.sum(dim=1)
        movement[fixed] = 0.0
        x = x + movement
        delta = x[:, None, :] - x[None, :, :]
        distance = torch.linalg.vector_norm(delta, dim=-1)
        direction = _unit_directions(delta, distance)
        overlap = (float(clash_distance) - distance).clamp(min=0.0)
        overlap = overlap.masked_fill(eye, 0.0)
        if same is not None:
            overlap = overlap.masked_fill(same, 0.0)

    return x.to(coords.dtype), iterations


def relax_clashes(
    coords: Tensor,
    fixed_mask: Tensor | None = None,
    clash_distance: float = 2.0,
    max_iterations: int = 100,
    step_size: float = 0.05,
    residue_index: Tensor | None = None,
) -> Tensor:
    """Push overlapping atoms apart, leaving the fixed residues untouched.

    The relaxation is a deterministic repulsive displacement: every pair that
    is closer than ``clash_distance`` contributes a displacement along the
    line joining the two atoms.  Fixed residues never receive a displacement,
    so they come out bit-identical to the input.

    Parameters
    ----------
    coords:
        Coordinates, ``(N, 3)``.
    fixed_mask:
        1 for the residues that must not move.
    clash_distance:
        Minimum allowed distance between two atoms.
    max_iterations:
        Maximum number of sweeps.
    step_size:
        Fraction of the overlap that is corrected per sweep.
    residue_index:
        Optional residue index per atom; pairs inside one residue are ignored.

    Returns
    -------
    Tensor
        Relaxed coordinates with the same shape as ``coords``.
    """
    relaxed, _iterations = _relax_clashes(
        coords,
        fixed_mask=fixed_mask,
        clash_distance=clash_distance,
        max_iterations=max_iterations,
        step_size=step_size,
        residue_index=residue_index,
    )
    return relaxed


class PackingBackendBase:
    """Base class for the external packers."""

    name = "base"

    @classmethod
    def available(cls) -> bool:
        """Whether the backend can run in this environment."""
        return False

    def pack(self, request: PackingRequest) -> Tensor:
        """Return packed coordinates for ``request``."""
        raise BackendUnavailableError(
            f"the {self.name} backend is a stub and has no packer installed"
        )


class MockPackingBackend(PackingBackendBase):
    """Dependency-free deterministic stand-in for a real packer."""

    name = "mock"

    @classmethod
    def available(cls) -> bool:
        return True

    def pack(self, request: PackingRequest) -> Tensor:
        cfg = request.config
        return relax_clashes(
            request.backbone_xyz,
            fixed_mask=request.fixed_mask,
            clash_distance=cfg.clash_distance,
            max_iterations=cfg.max_iterations,
            step_size=cfg.step_size,
        )


class RosettaBackend(PackingBackendBase):
    """PyRosetta integration stub.

    Install PyRosetta and hand the backend a callable that runs
    ``fastrelax`` or ``fixbb``; the callable receives the coordinates and the
    residue-type indices and returns new coordinates.
    """

    name = "rosetta"

    def __init__(self, packer: Callable[[PackingRequest], Tensor] | None = None) -> None:
        self._packer = packer

    @classmethod
    def available(cls) -> bool:
        try:  # pragma: no cover - depends on the environment
            import pyrosetta  # noqa: F401
        except Exception:
            return False
        return True

    def pack(self, request: PackingRequest) -> Tensor:
        if self._packer is None:
            raise BackendUnavailableError(
                "PyRosetta is not wired up: install it and pass a packer callable to "
                "RosettaBackend(packer=...), e.g. one that runs fastrelax or fixbb."
            )
        return self._packer(request)


class OpenMMBackend(PackingBackendBase):
    """OpenMM integration stub (energy minimisation of the built system)."""

    name = "openmm"

    def __init__(self, minimiser: Callable[[PackingRequest], Tensor] | None = None) -> None:
        self._minimiser = minimiser

    @classmethod
    def available(cls) -> bool:
        try:  # pragma: no cover - depends on the environment
            import openmm  # noqa: F401
        except Exception:
            return False
        return True

    def pack(self, request: PackingRequest) -> Tensor:
        if self._minimiser is None:
            raise BackendUnavailableError(
                "OpenMM is not wired up: install it and pass a minimiser callable to "
                "OpenMMBackend(minimiser=...)."
            )
        return self._minimiser(request)


class ESMFoldBackend(PackingBackendBase):
    """ESMFold integration stub.

    ESMFold predicts a structure from a sequence, so it is used here to
    rebuild side chains for a designed region: a predictor callable receives
    the backbone and returns coordinates.
    """

    name = "esmfold"

    def __init__(self, predictor: Callable[[PackingRequest], Tensor] | None = None) -> None:
        self._predictor = predictor

    @classmethod
    def available(cls) -> bool:
        return False

    def pack(self, request: PackingRequest) -> Tensor:
        if self._predictor is None:
            raise BackendUnavailableError(
                "ESMFold is not wired up: pass a predictor callable to "
                "ESMFoldBackend(predictor=...) that runs the model and returns "
                "coordinates."
            )
        return self._predictor(request)


def build_backend(backend: PackingBackend | str) -> PackingBackendBase | None:
    """Instantiate the backend named by ``backend`` (``None`` for ``none``)."""
    if isinstance(backend, str):
        backend = PackingBackend(backend)
    if backend is PackingBackend.NONE:
        return None
    if backend is PackingBackend.MOCK:
        return MockPackingBackend()
    if backend is PackingBackend.ROSETTA:
        return RosettaBackend()
    if backend is PackingBackend.OPENMM:
        return OpenMMBackend()
    if backend is PackingBackend.ESMFOLD:
        return ESMFoldBackend()
    raise ValueError(f"unknown packing backend: {backend!r}")


class SideChainPackingHook(nn.Module):
    """Modular post-processing: refine the generated backbone, then pack it.

    The hook is deliberately split into three seams:

    * :meth:`prepare` - turn generated coordinates into a
      :class:`PackingRequest`.
    * :meth:`refine` - the built-in, dependency-free clash relaxation.
    * :meth:`forward` - refine, then optionally call the backend and keep the
      better of the two results.

    Example
    -------
    >>> hook = SideChainPackingHook(PackingConfig(backend=PackingBackend.MOCK))
    >>> request = hook.prepare(backbone, residue_types, fixed_mask)
    >>> result = hook(request)
    >>> result.improved
    True
    """

    def __init__(
        self,
        cfg: PackingConfig | None = None,
        backend: PackingBackendBase | None = None,
    ) -> None:
        super().__init__()
        self.cfg = cfg or PackingConfig()
        self.backend = backend if backend is not None else build_backend(self.cfg.backend)

    # ---------------------------------------------------------------- prepare

    def prepare(
        self,
        backbone_xyz: Tensor,
        residue_types: Tensor | None = None,
        fixed_mask: Tensor | None = None,
        metadata: dict | None = None,
    ) -> PackingRequest:
        """Package generated coordinates for packing.

        Parameters
        ----------
        backbone_xyz:
            Generated backbone coordinates, ``(N, 3)``.
        residue_types:
            Residue-type index per residue, ``(N,)``; zeros when unknown.
        fixed_mask:
            1 for the residues that must not move, ``(N,)``.
        metadata:
            Free-form context passed through to the backend.
        """
        coords = backbone_xyz.to(torch.float64)
        if coords.dim() != 2 or coords.shape[1] != 3:
            raise ValueError(
                f"backbone_xyz must have shape (n_residues, 3); got {tuple(backbone_xyz.shape)}"
            )
        n = coords.shape[0]
        if residue_types is None:
            residue_types = torch.zeros(n, dtype=torch.long)
        residue_types = torch.as_tensor(residue_types).reshape(-1).to(torch.long)
        if residue_types.shape[0] != n:
            raise ValueError("residue_types must have one entry per residue")
        if fixed_mask is None:
            fixed_mask = torch.zeros(n, dtype=torch.bool)
        fixed_mask = torch.as_tensor(fixed_mask).reshape(-1).to(torch.bool)
        if fixed_mask.shape[0] != n:
            raise ValueError("fixed_mask must have one entry per residue")
        return PackingRequest(
            backbone_xyz=coords,
            residue_types=residue_types,
            fixed_mask=fixed_mask,
            config=self.cfg,
            metadata=dict(metadata or {}),
        )

    # ----------------------------------------------------------------- refine

    def refine(self, request: PackingRequest) -> tuple[Tensor, int]:
        """Run the built-in clash relaxation.

        Returns
        -------
        (coordinates, iterations)
        """
        coords = request.backbone_xyz
        fixed = request.fixed_mask if self.cfg.keep_fixed_residues else None
        before = clash_report(
            coords,
            _metric_config(self.cfg),
            residue_index=None,
        ).penalty
        relaxed, iterations = _relax_clashes(
            coords,
            fixed_mask=fixed,
            clash_distance=self.cfg.clash_distance,
            max_iterations=self.cfg.max_iterations,
            step_size=self.cfg.step_size,
        )
        after = clash_report(relaxed, _metric_config(self.cfg)).penalty
        if before <= after:
            iterations = 0
        return relaxed, iterations

    # ---------------------------------------------------------------- forward

    def forward(self, request: PackingRequest) -> PackingResult:
        """Refine, pack, and keep the better result."""
        metric_cfg = _metric_config(self.cfg)
        initial = clash_report(request.backbone_xyz, metric_cfg).penalty

        refined, iterations = self.refine(request)
        refined_penalty = clash_report(refined, metric_cfg).penalty

        # Only keep the refinement if it actually helped.
        if initial - refined_penalty < self.cfg.min_improvement:
            refined = request.backbone_xyz
            refined_penalty = initial
            iterations = 0

        backend_used: str | None = None
        accepted = False
        final = refined
        final_penalty = refined_penalty

        if self.backend is not None and self.backend.available():
            backend_request = replace(request, backbone_xyz=refined)
            try:
                packed = self.backend.pack(backend_request)
            except BackendUnavailableError:
                packed = None
            if packed is not None:
                packed = packed.to(torch.float64)
                packed_penalty = clash_report(packed, metric_cfg).penalty
                backend_used = self.backend.name
                if (not self.cfg.reject_worse) or packed_penalty <= final_penalty:
                    final = packed
                    final_penalty = packed_penalty
                    accepted = True

        return PackingResult(
            coordinates=final,
            initial_clash_penalty=initial,
            clash_penalty=final_penalty,
            iterations=iterations,
            backend_used=backend_used,
            accepted_backend=accepted,
        )


def _metric_config(cfg: PackingConfig):
    """Build the metric configuration implied by a packing configuration."""
    from .config import MetricConfig

    return MetricConfig(clash_distance=cfg.clash_distance)
