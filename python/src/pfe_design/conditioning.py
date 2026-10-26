"""Target-specific pocket conditioning.

The conditioner turns two things into the context that the vector field needs:

* the target protein's backbone coordinates ``X_target`` (B, Nt, 3), and
* a binary mask ``M`` (B, N) over the *complex*, where ``M[i] == 1`` marks a
  residue that is free to move (the designed binder) and ``M[i] == 0`` marks a
  residue that is fixed (the target).

Mask convention
---------------
``M`` is the **dynamic** mask.  Its complement is the fixed mask, and the
fixed target coordinates are carried along in :attr:`ConditioningOutput.fixed_xyz`
so that an integrator can re-impose them exactly after every step.

Complex layout
--------------
The target occupies the first ``Nt`` positions of the complex and the designed
region follows::

    complex index:  0 ............ Nt-1 | Nt ............ N-1
                    target (fixed)      | design (dynamic)

Pocket definition
-----------------
The pocket is the set of target residues that line the binding site.  It is
resolved in this order:

1. an explicitly supplied ``pocket_center``,
2. otherwise the target residues that lie within ``pocket_cutoff`` of the
   current design coordinates (the design region defines the site),
3. otherwise the centroid of the modelled target residues.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import torch
from torch import Tensor, nn

from .config import ModelConfig, TargetConditioningConfig

__all__ = [
    "ConditioningOutput",
    "PocketConditioner",
    "build_complex_mask",
    "dynamic_mask_from_layout",
]

# Width of the raw geometric feature vector fed to the context MLP:
#   relative position to pocket centre (3), normalised distance to centre (1),
#   local target density (1), distance to the nearest target residue (1),
#   dynamic flag (1), pocket-membership flag (1).
_FEATURE_DIM = 8


@dataclass(frozen=True)
class ConditioningOutput:
    """Everything the vector field and the integrator need after conditioning.

    Attributes
    ----------
    context:
        Per-residue context vectors, ``(B, N, feature_dim)``.
    pocket_center:
        Centre of the binding pocket, ``(B, 3)``.
    pocket_radius:
        Radius of the binding pocket, ``(B,)`` - the largest distance from the
        centre to a residue lining the pocket.
    dynamic_mask:
        ``M``: 1 for a designed residue, 0 for a fixed target residue, ``(B, N)``.
    fixed_mask:
        ``1 - M``, ``(B, N)``.
    valid_mask:
        1 for a modelled residue, ``(B, N)``.
    fixed_xyz:
        Coordinates of the fixed residues inside the complex layout; rows for
        dynamic residues are zero, ``(B, N, 3)``.
    pocket_mask:
        1 for the target residues lining the pocket, ``(B, Nt)``.
    """

    context: Tensor
    pocket_center: Tensor
    pocket_radius: Tensor
    dynamic_mask: Tensor
    fixed_mask: Tensor
    valid_mask: Tensor
    fixed_xyz: Tensor
    pocket_mask: Tensor

    @property
    def n_residues(self) -> int:
        """Number of residues in the complex."""
        return int(self.context.shape[1])

    @property
    def n_target(self) -> int:
        """Number of target residues."""
        return int(self.pocket_mask.shape[1])

    @property
    def n_design(self) -> int:
        """Number of designed residues."""
        return self.n_residues - self.n_target

    def to(self, device: torch.device | str) -> ConditioningOutput:
        """Return a copy with every tensor moved to ``device``."""
        moved = {}
        for name in (
            "context",
            "pocket_center",
            "pocket_radius",
            "dynamic_mask",
            "fixed_mask",
            "valid_mask",
            "fixed_xyz",
            "pocket_mask",
        ):
            moved[name] = getattr(self, name).to(device)
        return replace(self, **moved)


def build_complex_mask(
    n_target: int,
    n_design: int,
    batch_size: int = 1,
    device: torch.device | str = "cpu",
    dtype: torch.dtype = torch.float32,
) -> Tensor:
    """Build the binary complex mask ``M``.

    Parameters
    ----------
    n_target:
        Number of fixed target residues.
    n_design:
        Number of designed (dynamic) residues.
    batch_size:
        Number of complexes in the batch.
    device, dtype:
        Placement of the returned tensor.

    Returns
    -------
    Tensor
        ``(batch_size, n_target + n_design)`` mask with 1 on the designed
        residues (the last ``n_design`` columns) and 0 on the target residues.
    """
    if n_target < 0 or n_design < 0:
        raise ValueError("n_target and n_design must be non-negative")
    total = n_target + n_design
    mask = torch.zeros((batch_size, total), device=device, dtype=dtype)
    if n_design:
        mask[:, n_target:] = 1.0
    return mask


def dynamic_mask_from_layout(
    n_target: int,
    n_design: int,
    batch_size: int = 1,
    device: torch.device | str = "cpu",
    dtype: torch.dtype = torch.float32,
) -> Tensor:
    """Alias of :func:`build_complex_mask` kept for readability at call sites."""
    return build_complex_mask(n_target, n_design, batch_size, device, dtype)


class PocketConditioner(nn.Module):
    """Build per-residue context from the target structure and the mask ``M``.

    The conditioner is pure geometry: it never looks at the sequence, so it can
    be reused for any target and for de novo binders with no prior structure.
    """

    def __init__(
        self,
        cfg: TargetConditioningConfig | None = None,
        model_cfg: ModelConfig | None = None,
    ) -> None:
        super().__init__()
        self.cfg = cfg or TargetConditioningConfig()
        self.model_cfg = model_cfg or ModelConfig()
        self.feature_mlp = nn.Sequential(
            nn.Linear(_FEATURE_DIM, self.cfg.feature_dim),
            nn.SiLU(),
            nn.Linear(self.cfg.feature_dim, self.cfg.feature_dim),
        )
        self.dropout = nn.Dropout(self.cfg.dropout)

    # ------------------------------------------------------------------ utils

    @staticmethod
    def _masked_centroid(coords: Tensor, valid: Tensor) -> Tensor:
        """Centroid of the valid rows of ``coords`` (B, N, 3)."""
        weights = valid.to(coords.dtype).unsqueeze(-1)
        # weights.sum(dim=1) is (B, 1), which divides a (B, 3) centroid.
        total = weights.sum(dim=1).clamp(min=1.0)
        return (coords * weights).sum(dim=1) / total

    def _resolve_pocket(
        self,
        target_xyz: Tensor,
        target_valid: Tensor,
        dynamic_mask: Tensor,
        complex_xyz: Tensor | None,
        pocket_center: Tensor | None,
    ) -> tuple[Tensor, Tensor]:
        """Return ``(pocket_center, pocket_mask)`` for the complex.

        ``pocket_mask`` is ``(B, Nt)``: 1 for the target residues lining the
        pocket.
        """
        batch, n_target, _ = target_xyz.shape
        device, dtype = target_xyz.device, target_xyz.dtype
        cutoff = float(self.cfg.pocket_cutoff)

        if pocket_center is not None:
            center = pocket_center.to(device=device, dtype=dtype)
            if center.dim() == 1:
                center = center.unsqueeze(0).expand(batch, 3)
            if center.shape != (batch, 3):
                raise ValueError(
                    f"pocket_center must have shape (batch, 3); got {tuple(pocket_center.shape)}"
                )
            dist = torch.linalg.vector_norm(target_xyz - center.unsqueeze(1), dim=-1)
            mask = (dist < cutoff) & target_valid.bool()
        elif complex_xyz is not None and complex_xyz.shape[1] > n_target:
            # The design region defines the site: a target residue lines the
            # pocket when it lies within the cutoff of any designed residue.
            # The fallback for a complex whose design region is nowhere near
            # its target is applied per structure, so one member of a batch
            # cannot change the pocket of another.
            design_xyz = complex_xyz[:, n_target:]
            has_design = dynamic_mask[:, n_target:].any(dim=-1)
            target_d = torch.cdist(target_xyz, design_xyz)
            mask = ((target_d < cutoff).any(dim=-1)) & target_valid.bool()
            mask = mask & has_design.unsqueeze(-1)

            weights = mask.to(dtype).unsqueeze(-1)
            sums = weights.sum(dim=1)
            empty = sums.squeeze(-1) <= 0
            pocket_centroid = (target_xyz * weights).sum(dim=1) / sums.clamp(min=1.0)
            fallback = self._masked_centroid(target_xyz, target_valid)
            center = torch.where(empty.unsqueeze(-1), fallback, pocket_centroid)
            # A complex with no pocket residue takes every modelled target
            # residue, so the pocket is never empty.
            mask = torch.where(empty.unsqueeze(-1), target_valid.bool(), mask)
        else:
            center = self._masked_centroid(target_xyz, target_valid)
            mask = target_valid.bool()

        # Cap the pocket at the configured number of residues, closest first.
        cap = self.cfg.max_pocket_residues
        if cap is not None and n_target > cap:
            dist_to_center = torch.linalg.vector_norm(
                target_xyz - center.unsqueeze(1), dim=-1
            )
            dist_to_center = dist_to_center.masked_fill(~mask, float("inf"))
            order = dist_to_center.argsort(dim=-1)
            keep = order[:, :cap]
            trimmed = torch.zeros_like(mask)
            trimmed.scatter_(1, keep, True)
            mask = mask & trimmed

        return center, mask

    # ---------------------------------------------------------------- forward

    def forward(
        self,
        target_xyz: Tensor,
        target_valid: Tensor,
        dynamic_mask: Tensor,
        complex_xyz: Tensor | None = None,
        pocket_center: Tensor | None = None,
    ) -> ConditioningOutput:
        """Condition on the target structure.

        Parameters
        ----------
        target_xyz:
            Target backbone coordinates, ``(B, Nt, 3)``.
        target_valid:
            1 for a modelled target residue, ``(B, Nt)``.
        dynamic_mask:
            The binary mask ``M`` over the complex, ``(B, Nt + Nd)``; 1 marks a
            designed residue.
        complex_xyz:
            Optional current coordinates of the whole complex, ``(B, N, 3)``.
            When omitted, the designed residues are assumed to start at the
            pocket centre.
        pocket_center:
            Optional explicit pocket centre, ``(B, 3)``.

        Returns
        -------
        ConditioningOutput
        """
        if target_xyz.dim() != 3 or target_xyz.shape[-1] != 3:
            raise ValueError(
                f"target_xyz must have shape (batch, n_target, 3); got {tuple(target_xyz.shape)}"
            )
        batch, n_target, _ = target_xyz.shape
        if target_valid.shape != (batch, n_target):
            raise ValueError(
                f"target_valid must have shape {(batch, n_target)}; got {tuple(target_valid.shape)}"
            )
        if dynamic_mask.dim() != 2 or dynamic_mask.shape[0] != batch:
            raise ValueError(
                f"dynamic_mask must have shape (batch, n_complex); got {tuple(dynamic_mask.shape)}"
            )
        n_complex = dynamic_mask.shape[1]
        if n_complex < n_target:
            raise ValueError(
                f"the complex mask has {n_complex} columns but the target has {n_target} residues"
            )
        if bool((dynamic_mask[:, :n_target] != 0).any()):
            raise ValueError(
                "the first n_target columns of the mask are the target and must be 0 (fixed)"
            )
        # The mask has to be a clean binary partition. A fractional value would
        # be read as a weight by the velocity mask and as a flag by the
        # coordinate write-back, and the two readings disagree.
        if not bool(((dynamic_mask == 0) | (dynamic_mask == 1)).all()):
            raise ValueError(
                "the mask must be binary (1 = designed and free to move, "
                "0 = fixed target); a fractional value is read as a weight by "
                "one constraint and as a flag by the other"
            )

        device = target_xyz.device
        dtype = target_xyz.dtype
        if self.cfg.detach_target:
            target_xyz = target_xyz.detach()

        dynamic = dynamic_mask.to(dtype)
        fixed = 1.0 - dynamic
        # A designed residue is always valid; a target residue is valid only
        # when the structure actually models it.
        valid = dynamic.clone()
        valid[:, :n_target] = fixed[:, :n_target] * target_valid.to(dtype)

        center, pocket_mask = self._resolve_pocket(
            target_xyz, target_valid.to(dtype), dynamic, complex_xyz, pocket_center
        )

        # Pocket radius: the furthest residue that lines the pocket.
        dist_center = torch.linalg.vector_norm(target_xyz - center.unsqueeze(1), dim=-1)
        pocket_radius = (dist_center * pocket_mask.to(dtype)).max(dim=1).values
        pocket_radius = pocket_radius.clamp(min=1.0)
        if not self.cfg.normalize_features:
            pocket_radius = torch.ones_like(pocket_radius)

        # Coordinates of the whole complex; the design region defaults to the
        # pocket centre when no complex coordinates are supplied.
        if complex_xyz is None:
            complex_xyz = target_xyz.new_zeros((batch, n_complex, 3))
            complex_xyz = complex_xyz + center.unsqueeze(1)
            complex_xyz[:, :n_target] = target_xyz
        if complex_xyz.shape != (batch, n_complex, 3):
            raise ValueError(
                f"complex_xyz must have shape {(batch, n_complex, 3)}; "
                f"got {tuple(complex_xyz.shape)}"
            )

        # Distances from every complex residue to every target residue.
        pair = torch.cdist(complex_xyz, target_xyz)  # (B, N, Nt)
        target_valid_f = target_valid.to(dtype)
        pair = pair.masked_fill(
            (~target_valid_f.bool()).unsqueeze(1), float("inf")
        )

        cutoff = float(self.cfg.pocket_cutoff)
        density = ((pair < cutoff).to(dtype).sum(dim=-1)) / max(1, n_target)

        nearest = pair.min(dim=-1).values
        nearest = torch.where(fixed > 0, torch.zeros_like(nearest), nearest)
        nearest = torch.nan_to_num(nearest, posinf=0.0, neginf=0.0)

        rel = (complex_xyz - center.unsqueeze(1)) / pocket_radius[:, None, None]
        dist_norm = torch.linalg.vector_norm(rel, dim=-1, keepdim=True)

        pocket_flag = torch.zeros((batch, n_complex), device=device, dtype=dtype)
        if self.cfg.use_pocket_membership:
            pocket_flag[:, :n_target] = pocket_mask.to(dtype)

        features = torch.cat(
            [
                rel,
                dist_norm,
                density.unsqueeze(-1),
                (nearest / pocket_radius[:, None]).unsqueeze(-1),
                dynamic.unsqueeze(-1),
                pocket_flag.unsqueeze(-1),
            ],
            dim=-1,
        )
        context = self.dropout(self.feature_mlp(features))

        fixed_xyz = target_xyz.new_zeros((batch, n_complex, 3))
        fixed_xyz[:, :n_target] = target_xyz

        return ConditioningOutput(
            context=context,
            pocket_center=center,
            pocket_radius=pocket_radius,
            dynamic_mask=dynamic,
            fixed_mask=fixed,
            valid_mask=valid,
            fixed_xyz=fixed_xyz,
            pocket_mask=pocket_mask.to(dtype),
        )
