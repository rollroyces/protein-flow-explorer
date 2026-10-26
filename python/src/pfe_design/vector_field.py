"""The backbone vector field, conditioned on a target pocket.

This is the network the objective calls "the vector field network": it takes a
set of backbone coordinates and predicts the flow-matching velocity at every
residue.  The change this module implements is that the network now also
accepts the target conditioning - the per-residue context produced by
:class:`~pfe_design.conditioning.PocketConditioner` - so the same network can
be used for de novo generation and for a design that is steered towards a
specific pocket.

Architecture
------------
Each block is a pair-biased self-attention followed by a feed-forward layer:

* pair features are the radial basis encoding of the inter-residue distances
  of the current coordinates,
* the pair features bias the attention logits, so the network can express
  "residues that are close in space interact",
* invalid (unmodelled) residues are masked out of the attention keys.

The output layer is zero-initialised, so an untrained field predicts zero
velocity and a trajectory starts from the input configuration.  That makes the
masking guarantees in :mod:`pfe_design.integration` testable in isolation from
the quality of the learned field.
"""

from __future__ import annotations

import math

import torch
from torch import Tensor, nn

from .conditioning import ConditioningOutput, PocketConditioner
from .config import ModelConfig, TargetConditioningConfig

__all__ = ["BackboneVelocityField", "PairwiseBiasAttention", "rbf_features"]


def rbf_features(
    distances: Tensor,
    n_bins: int,
    max_distance: float,
    dtype: torch.dtype | None = None,
) -> Tensor:
    """Radial basis encoding of a distance tensor.

    Parameters
    ----------
    distances:
        Distances of any shape ``(..., )``.
    n_bins:
        Number of basis functions.
    max_distance:
        Distance covered by the last basis function (Angstrom).
    dtype:
        The dtype of the result.  It is honoured: a float64 distance tensor
        would otherwise promote the basis to float64 and take the layer's
        float32 weights with it.

    Returns
    -------
    Tensor
        ``(..., n_bins)`` encoding.
    """
    dtype = dtype or distances.dtype
    mu = torch.linspace(0.0, max_distance, n_bins, device=distances.device, dtype=dtype)
    sigma = max_distance / max(1, n_bins - 1)
    return torch.exp(-(((distances.to(dtype).unsqueeze(-1) - mu) / sigma) ** 2))


class PairwiseBiasAttention(nn.Module):
    """Multi-head self-attention whose logits are biased by pair features."""

    def __init__(
        self,
        hidden_dim: int,
        n_heads: int,
        n_rbf_bins: int,
        dropout: float = 0.0,
    ) -> None:
        super().__init__()
        if hidden_dim % n_heads != 0:
            raise ValueError(
                f"hidden_dim ({hidden_dim}) must be divisible by n_heads ({n_heads})"
            )
        self.hidden_dim = hidden_dim
        self.n_heads = n_heads
        self.head_dim = hidden_dim // n_heads
        self.norm = nn.LayerNorm(hidden_dim)
        self.qkv = nn.Linear(hidden_dim, 3 * hidden_dim)
        self.out = nn.Linear(hidden_dim, hidden_dim)
        self.pair_bias = nn.Linear(n_rbf_bins, n_heads)
        self.pair_value = nn.Linear(n_rbf_bins, hidden_dim)
        self.dropout = nn.Dropout(dropout)

    def forward(
        self,
        h: Tensor,
        rbf: Tensor,
        key_mask: Tensor | None = None,
        pair_mask: Tensor | None = None,
    ) -> Tensor:
        """Attend over ``h`` (B, N, hidden) with pair features ``rbf`` (B, N, N, bins).

        ``pair_mask`` marks the residue pairs that both exist (B, N, N).  It
        keeps an unmodelled residue out of the pair features entirely, not
        just out of the attention keys: a residue that is not in the structure
        has no coordinates, so it must not shape the velocities of the
        residues that are.
        """
        batch, n_res, _ = h.shape
        qkv = self.qkv(self.norm(h))
        qkv = qkv.reshape(batch, n_res, 3, self.n_heads, self.head_dim)
        q, k, v = qkv[:, :, 0], qkv[:, :, 1], qkv[:, :, 2]
        q = q.permute(0, 2, 1, 3)
        k = k.permute(0, 2, 1, 3)
        v = v.permute(0, 2, 1, 3)

        logits = torch.matmul(q, k.transpose(-1, -2)) / math.sqrt(self.head_dim)
        logits = logits + self.pair_bias(rbf).permute(0, 3, 1, 2)
        if key_mask is not None:
            logits = logits.masked_fill(~key_mask[:, None, None, :], float("-inf"))
        attn = torch.nan_to_num(logits.softmax(dim=-1))
        attn = self.dropout(attn)

        out = torch.matmul(attn, v).permute(0, 2, 1, 3).reshape(batch, n_res, self.hidden_dim)
        out = self.out(out)
        # Pair-conditioned value term: one channel per hidden unit, averaged
        # over the neighbours of each query residue.  The average is taken
        # over the existing neighbours only, so an unmodelled residue does
        # not dilute the term either - and the invalid pairs are masked
        # before the sum, so even the layer's own bias cannot leak through
        # and make the number of unmodelled residues matter.
        pair_value = self.pair_value(rbf)
        if pair_mask is not None:
            pair_value = pair_value * pair_mask.unsqueeze(-1).to(pair_value.dtype)
            counts = pair_mask.sum(dim=2).clamp(min=1).to(pair_value.dtype)
        else:
            counts = torch.full(
                (batch, n_res),
                float(n_res),
                dtype=pair_value.dtype,
                device=pair_value.device,
            )
        out = out + pair_value.sum(dim=2) / counts.unsqueeze(-1)
        return out


class _Block(nn.Module):
    """One pair-biased attention block followed by a feed-forward layer."""

    def __init__(self, hidden_dim: int, n_heads: int, n_rbf_bins: int, dropout: float) -> None:
        super().__init__()
        self.attn = PairwiseBiasAttention(hidden_dim, n_heads, n_rbf_bins, dropout)
        self.attn_norm = nn.LayerNorm(hidden_dim)
        self.ff_norm = nn.LayerNorm(hidden_dim)
        self.ff = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim * 4),
            nn.SiLU(),
            nn.Linear(hidden_dim * 4, hidden_dim),
        )
        self.dropout = nn.Dropout(dropout)

    def forward(
        self, h: Tensor, rbf: Tensor, key_mask: Tensor | None, pair_mask: Tensor | None = None
    ) -> Tensor:
        h = h + self.dropout(self.attn(self.attn_norm(h), rbf, key_mask, pair_mask))
        h = h + self.dropout(self.ff(self.ff_norm(h)))
        return h


class BackboneVelocityField(nn.Module):
    """Predicts the flow-matching velocity of every residue in the complex.

    Parameters
    ----------
    model_cfg:
        Shape of the network.
    cond_cfg:
        Width of the context vector the conditioner produces; must match the
        conditioner used at run time.
    """

    def __init__(
        self,
        model_cfg: ModelConfig | None = None,
        cond_cfg: TargetConditioningConfig | None = None,
    ) -> None:
        super().__init__()
        self.cfg = model_cfg or ModelConfig()
        self.cond_cfg = cond_cfg or TargetConditioningConfig()
        hidden = self.cfg.hidden_dim

        # The conditioner is owned by the field, so the network can be handed
        # the target structure and the mask directly.
        self.conditioner = PocketConditioner(self.cond_cfg, self.cfg)
        self.type_embedding = nn.Embedding(self.cfg.n_residue_types, hidden)
        self.context_projection = nn.Linear(self.cond_cfg.feature_dim, hidden)
        self.input_norm = nn.LayerNorm(hidden)
        self.blocks = nn.ModuleList(
            _Block(hidden, self.cfg.n_heads, self.cfg.n_rbf_bins, self.cfg.dropout)
            for _ in range(self.cfg.n_layers)
        )
        self.output_norm = nn.LayerNorm(hidden)
        self.to_velocity = nn.Linear(hidden, 3)
        # Zero-initialised: an untrained field is the identity flow.
        nn.init.zeros_(self.to_velocity.weight)
        nn.init.zeros_(self.to_velocity.bias)

    def condition(
        self,
        x: Tensor,
        target_xyz: Tensor | None = None,
        target_valid: Tensor | None = None,
        dynamic_mask: Tensor | None = None,
        pocket_center: Tensor | None = None,
    ) -> ConditioningOutput:
        """Build the conditioning from the target structure and the mask ``M``.

        This is the seam through which the network receives the target: the
        target backbone coordinates and the binary mask go in, and the
        per-residue context comes out.
        """
        if target_xyz is None or dynamic_mask is None:
            raise ValueError(
                "a ConditioningOutput is required: the field is target-conditioned "
                "(or pass target_xyz and dynamic_mask so it can build one)"
            )
        if target_valid is None:
            target_valid = torch.ones(
                target_xyz.shape[:2], device=target_xyz.device, dtype=target_xyz.dtype
            )
        return self.conditioner(
            target_xyz=target_xyz,
            target_valid=target_valid,
            dynamic_mask=dynamic_mask,
            complex_xyz=x,
            pocket_center=pocket_center,
        )

    def forward(
        self,
        x: Tensor,
        residue_types: Tensor | None = None,
        cond: ConditioningOutput | None = None,
        *,
        target_xyz: Tensor | None = None,
        target_valid: Tensor | None = None,
        dynamic_mask: Tensor | None = None,
        pocket_center: Tensor | None = None,
    ) -> Tensor:
        """Predict the velocity at the coordinates ``x``.

        The field is target-conditioned: it needs either a prepared
        :class:`~pfe_design.conditioning.ConditioningOutput`, or the target
        backbone coordinates and the binary mask ``M``, from which it builds
        the conditioning itself.

        Parameters
        ----------
        x:
            Current coordinates of the complex, ``(B, N, 3)``.
        residue_types:
            Residue-type indices, ``(B, N)``; 0 means "unknown".
        cond:
            Conditioning produced by
            :class:`~pfe_design.conditioning.PocketConditioner`.
        target_xyz:
            Target backbone coordinates, ``(B, Nt, 3)``.  Used when ``cond`` is
            not given.
        target_valid:
            1 for a modelled target residue, ``(B, Nt)``.  Defaults to all
            modelled.
        dynamic_mask:
            The binary mask ``M`` over the complex, ``(B, N)``; 1 marks a
            designed residue.  Used when ``cond`` is not given.
        pocket_center:
            Optional pocket centre, ``(B, 3)``.

        Returns
        -------
        Tensor
            Velocities, ``(B, N, 3)``, in the same units as ``x``.
        """
        if x.dim() != 3 or x.shape[-1] != 3:
            raise ValueError(f"x must have shape (batch, n_residues, 3); got {tuple(x.shape)}")
        batch, n_res, _ = x.shape

        if residue_types is None:
            residue_types = torch.zeros((batch, n_res), dtype=torch.long, device=x.device)
        if residue_types.shape != (batch, n_res):
            raise ValueError(
                f"residue_types must have shape {(batch, n_res)}; got {tuple(residue_types.shape)}"
            )
        if int(residue_types.max()) >= self.cfg.n_residue_types:
            raise ValueError(
                f"residue type index {int(residue_types.max())} exceeds "
                f"n_residue_types ({self.cfg.n_residue_types})"
            )

        if cond is None:
            cond = self.condition(
                x,
                target_xyz=target_xyz,
                target_valid=target_valid,
                dynamic_mask=dynamic_mask,
                pocket_center=pocket_center,
            )
        if cond.context.shape[:2] != (batch, n_res):
            raise ValueError(
                f"conditioning context has shape {tuple(cond.context.shape)}, "
                f"expected ({batch}, {n_res}, ...)"
            )
        if cond.context.shape[-1] != self.cond_cfg.feature_dim:
            raise ValueError(
                f"conditioning width {cond.context.shape[-1]} does not match "
                f"feature_dim {self.cond_cfg.feature_dim}"
            )

        h = self.type_embedding(residue_types)
        h = h + self.context_projection(cond.context.to(h.dtype))
        h = self.input_norm(h)

        distances = torch.cdist(x, x)
        rbf = rbf_features(distances, self.cfg.n_rbf_bins, self.cfg.rbf_max, dtype=h.dtype)

        key_mask = cond.valid_mask.to(torch.bool)
        # A residue that is not in the structure has no pair features either.
        pair_mask = key_mask.unsqueeze(-1) & key_mask.unsqueeze(-2)
        rbf = rbf * pair_mask.unsqueeze(-1).to(rbf.dtype)
        for block in self.blocks:
            h = block(h, rbf, key_mask, pair_mask)

        return self.to_velocity(self.output_norm(h))
