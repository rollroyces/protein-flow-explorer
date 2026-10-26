"""Pydantic configuration schemas for the target-conditioned design pipeline.

Every stage of the pipeline - conditioning, integration, evaluation and packing -
is driven by a frozen, validated schema.  A run is therefore fully described by
its configuration objects, and an invalid combination (for example a head count
that does not divide the hidden dimension) is rejected before any tensor is
allocated.

The schemas are deliberately independent of each other: a
:class:`FlowIntegrationConfig` can be used on its own to integrate a velocity
field that was built elsewhere.
"""

from __future__ import annotations

from enum import Enum
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

__all__ = [
    "Solver",
    "PackingBackend",
    "BaseSchema",
    "ModelConfig",
    "TargetConditioningConfig",
    "FlowIntegrationConfig",
    "MetricConfig",
    "PackingConfig",
]


class Solver(str, Enum):
    """Integrator used to advance the flow-matching trajectory."""

    EULER = "euler"
    HEUN = "heun"


class PackingBackend(str, Enum):
    """External side-chain packer that can be driven by the packing hook."""

    NONE = "none"
    MOCK = "mock"
    ROSETTA = "rosetta"
    OPENMM = "openmm"
    ESMFOLD = "esmfold"


class BaseSchema(BaseModel):
    """Base class for every configuration schema.

    ``frozen`` makes instances hashable and immutable, so a configuration can be
    attached to a run's provenance record.  ``extra="forbid"`` turns a typo in a
    keyword into an error rather than a silently ignored setting.
    """

    model_config = ConfigDict(
        frozen=True,
        extra="forbid",
        validate_assignment=True,
        arbitrary_types_allowed=False,
    )


class ModelConfig(BaseSchema):
    """Shape of the backbone vector-field network.

    The network is intentionally small by default: it is a reference
    implementation of the conditioning contract, sized so that the test-suite
    runs on CPU in seconds while still exercising every code path.
    """

    hidden_dim: int = Field(default=128, gt=0, le=4096)
    n_layers: int = Field(default=4, ge=1, le=64)
    n_heads: int = Field(default=4, ge=1, le=64)
    n_rbf_bins: int = Field(default=16, ge=2, le=128)
    rbf_max: float = Field(default=20.0, gt=0.0, le=200.0)
    dropout: float = Field(default=0.0, ge=0.0, lt=1.0)
    n_residue_types: int = Field(default=22, ge=2, le=64)
    """Number of residue-type embedding slots (index 0 is reserved for unknown)."""

    @model_validator(mode="after")
    def _validate_heads(self) -> ModelConfig:
        if self.hidden_dim % self.n_heads != 0:
            raise ValueError(
                f"hidden_dim ({self.hidden_dim}) must be divisible by "
                f"n_heads ({self.n_heads})"
            )
        return self

    @property
    def head_dim(self) -> int:
        """Dimension of a single attention head."""
        return self.hidden_dim // self.n_heads


class TargetConditioningConfig(BaseSchema):
    """How the target protein and its pocket are turned into context features.

    ``pocket_cutoff`` is the distance (Angstrom) within which a target residue
    is considered to line the pocket.  ``feature_dim`` is the width of the
    per-residue context vector handed to the vector field.
    """

    pocket_cutoff: float = Field(default=12.0, gt=0.0, le=100.0)
    """Distance within which a target residue lines the pocket (Angstrom)."""

    max_pocket_residues: int | None = Field(default=None, ge=1, le=4096)
    """Cap on how many target residues may line the pocket (closest first)."""

    feature_dim: int = Field(default=32, ge=1, le=1024)
    dropout: float = Field(default=0.0, ge=0.0, lt=1.0)
    normalize_features: bool = True
    """Scale geometric features by the pocket radius instead of raw Angstroms."""

    detach_target: bool = True
    """Detach the target coordinates, so no gradient flows back into the target."""

    use_pocket_membership: bool = True
    """Include a per-residue flag marking the residues that line the pocket."""


class FlowIntegrationConfig(BaseSchema):
    """Settings for the constrained ODE / SDE integration of the vector field.

    The hard constraints (zero velocity and zero displacement on fixed target
    residues) are always applied; ``enforce_hard_constraints`` only controls
    whether the target coordinates are additionally re-imposed exactly after
    each step, which protects against floating-point drift.
    """

    n_steps: int = Field(default=100, ge=1, le=1_000_000)
    step_size: float = Field(default=0.01, gt=0.0, le=10.0)
    solver: Solver = Solver.EULER
    stochastic: bool = False
    noise_scale: float = Field(default=0.0, ge=0.0, le=10.0)
    """Standard deviation of the diffusion term when ``stochastic`` is enabled."""

    seed: int | None = Field(default=None, ge=0)
    device: str = Field(default="cpu")
    dtype: Literal["float32", "float64"] = "float32"
    enforce_hard_constraints: bool = True
    max_displacement: float | None = Field(default=None, gt=0.0)
    """Optional cap on the per-step displacement of any residue (Angstrom)."""

    @model_validator(mode="after")
    def _validate_noise(self) -> FlowIntegrationConfig:
        if self.stochastic and self.noise_scale <= 0.0:
            raise ValueError("stochastic integration requires noise_scale > 0")
        return self


class MetricConfig(BaseSchema):
    """Thresholds and sampling parameters for the geometry evaluation metrics."""

    clash_distance: float = Field(default=2.0, gt=0.0, le=10.0)
    """Any inter-atomic distance below this is counted as a steric clash."""

    clash_mode: Literal["distance", "vdw"] = "distance"
    """``distance``: fixed threshold.  ``vdw``: sum of van der Waals radii."""

    vdw_tolerance: float = Field(default=0.4, ge=0.0, le=2.0)
    """Overlap subtracted from the summed radii in ``vdw`` mode."""

    exclude_intra_residue_pairs: bool = True
    """Ignore pairs of atoms that belong to the same residue.

    A covalent bond is shorter than any clash threshold, so counting it would
    report a clash for every residue in every structure.  Only has an effect
    when residue indices are supplied.
    """

    bsa_probe_radius: float = Field(default=1.4, gt=0.0, le=5.0)
    """Radius of the solvent probe used for the Shrake-Rupley surface (Angstrom)."""

    bsa_n_points: int = Field(default=100, ge=8, le=10_000)
    """Sphere points per atom for the Shrake-Rupley surface sampling."""


class PackingConfig(BaseSchema):
    """Settings for the post-processing refinement and packing hook."""

    backend: PackingBackend = PackingBackend.MOCK
    max_iterations: int = Field(default=100, ge=1, le=10_000)
    step_size: float = Field(default=0.05, gt=0.0, le=10.0)
    clash_distance: float = Field(default=2.0, gt=0.0, le=10.0)
    """Distance below which the refinement pushes atoms apart."""

    min_improvement: float = Field(default=1e-6, ge=0.0)
    """A refinement is only accepted if the clash penalty falls by at least this."""

    keep_fixed_residues: bool = True
    """Fixed target residues stay bit-identical through the refinement."""

    reject_worse: bool = True
    """Discard a backend result that is worse than the input (by clash penalty)."""
