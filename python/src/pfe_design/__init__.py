"""Target-specific pocket conditioning and geometry evaluation for protein design.

The package provides the conditioning half of a flow-matching design pipeline:

* :class:`~pfe_design.conditioning.PocketConditioner` turns a target structure
  and a binary complex mask into per-residue context.
* :class:`~pfe_design.vector_field.BackboneVelocityField` predicts the
  flow-matching velocity, conditioned on that context.
* :func:`~pfe_design.integration.integrate` advances the trajectory with the
  fixed target residues held exactly in place.
* :mod:`pfe_design.metrics` measures the result: iRMSD, ligand RMSD, buried
  surface area and steric clashes.
* :class:`~pfe_design.packing.SideChainPackingHook` refines the generated
  backbone and hands it to a side-chain packer.
"""

from __future__ import annotations

from .conditioning import (
    ConditioningOutput,
    PocketConditioner,
    build_complex_mask,
    dynamic_mask_from_layout,
)
from .config import (
    BaseSchema,
    FlowIntegrationConfig,
    MetricConfig,
    ModelConfig,
    PackingBackend,
    PackingConfig,
    Solver,
    TargetConditioningConfig,
)
from .integration import (
    IntegrationResult,
    VelocityModel,
    apply_mask_constraints,
    integrate,
    sample_binder,
)
from .metrics import (
    VDW_RADII,
    BSAResult,
    ClashReport,
    ComplexEvaluation,
    apply_transform,
    buried_surface_area,
    clash_penalty,
    clash_penalty_fn,
    clash_report,
    evaluate_complex,
    interface_rmsd,
    ligand_rmsd,
    sasa,
    superpose,
    vdw_radii,
)
from .packing import (
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
from .sequences import (
    AA1_ORDER,
    canonical_sequence,
    residue_types_from_sequence,
    sequence_from_residue_types,
)
from .vector_field import BackboneVelocityField, PairwiseBiasAttention, rbf_features

__version__ = "0.1.0"

__all__ = [
    "__version__",
    # configuration
    "BaseSchema",
    "FlowIntegrationConfig",
    "MetricConfig",
    "ModelConfig",
    "PackingBackend",
    "PackingConfig",
    "Solver",
    "TargetConditioningConfig",
    # conditioning
    "ConditioningOutput",
    "PocketConditioner",
    "build_complex_mask",
    "dynamic_mask_from_layout",
    # vector field and integration
    "BackboneVelocityField",
    "IntegrationResult",
    "PairwiseBiasAttention",
    "VelocityModel",
    "apply_mask_constraints",
    "integrate",
    "rbf_features",
    "sample_binder",
    # metrics
    "BSAResult",
    "ClashReport",
    "ComplexEvaluation",
    "VDW_RADII",
    "apply_transform",
    "buried_surface_area",
    "clash_penalty",
    "clash_penalty_fn",
    "clash_report",
    "evaluate_complex",
    "interface_rmsd",
    "ligand_rmsd",
    "sasa",
    "superpose",
    "vdw_radii",
    # sequences
    "AA1_ORDER",
    "canonical_sequence",
    "residue_types_from_sequence",
    "sequence_from_residue_types",
    # packing
    "BackendUnavailableError",
    "ESMFoldBackend",
    "MockPackingBackend",
    "OpenMMBackend",
    "PackingBackendBase",
    "PackingRequest",
    "PackingResult",
    "RosettaBackend",
    "SideChainPackingHook",
    "build_backend",
    "relax_clashes",
]
