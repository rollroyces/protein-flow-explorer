# pfe-design

Target-specific pocket conditioning and geometry evaluation for protein
design — the modelling half of `protein-flow-explorer`.

The static site in the repository root is an interactive structure viewer.
This package is the part that *designs* a protein: it conditions a
flow-matching backbone field on a chosen pocket of a target protein, integrates
the trajectory with the target held rigid, measures the resulting complex, and
prepares it for side-chain packing.

It is a Python package that lives beside the site, not inside it: the site
still has no build step and no runtime dependency.

## Install

```bash
pip install -e "python[dev]"      # from the repository root
```

Only `torch`, `numpy` and `pydantic` are required. `biopython` is optional
and is used as an independent cross-check of the surface-area calculation.

The test-suite runs straight from a checkout without installing anything:

```bash
cd python && python -m pytest
```

## The mask convention

Everything in this package hangs off one convention, so it is worth stating
once and precisely.

The complex is a single tensor of coordinates with shape `(B, N, 3)`. The
target occupies the first `Nt` positions and the designed region follows:

```
complex index:  0 ............ Nt-1 | Nt ............ N-1
                target (fixed)      | design (dynamic)
```

`M` is the **dynamic** mask: `M[i] == 1` marks a residue that is free to move,
`M[i] == 0` marks a fixed target residue. Its complement is the fixed mask.

```python
from pfe_design import build_complex_mask

M = build_complex_mask(n_target=312, n_design=64)   # 1 on the last 64 columns
```

## Conditioning

```python
import torch
from pfe_design import PocketConditioner, TargetConditioningConfig

conditioner = PocketConditioner(TargetConditioningConfig(feature_dim=32))
cond = conditioner(
    target_xyz=target_backbone,      # (B, Nt, 3)
    target_valid=target_valid,      # (B, Nt) 1 = modelled
    dynamic_mask=M,                  # (B, Nt + Nd)
    pocket_center=pocket_centre,     # (B, 3) the site to aim at
)
```

The pocket is the set of target residues lining the binding site. It is
resolved from an explicit centre, else from the design region's proximity to
the target, else from the target's own centroid. The conditioner is pure
geometry — it never looks at the sequence — so it works for de novo binders
with no prior structure.

The vector field can take the target and the mask directly, and builds the
conditioning itself:

```python
from pfe_design import BackboneVelocityField

field = BackboneVelocityField()
velocity = field(
    complex_xyz,
    residue_types,
    target_xyz=target_backbone,
    target_valid=target_valid,
    dynamic_mask=M,
    pocket_center=pocket_centre,
)
```

`pocket_center` is exactly what the static site's pocket detector reports, so
a cavity found in the viewer can be handed straight to the design package.

## Integration with the target held rigid

```python
from pfe_design import FlowIntegrationConfig, integrate

result = integrate(x0, model, cond, FlowIntegrationConfig(n_steps=100, step_size=0.01))
assert result.max_fixed_velocity == 0.0        # what was applied to the target
assert result.max_fixed_displacement == 0.0    # how far the target moved
print(result.max_raw_fixed_velocity)           # what the model wanted to do to it
```

Two mechanisms enforce the objective's constraint:

1. **Zero velocity** - the predicted velocity is multiplied by the dynamic
   mask before it is used, so a fixed residue never receives a velocity
   whatever the network predicts. This is the mechanism that does the work.
2. **Zero update** - after every step the fixed residues are written back to
   their original coordinates, so the fixed coordinates come out bit-identical
   to the start rather than merely equal to within the arithmetic of the
   update.

They are not independent: with a zero velocity the displacement is already
zero, so removing the write-back changes no observable number. It is kept so
that the guarantee does not depend on that identity holding - a different
float format, a non-finite upstream velocity, or a future solver that
updates the state some other way would all break the implication rather
than the guarantee.

`max_fixed_velocity` is zero by construction, because the mask zeroes the
velocity before it is used. The informative number is
`max_raw_fixed_velocity`: the velocity the model *wanted* on the fixed
residues before the mask was applied, which says how hard the constraint is
pulling.

The mask must be binary. A fractional value would be read as a *weight* by
the velocity mask and as a *flag* by the coordinate write-back, and the two
readings disagree, so both the conditioner and the integrator refuse one.

The integrator also refuses a starting state that does not match the
conditioning: a residue count that disagrees, or a target that is somewhere
else, because the write-back would move the target on the first step.

## Measuring the result

```python
from pfe_design import MetricConfig, buried_surface_area, clash_report, interface_rmsd, ligand_rmsd

irmsd = interface_rmsd(mobile, target, interface_mask, receptor_mask=receptor_mask)
ligand = ligand_rmsd(mobile, target, receptor_mask, ligand_mask)   # Angstrom
bsa = buried_surface_area(complex_xyz, radii, group_a, group_b)  # square Angstrom
clashes = clash_report(complex_xyz, MetricConfig(clash_distance=2.0))
```

The superposition helpers return the rotation and translation as well as the
aligned coordinates, so a ligand can be carried along with the receptor fit
instead of being fitted separately — which is what hides receptor motion in
naive docking scores.

`interface_rmsd` fits on "everything that is not the interface" unless a
`receptor_mask` is given. In a complex that also contains a ligand, pass the
receptor explicitly: otherwise the ligand's own atoms join the fit and tilt
it, and a displaced ligand changes the iRMSD of an interface that never moved.

The surface-area implementation is cross-checked against Biopython's
independent `Bio.PDB.SASA.ShrakeRupley`: the totals agree to 0.6 %, and the
RMS per-atom error shrinks as the sampling is refined, which is the signature
of sampling noise rather than a systematic error. That check runs only when
Biopython is installed; the rest of the suite needs nothing beyond the
declared dependencies.

The buried area is the burial of the two partners *against each other*: the
complex term is the surface of `A ∪ B`, not of every atom in the file. Ions,
solvent and a third chain are not part of that burial, and including them
would make the number negative in the limit - which is not a physical
quantity.

## Side-chain packing

```python
from pfe_design import PackingConfig, SideChainPackingHook

hook = SideChainPackingHook(PackingConfig(backend="mock"))
request = hook.prepare(backbone_xyz, residue_types, fixed_mask)
result = hook(request)
print(result.clash_penalty, result.improved)
```

The hook refines the generated backbone by pushing overlapping atoms apart
while the fixed target residues stay bit-identical, and can then hand the
result to an external packer. PyRosetta, OpenMM and ESMFold are integrated as
stubs with a real seam: each declares whether it is installed, and any of them
accepts an injected callable, so a site-specific installation can be dropped
in without touching this package.

## Layout

```
python/
├── pyproject.toml
├── src/pfe_design/
│   ├── config.py        Pydantic schemas: model, conditioning, integration, metrics, packing
│   ├── conditioning.py  PocketConditioner - X_target and M into per-residue context
│   ├── vector_field.py  BackboneVelocityField - the conditioned vector field
│   ├── integration.py   Constrained ODE/SDE integration
│   ├── metrics.py       iRMSD, ligand RMSD, BSA, steric clashes
│   └── packing.py       Refinement and side-chain packing hook
└── tests/
```
