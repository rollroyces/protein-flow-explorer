# Protein Flow Explorer

An interactive tool for visualising, measuring and explaining 3D protein structures —
built for teaching, extended for research.

> **Live:** https://rollroyces.github.io/protein-flow-explorer/
>
> **Audit:** [AUDIT.md](AUDIT.md) — structured code and architecture review
> (performance, architecture, robustness, DevOps) with a prioritised top 3.

## How to present a structure bioinformatically

The guiding principle: **a deposited structure is a linked record, not a picture.**

Most viewers present a model as an object to be admired. Bioinformatics presents it as
data to be interrogated, where every visual choice encodes a stated quantity and every
claim resolves to a citable identifier. The **Presentation** tab in the app implements
this directly; the method is worth stating independently of the code:

### 1. Present the evidence hierarchy, in order

A structure is nested data. Show each level with its own provenance, never flattened
into one image:

```
Entry           PDB ID, method, resolution      ← is the model right?
 └─ Entities     polymer entities, chains      ← what is in it?
     └─ Ligands  non-polymer entities         ← what binds it?
         └─ Residues  modelled, SS annotation  ← what is resolved?
             └─ Atoms  coordinates, B-factor  ← how confident are they?
```

Residue count, chain list and ligand identities each belong at their own level. A viewer
that shows "574 residues" without saying which chains are modelled, or which ligands
were co-crystallised, is withholding the provenance that makes the number meaningful.

### 2. Make colour quantitative and declared

Every colour must encode a named quantity, and the user must be told which:

| Colour scheme | Encodes | Use it to |
|---|---|---|
| By chain | polymer entity | locate a subunit across a complex |
| Secondary structure | helix / strand / coil | read fold architecture |
| Spectrum (N→C) | residue index | follow a sequence walk, spot domain boundaries |
| B-factor | atomic displacement (Å²) | read flexibility in the crystal |
| pLDDT | prediction confidence (0–100) | judge an AlphaFold model's trust |
| Conservation | MSA column identity | find functional sites |
| Hydrophobicity | Kyte–Doolittle | find a ligand-facing face |

Decorative colour is a lie: it implies a signal that is not there. The app renders this
as a clickable **Annotation key** — selecting a row applies that scheme to the viewer.

### 3. State two independent kinds of uncertainty

These are routinely conflated, and conflating them is how bad conclusions get made:

- **Is the model right?** — method and resolution bound how much you can trust
  coordinate positions. ≤ 1.8 Å: side-chain rotamers are usually reliable.
  > 3.0 Å: loops and distal side chains may be unmodelled entirely.
- **Is this region ordered?** — B-factor or pLDDT bound how much you can trust the
  *conformation* of a specific stretch. A well-resolved entry can still contain long
  disordered loops with well-measured high B-factors.

A conclusion drawn from a high-B-factor loop at 2.5 Å resolution fails on both counts.

### 4. Keep identifiers resolvable

Every claim must resolve to a stable ID: PDB entry, UniProt accession (via SIFTS),
chain and residue numbers, ClinVar variant. A finding that cannot be cited cannot be
reproduced. The **Cross-references** panel lists them for the loaded record, and the app
generates a one-line citation.

### 5. Quantify alongside the 3D view

The 3D pane is one projection. The payload is the numbers: RMSD per residue, pocket
volume and enclosure, contact lists with distances, physicochemical descriptors,
conservation scores, mutation impact. Both belong on screen together — the picture
without the numbers is not evidence.

## Features

### For students
- **Representation** — cartoon / surface / sticks / spacefill / ribbon / lines, with
  six colour schemes
- **Guided tours** — scripted narrated walkthroughs for the curated proteins
- **Trajectory playback** — hemoglobin T→R, adenylate kinase open→closed,
  calmodulin extended→compact, HIV-1 protease apo vs inhibitor-bound

### For researchers
- **Measurement** — click atoms for distances (Å), angles (°), dihedrals (τ)
- **Sequence panel** — click any residue to locate and highlight it in 3D; SS overlay
- **UniProt / AlphaFold lookup** — type `TP53` or `P04637`
- **Structure comparison** — superimpose on Cα, colour by per-residue RMSD
- **Conservation (MSA)** — ~10 UniProt homologs aligned in-browser
- **pLDDT colouring** — AlphaFold confidence from the B-factor column
- **Variant map** — ClinVar variants placed on the structure: canonical UniProt
  positions mapped through the construct's own sequence, wild-type checked
  against the residue, coloured by clinical significance (pathogenic, drug
  response, VUS…), with the population allele frequency from gnomAD
- **Tissue expression** — median expression across 54 human tissues from GTEx,
  grouped into compartments, with a specificity verdict: is the protein even
  made in the tissue where the disease happens?
- **Pocket detection** — Fpocket-style concavity scan, validated against ground truth
- **Ligand contacts** — H-bonds (3.5 Å), salt bridges (4.0 Å), π (5.5 Å), hydrophobic
- **Ligand properties** — MW, HBD/HBA, TPSA, cLogP, Lipinski Ro5, Veber
- **Mutation impact** — packing density and pocket membership for a point mutation
- **Export** — visible chains as PDB (PyMOL / ChimeraX), measurements as CSV
- **Share link** — the current view encoded in the URL

## Library

13 curated entries with educational context: hemoglobin, lysozyme, HIV-1 protease,
adenylate kinase, myoglobin, calmodulin, insulin, T4 lysozyme, p53, EGFR + erlotinib,
KRAS G12C + sotorasib, type II collagen, aggrecan.

Load anything else by PDB ID, gene symbol, UniProt accession, or by dropping a
`.pdb` / `.cif` file.

## Development

```bash
python3 -m http.server 8765   # then open http://localhost:8765
npm test                      # 1,316 assertions across 10 suites
npm run check                 # syntax-check every JS file
npx eslint@9 .                # lint (the config is in the repo, no dependency)

cd python && python -m pytest # the design package's own suite
```

No build step, no dependencies. Vanilla HTML/CSS/JS + 3Dmol.js.

The scientific core lives in `src/` as eight dependency-free modules —
`seq.js` (sequence parsing, alignment, MSA, conservation), `geometry.js`
(distance/angle/dihedral/RMSD), `pockets.js` (spatial hash, density gates,
ray casting), `colors.js` (the colour ramps behind every scheme in the
colour dropdown), `variants.js` (ClinVar record parsing, significance
classification, the UniProt→PDB residue map), `expression.js` (GTEx tissue
expression and specificity), `validate.js` (input and payload validation) and
`annotations.js` (the annotation key, cross-checked against the palettes the
app actually renders). They load as plain scripts and publish themselves on
the `PFE` namespace, so they can be `require()`d and unit-tested with no DOM
stub at all. See [AUDIT.md](AUDIT.md) for the architecture review that
motivated the split.

## Design package

The site reads structures; it cannot design anything. `python/` is the other
half — a PyTorch package that conditions a flow-matching backbone field on a
chosen pocket of a target protein, integrates the trajectory with the target
held rigid, measures the complex, and prepares it for side-chain packing.

It is a separate package with its own `pyproject.toml`, its own pytest suite
and its own CI job. The site above is unchanged by it: still no build step,
still no runtime dependency.

```python
from pfe_design import (
    build_complex_mask,
    PocketConditioner, TargetConditioningConfig,
    BackboneVelocityField, ModelConfig,
    FlowIntegrationConfig, integrate,
    interface_rmsd, buried_surface_area, clash_report,
    SideChainPackingHook,
)

M = build_complex_mask(n_target=312, n_design=64)   # 1 = dynamic (designed)
cond = PocketConditioner(TargetConditioningConfig())(
    target_xyz=target_backbone, target_valid=target_valid,
    dynamic_mask=M, pocket_center=site,
)
result = integrate(x0, model, cond, FlowIntegrationConfig(n_steps=100))
assert result.max_fixed_velocity == 0.0        # the target never moves
```

The design region is steered by a binary mask `M` over the complex, and the
fixed target residues are held by two mechanisms: their velocity is zeroed, and
their coordinates are written back after every integration step, so the target
cannot drift away from where it started. The pocket centre is exactly what the
site's pocket detector reports, so a cavity found in the viewer can be handed
straight to the design package. The vector field is an untrained reference
implementation — the point of the package is the conditioning contract and the
measurement, not a trained model. See [python/README.md](python/README.md) for
the details and [AUDIT.md § 12](AUDIT.md) for the reasoning behind each metric
definition, including what an adversarial review found and what it cost.

## Notes on method

- **Pocket detection** follows Fpocket's concavity principle rather than raw density:
  a pocket must sit in a genuine void, have walls, be globally enclosed, and clear a
  32-ray enclosure test that separates concave pockets from convex grooves. Validated
  on EGFR/erlotinib (1M17), where the top-ranked pocket is fully occupied by the
  crystallographic ligand. Reported volumes are grid estimates — use them to rank,
  not as absolute numbers.
- **Packing density** uses the conventional CA-to-CA definition within 8 Å. Thresholds
  are calibrated against oxymyoglobin (1MBO): buried ≥ 10, moderate ≥ 6.
- **R-factors** are read from the `refine` / `refine_hist` blocks under capital-R
  mmCIF keys; entries that genuinely lack them show an em dash.
- **Hydrogen bonds** are geometric N/O–N/O contacts under 3.5 Å, not computed with
  explicit hydrogens. Treat as a shortlist, not a definitive list.

## Licence

Educational tool, not for clinical use. Structure data from
[RCSB PDB](https://www.rcsb.org) and [AlphaFold DB](https://alphafold.ebi.ac.uk);
identifiers from [UniProt](https://www.uniprot.org) and
[ClinVar](https://www.ncbi.nlm.nih.gov/clinvar). Rendering by
[3Dmol.js](https://3dmol.org).
