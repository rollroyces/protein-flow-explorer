# Protein Flow Explorer

An interactive web tool for visualizing and explaining 3D protein structures. Built for
teaching, presenting, and exploring molecular biology.

> **Live demo:** [https://rollroyces.github.io/protein-flow-explorer](https://rollroyces.github.io/protein-flow-explorer)

## What it does

Protein Flow Explorer makes it easy to **understand and explain** 3D biological structures.
It bundles three complementary modes in a single web app:

1. **Interactive Controls** — switch between cartoon / surface / sticks / spacefill
   representations and colour schemes (by chain, secondary structure, hydrophobicity,
   B-factor, …) to inspect any structure from any angle.

2. **Guided Tour** — a scripted, step-by-step walkthrough for each curated protein.
   Narration cards explain *why* each feature matters, while the camera zooms in on
   active sites, heme groups, or domain interfaces.

3. **Trajectory Playback** — slide or auto-play between two PDB states of the same protein
   (hemoglobin T→R, kinase open→closed, calmodulin extended→compact, HIV-1 protease
   apo vs inhibitor-bound). Drag the slider to scrub.

## Curated library

8 teaching-grade proteins are preloaded with educational context:

| Protein                  | PDB IDs         | Tour                              | Trajectory                      |
| ----------------------- | --------------- | --------------------------------- | ------------------------------- |
| Hemoglobin              | 2DN2, 1HHO      | Oxygen-binding journey            | T → R allosteric transition     |
| Hen Egg-White Lysozyme  | 1LYZ            | Catalytic mechanism               | —                               |
| HIV-1 Protease          | 1HVR, 2PC1      | Drug-target walkthrough           | Apo vs inhibitor-bound          |
| Adenylate Kinase        | 4AKE, 1AKE      | Induced-fit domain motion         | Open → closed                   |
| Sperm Whale Myoglobin   | 1MBO            | First solved structure            | —                               |
| Calmodulin              | 1CLL, 1CDL      | Calcium-triggered switch          | Extended → compact              |
| Insulin Hexamer         | 4INS            | From hexamer to hormone           | —                               |
| T4 Lysozyme             | 1L63            | Folding stability sandbox         | —                               |

All structures are loaded live from the open-access
[RCSB Protein Data Bank](https://www.rcsb.org) — no server-side storage required.

## Tech

- **Rendering:** [3Dmol.js](https://3dmol.org) (WebGL)
- **Data:** RCSB PDB
- **Frontend:** Vanilla HTML + CSS + JavaScript, no build step
- **Hosting:** GitHub Pages

## Local development

```bash
# Clone or download this folder
cd protein-viz

# Serve any way you like
python3 -m http.server 8765
# then open http://localhost:8765
```

That's it — no `npm install`, no bundler.

## File structure

```
protein-viz/
├── index.html        # Single-page app entry
├── styles.css        # All styles
├── app.js            # Application logic (viewer, tours, trajectories)
├── data/
│   └── proteins.js   # Protein library, tours, trajectories
└── README.md
```

## License

Open-access teaching tool. PDB structures are public domain; 3Dmol.js is BSD-licensed.
Please cite RCSB PDB and 3Dmol.js when reusing this work.