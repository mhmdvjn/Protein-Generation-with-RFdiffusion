# Protein Generation with RFdiffusion


./test/scripts/run_pipeline.sh

Student project repository for **RFdiffusion**-based de novo protein generation and design.
University course — Bioinformatics / Computational Biology, Semester 1.

RFdiffusion is a diffusion (flow-matching) generative model that designs novel protein backbones conditioned on structural and functional constraints, introduced by the Baker lab (UW Medicine Institute for Protein Design, RosettaCommons). This repo contains four guided student projects plus local copies of the upstream RFdiffusion implementations (RFdiffusion2 and RFdiffusion3 / RC-Foundry) for reference and inference.

## Repository structure

| Path | Contents |
|------|----------|
| `materials/` | Supplementary reading, paper PDFs, lecture notes, and reference docs (to be populated) |
| `project-1/` | Hotspot sensitivity: does hotspot choice change design quality? |
| `project-2/` | Atomic constraints: ablation of atom-level motif constraints (geometry, clashes, diversity) |
| `project-3/` | End-to-end design and validation pipeline: RFD3 → MPNN → RF3 |
| `project-4/` | Targeted protein design for one target family; ranked, validated candidates |
| `RFdiffusion2/` | Upstream RFdiffusion2 code (RosettaCommons). Apptainer/Singularity image-based inference; atom-level enzyme active-site scaffolding. See `RFdiffusion2/README.md` |
| `RFdiffusion3/` | Upstream RFdiffusion3 / RC-Foundry code. RFD3 (design), RF3 (folding), ProteinMPNN/LigandMPNN (inverse folding). See `RFdiffusion3/README.md` |

## Background: RFdiffusion

RFdiffusion generates protein backbones by iteratively denoising a random structure under user-specified constraints (motif scaffolding, binder design, symmetry, hotspot constraints, etc.). The diffusion process is followed by sequence design with an inverse-folding model (e.g. ProteinMPNN / LigandMPNN) and structure validation with a folding model (e.g. RF3 / Chai-1 / AlphaFold3).

Key papers:

- **RFdiffusion (original):** Watson et al., *De novo design of protein structure and function with RFdiffusion*, Nature, 2023.
- **RFdiffusion2:** Ahern et al., *Atom level enzyme active site scaffolding using RFdiffusion2*, bioRxiv, 2025.
- **RFdiffusion3:** Butcher et al., *De novo Design of All-atom Biomolecular Interactions with RFdiffusion3*, bioRxiv / Science, 2025–2026.
- **ProteinMPNN:** Dauparas et al., *Robust deep learning-based protein sequence design using ProteinMPNN*, Science, 2022.
- **LigandMPNN:** Dauparas et al., *Atomic context-conditioned protein sequence design using LigandMPNN*, Nature Methods, 2025.
- **RosettaFold3:** RF3 structure prediction for validation of designs.

See each project directory for its own `README.md` with the research question, expected output, and supplementary info.

## Projects

### Project 1 — Hotspot Sensitivity

**Question:** Does hotspot choice change design quality?

**Pipeline:** Generate matched designs using different hotspot residue selections while holding the sampling budget equal; compare designs on structure and sequence metrics.

**Output:** Matched designs + metric comparison (e.g. backbone RMSD to motif, pLDDT/IPAE, sequence diversity, ligand-contact recovery).

**Control:** Equal sampling budget across conditions (same number of diffusion trajectories and MPNN sequences per hotspot).

**Supplementary info:**
- **Hotspots:** Residues (or atoms) specified as functional constraints that the generator must scaffold or preserve. In RFdiffusion these appear as "contigs"/motif constraints; in RFdiffusion2/RFD3 they can be atom-level tip atoms or side-chain rotamer constraints.
- **Why control the budget?** Design quality scales with compute (more trajectories = better samples). Without an equal-budget control, differences between hotspots could be caused by sampling effort rather than the hotspot choice itself.
- **Metrics to report:** motif recapitulation RMSD, all-atom clash score, per-residue pLDDT, interface RMSD (for binders), sequence identity/diversity among the top designs, and optional experimental validation (e.g. expression/solubility) — not required for this project.
- **Deliverable:** Code used to generate and compare designs, the methods description, and figures (trajectory overlay, metric boxplots, PyMOL renders).

### Project 2 — Atomic Constraints

**Question:** Does adding atom-level constraints improve designs?

**Compare:** Designs generated with vs. without atomic constraints on the same target.

**Measure:** Geometry (bond/angle/torsion plausibility), steric clashes, sequence/structure diversity.

**Supplementary info:**
- **Atom-level constraints:** In RFdiffusion2 and RFD3, individual atoms (e.g. catalytic triad atoms, ligand heavy atoms, metal-coordinating atoms) can be fixed as "atomized" regions while the rest of the backbone diffuses. This constrains active-site geometry directly instead of relying on residue-level hotspots.
- **Ablation design:** Keep target, sampling budget, noise schedule, and MPNN settings identical; toggle only the atomic constraint set.
- **Geometry metrics:** Ramachandran outliers, rotamer outliers, bond-length/angle deviations (e.g. via MolProbity or OpenFold/AF3 geometry checks).
- **Clash metrics:** Number and severity of steric overlaps (< 2.0–2.4 Å non-bonded contacts), especially at motif–backbone and motif–ligand interfaces.
- **Diversity metrics:** Pairwise C-alpha RMSD / TM-score among generated designs; sequence entropy at design positions.
- **Status:** Ablation study only — not a clinical claim. Do not report binding affinities or efficacy unless experimentally measured.

### Project 3 — Design and Validate

**Pipeline:** RFD3 → MPNN → RF3

**Result:** Reproducible validation report.

**Deliver:** Code + methods + figures.

**Supplementary info:**
- **RFD3 (RFdiffusion3):** All-atom generative model producing constrained biomolecular structures (proteins, and with RFdiffusion3NA, nucleic acids). Outputs a backbone + side-chain trajectory under complex constraints.
- **MPNN (ProteinMPNN / LigandMPNN):** Inverse-folding step that designs a sequence compatible with the generated backbone, optionally conditioned on ligands or motif rotamers (LigandMPNN).
- **RF3 (RosettaFold3):** Structure prediction model used to re-fold the MPNN-designed sequence and validate that the predicted structure recapitulates the designed backbone and constraints (compare predicted vs. designed C-alpha RMSD, pLDDT, IPAE).
- **Validation report should include:** design inputs (constraints/contigs), generation parameters (seed, noise scale, number of samples), MPNN sequences, RF3 predicted structures, recapitulation metrics, clash geometry, and PyMOL figures of the diffusion trajectory and final models.
- **Reproducibility:** Pin random seeds, software versions (conda/env or Apptainer image), model checkpoint versions, and record all CLI arguments in the report.

### Project 4 — Targeted Protein Design

**Focus:** One target family (e.g. a kinase, a binder target, an enzyme family).

**Deliver:** Ranked, validated candidates.

**Supplementary info:**
- **Target family:** Choose one biologically defined family and define a consistent design task (e.g. scaffold a known active-site motif, design a binder to a known epitope, or generate a symmetric oligomer).
- **Ranking:** Rank candidates using orthogonal metrics only — e.g. motif recapitulation RMSD, RF3/Chai-1 predicted confidence (pLDDT/IPAE), clash score, design diversity. Do not invent binding affinities, IC50 values, or efficacy scores.
- **Validation:** Report which candidates passed structural validation (geometry + recapitulation) and provide structures (PDB/CIF) for the top-ranked candidates.
- **Rule:** No invented affinities or efficacy scores. Any functional claim requires wet-lab measurement and is outside the scope of this computational project.

## Getting started

The upstream implementations carry their own installation instructions:

- **RFdiffusion2** (`RFdiffusion2/`): Apptainer/Singularity image (`bakerlab_rf_diffusion_aa.sif`), model weights downloaded via `python setup.py`, inference via `rf_diffusion/benchmark/pipeline.py`. See `RFdiffusion2/README.md` and https://rosettacommons.github.io/RFdiffusion2/.
- **RFdiffusion3 / RC-Foundry** (`RFdiffusion3/`): `pip install "rc-foundry[all]"`, weights via `foundry install base-models`. Supports Linux (CUDA), Intel XPU, and macOS Apple Silicon (MPS, community fork). See `RFdiffusion3/README.md`.

General requirements: NVIDIA GPU strongly recommended for inference; CPU-only runs are extremely slow. Python 3.10+, conda/mamba or uv, and Apptainer/Singularity for RFdiffusion2.

## Materials

The `materials/` folder is intended for supplementary course material: paper PDFs, lecture slides, constraint/contig syntax references, metric definitions, and PyMOL visualization scripts. It is currently empty — populate it as the course progresses.

## Citation

If you use the upstream code in this repo, cite the corresponding work:

```bibtex
@article{watson2023de,
  title={De novo design of protein structure and function with RFdiffusion},
  author={Watson, Joseph L and Juergens, David and Bennett, Nathaniel R and Trippe, Brian L and Yim, Jason and Eisenach, Helen E and Ahern, Woody and Hersey, Andrew J and Woetzi, Emma and Han, Jue and others},
  journal={Nature},
  volume={620},
  number={7976},
  pages={1089--1100},
  year={2023},
  publisher={Nature Publishing Group UK London}
}

@article{ahern2025atom,
  title={Atom level enzyme active site scaffolding using RFdiffusion2},
  author={Ahern, Woody and Yim, Jason and Tischer, Doug and Salike, Saman and Woodbury, Seth M and Kim, Donghyo and Kalvet, Indrek and Kipnis, Yakov and Coventry, Brian and Altae-Tran, Han Raut and others},
  journal={bioRxiv},
  pages={2025--04},
  year={2025},
  publisher={Cold Spring Harbor Laboratory}
}

@article{butcher2025_rfdiffusion3,
  author = {Butcher, Jasper and Krishna, Rohith and Mitra, Raktim and Brent, Rafael Isaac and Li, Yanjing and Corley, Nathaniel and Kim, Paul T and Funk, Jonathan and Mathis, Simon Valentin and Salike, Saman and Muraishi, Aiko and Eisenach, Helen and Thompson, Tuscan Rock and Chen, Jie and Politanska, Yuliya and Sehgal, Enisha and Coventry, Brian and Zhang, Odin and Qiang, Bo and Didi, Kieran and Kazman, Maxwell and DiMaio, Frank and Baker, David},
  title = {De novo Design of All-atom Biomolecular Interactions with RFdiffusion3},
  elocation-id = {2025.09.18.676967},
  year = {2025},
  doi = {10.1101/2025.09.18.676967},
  publisher = {Cold Spring Harbor Laboratory},
  journal = {bioRxiv}
}

@article{dauparas2022robust,
  title={Robust deep learning--based protein sequence design using ProteinMPNN},
  author={Dauparas, Justas and Anishchenko, Ivan and Bennett, Nathaniel and Bai, Hua and Ragotte, Robert J and Milles, Lukas F and Wicky, Basile IM and Courbet, Alexis and de Haas, Rob J and Bethel, Neville and others},
  journal={Science},
  volume={378},
  number={6615},
  pages={49--56},
  year={2022},
  publisher={American Association for the Advancement of Science}
}
```

## Disclaimer

This repository is for educational coursework. Generated protein designs are computational predictions only. Do not use them for therapeutic, diagnostic, or clinical purposes without experimental validation and appropriate regulatory review. Projects 2–4 explicitly prohibit reporting invented affinities or efficacy scores.

## License

Upstream subdirectories (`RFdiffusion2/`, `RFdiffusion3/`) retain their original licenses (see their `LICENSE.md` files). The student project content (`project-1/` through `project-4/`, `materials/`, and this README) is provided for course use.
