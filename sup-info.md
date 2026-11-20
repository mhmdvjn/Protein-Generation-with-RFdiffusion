
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

---

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

---

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

---

### Project 4 — Targeted Protein Design

**Focus:** One target family (e.g. a kinase, a binder target, an enzyme family).

**Deliver:** Ranked, validated candidates.

**Supplementary info:**

- **Target family:** Choose one biologically defined family and define a consistent design task (e.g. scaffold a known active-site motif, design a binder to a known epitope, or generate a symmetric oligomer).
- **Ranking:** Rank candidates using orthogonal metrics only — e.g. motif recapitulation RMSD, RF3/Chai-1 predicted confidence (pLDDT/IPAE), clash score, design diversity. Do not invent binding affinities, IC50 values, or efficacy scores.
- **Validation:** Report which candidates passed structural validation (geometry + recapitulation) and provide structures (PDB/CIF) for the top-ranked candidates.
- **Rule:** No invented affinities or efficacy scores. Any functional claim requires wet-lab measurement and is outside the scope of this computational project.
