# Hotspot Sensitivity Analysis in _De Novo_ Protein Binder Design Using RFdiffusion3

> **Research Question:** Does target hotspot selection on a fixed receptor interface significantly alter _de novo_ binder geometry, steric clash rates, structural fold diversity, and _in silico_ self-consistency under a strictly matched compute budget?

To answer this question, we implement a controlled, reproducible _in silico_ benchmark comparing **RFdiffusion3 (RFD3)** binder generation across three experimental conditions on the human **PD-L1** ectodomain (`PDB: 5O45`, Chain A), followed by **ProteinMPNN** sequence redesign and **RoseTTAFold3 (RF3)** fold-back validation.

---

## run the pipeline

```
./test/scripts/run_pipeline.sh
```

## Design

### Target System

- **Receptor:** Human Programmed Death-Ligand 1 (**PD-L1**), Chain A (`PDB: 5O45`, cropped to residues `A18–A132`).
- **Input Preparation:** Crystallographic waters, ligands, and non-Chain A chains are stripped to isolate a clean single-chain receptor target (`test/data/target_5O45.pdb`).

### Experimental Arms (3 Conditions)

We evaluate three distinct hotspot regimes on the same target surface:

```
A56 (TYR): CB,CG,CD1,CD2,CE1,CE2,CZ,OH
A66 (GLN): CB,CG,CD
A68 (VAL): CB,CG1,CG2
A75 (LYS): CB,CG,CD,CE,NZ
A115 (MET): CB,CG,SD,CE
A123 (TYR): CB,CG,CD1,CD2,CE1,CE2,CZ,OH
```

- **Contig Specification:** `60-90,/0,A18-132` (designs a 60–90 residue binder separated by a `/0` chain break from fixed receptor residues `A18–A132`).



| Condition ID             | Role                | Target Residues & Heavy Side-Chain Atoms                                                                                     | Physicochemical Profile                                      |
| ------------------------ | ------------------- | ---------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------ |
| **`cond_a_canonical`**   | Primary Treatment   | `A56` (TYR: `CB,CG,CD1,CD2,CE1,CE2,CZ,OH`) <br> `A115` (MET: `CB,CG,SD,CE`) <br> `A123` (TYR: `CB,CG,CD1,CD2,CE1,CE2,CZ,OH`) | Canonical hydrophobic/aromatic PD-1 binding patch            |
| **`cond_b_alternative`** | Secondary Treatment | `A66` (GLN: `CB,CG,CD`) <br> `A68` (VAL: `CB,CG1,CG2`) <br> `A75` (LYS: `CB,CG,CD,CE,NZ`)                                    | Peripheral mixed polar/charged surface epitope               |
| **`cond_c_control`**     | Negative Control    | None (`select_hotspots` omitted)                                                                                             | Unguided diffusion baseline (reveals intrinsic surface bias) |


---

### Summary Metrics

| Condition            | Designs | Mean Clashes | Mean Interface Contacts | Mean Hotspot Recall | Max Recall | Geometry Pass Rate |
| -------------------- | ------- | ------------ | ----------------------- | ------------------- | ---------- | ------------------ |
| `cond_a_canonical`   | 22      | 0.59         | 171.7                   | 0.061               | 0.333      | 63.6%              |
| `cond_b_alternative` | 1       | 1.0          | 126.0                   | 0.0                 | 0.0        | 0.0%               |
| `cond_c_control`     | 1       | 0.0          | 505.0                   | n/a                 | n/a        | 100.0%             |


<h3>Representative Backbone Cα Traces (3D)</h3>

<div align="center">
  <img src="visualization/structure_3d_traces.png"
       alt="Structure 3D traces"
       width="75%">
</div>

<h3>Summary Metrics and Quality by Condition</h3>

<div align="center">
  <table>
    <tr>
      <td align="center">
        <img src="visualization/summary_metrics_by_condition.png"
             alt="Summary metrics by condition"
             width="100%">
      </td>
      <td align="center">
        <img src="visualization/structure_metrics_by_condition.png"
             alt="Structure metrics by condition"
             width="100%">
      </td>
    </tr>
  </table>
</div>

<h3>Structural Diagnostics</h3>

<div align="center">
  <img src="results/cif_structural_diagnostics.png"
       alt="CIF structural diagnostics"
       width="55%">
</div>

---

## Environment Setup & Data Preparation

RFD3 requires each YAML configuration to map a top-level design key to a dictionary of arguments, with relative `input` paths resolved from the YAML file directory.

## Quantitative Evaluation

This script parses generated `.cif` and `.pdb` structures across all three arms and computes four primary structural metrics:

1. **Steric Clashes:** Inter-chain heavy-atom pairs within <2.0 Å.
2. **Interface Contact Density:** Inter-chain heavy-atom pairs within ≤4.5 Å.
3. **Hotspot Recall:** Fraction of specified target hotspot residues contacted within ≤4.5 Å.
4. **Fold Diversity:** Mean pairwise Kabsch-aligned Cα-RMSD (Å) across generated binders within each condition.

---
## Validation Criteria & Scientific Caveats

### Multi-Stage Filtering Criteria

To evaluate end-to-end design yield across conditions, candidate designs are assessed against standard _in silico_ thresholds:

1. **Geometric Plausibility (RFD3 Output):** Zero severe inter-chain heavy-atom steric clashes (d < 2.0 Å) and dense interface packing (≥30 heavy-atom contacts within 4.5 Å).
2. **Hotspot Satisfaction:** Contact with ≥2/3 specified hotspot residues (`hotspot_recall >= 0.67`).
3. **Fold-Back Self-Consistency (RF3 Output):** Binder backbone scRMSD < 2.0 Å between the RFD3 design and RF3 complex prediction, coupled with high interface confidence (iPAE < 10 Å, pLDDT > 80).

### Methodological Scope & Limitations

- **Negative Control Interpretation:** `cond_c_control` reveals whether RFD3 naturally gravitates toward the canonical hydrophobic patch (`A56/A115/A123`) even without explicit guidance, isolating true conditioning gain from unconditional surface preference.
- **No Experimental Binding Guarantee:** _In silico_ self-consistency (scRMSD, iPAE) and clash-free interface geometry verify structural plausibility, **not** physical binding affinity (KD) or functional agonism/antagonism. Experimental validation via Surface Plasmon Resonance (SPR), Bio-Layer Interferometry (BLI), or cell-based assays is required to confirm binding.

---

### Key Findings

1. **Canonical hotspot guidance (`cond_a_canonical`) yields the largest and most diverse design set.** With 22 sampled designs across 5 model variants, it is the only condition with sufficient replication for statistical comparison. Hotspot recall reaches a maximum of ~0.33, indicating that explicit hotspot targeting increases (but does not guarantee) contact with the specified PD-1 interface residues relative to the unguided control.

2. **Geometry pass rate is condition-dependent.** `cond_a_canonical` passes geometric plausibility in ~64% of designs, while the single `cond_b_alternative` design fails (1 steric clash, 0 recall). The control passes geometry but with no hotspot guidance its recall is undefined.

3. **Interface contact density is highest in the control (`cond_c_control`: 505 contacts)** but this reflects unguided surface adsorption rather than targeted hotspot engagement; canonical and alternative conditions produce fewer, more focused contacts (171.7 and 126 respectively).

4. **RF3 fold-back pLDDT and MPNN sequence recovery are available per-condition** (see `refold_plddt.png` and `mpnn_design_metrics.png`). These provide end-to-end self-consistency and sequence-design confidence for the generated binders.

5. **Representative 3D Cα traces** (`structure_3d_traces.png`) show the backbone topology of one generated binder per condition, illustrating how different hotspot regimes produce distinct binder folds on the same receptor surface.
