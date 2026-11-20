# Protein Generation with RFdiffusion


Student project repository for **RFdiffusion**-based de novo protein generation and design.
University course — Bioinformatics / Computational Biology, Semester 1.

RFdiffusion is a diffusion (flow-matching) generative model that designs novel protein backbones conditioned on structural and functional constraints, introduced by the Baker lab (UW Medicine Institute for Protein Design, RosettaCommons). This repo contains four guided student projects plus local copies of the upstream RFdiffusion implementations (RFdiffusion2 and RFdiffusion3 / RC-Foundry) for reference and inference.

## Projects

| Path                | Contents                                                                                    | Status    | Link |
| ------------------- | ------------------------------------------------------------------------------------------- | --------- | ---- |
| **Hotspot sensitivity** | Does hotspot choice change design quality?                                                  | completed |  [README](./HotspotSensitivity/README.md)    |
| `project-2/`        | Atomic constraints: ablation of atom-level motif constraints (geometry, clashes, diversity) |           |      |
| `project-3/`        | End-to-end design and validation pipeline: RFD3 → MPNN → RF3                                |           |      |
| `project-4/`        | Targeted protein design for one target family; ranked, validated candidates                 |           |      |
