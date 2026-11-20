import csv
import gzip
import json
import os
import re
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RESULTS_DIR = Path("/Users/mhmdvjn/1.work/Uni/sem1/bio/projects/Protein-Generation/HotspotSensitivity/results")
OUT_DIR = Path("/Users/mhmdvjn/1.work/Uni/sem1/bio/projects/Protein-Generation/HotspotSensitivity/visualization")
OUT_DIR.mkdir(parents=True, exist_ok=True)

COND_ORDER = ["cond_a_canonical", "cond_b_alternative", "cond_c_control"]
COND_LABELS = {
    "cond_a_canonical": "Condition A\nCanonical",
    "cond_b_alternative": "Condition B\nAlternative",
    "cond_c_control": "Condition C\nControl",
}
COLORS = {"cond_a_canonical": "#4c78a8", "cond_b_alternative": "#f58518", "cond_c_control": "#54a24b"}


def load_summary_metrics():
    rows = []
    with open(RESULTS_DIR / "summary_metrics.csv", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["binder_length"] = int(row["binder_length"])
            row["steric_clashes"] = int(row["steric_clashes"])
            row["interface_contacts"] = int(row["interface_contacts"])
            row["hotspot_recall"] = float(row["hotspot_recall"]) if row["hotspot_recall"] else np.nan
            row["geometry_pass"] = row["geometry_pass"].lower() == "true"
            rows.append(row)
    return rows


def load_structure_jsons():
    records = []
    for json_path in sorted((RESULTS_DIR / "structures").rglob("*.json")):
        with open(json_path) as f:
            data = json.load(f)
        spec = data.get("specification", {})
        extra = spec.get("extra", {})
        hotspots = spec.get("select_hotspots", {}) or {}
        rec = {
            "design_id": extra.get("example_id", json_path.stem),
            "condition": extra.get("example", ""),
            "num_residues": data["metrics"].get("num_residues"),
            "max_ca_deviation": data["metrics"].get("max_ca_deviation"),
            "n_chainbreaks": data["metrics"].get("n_chainbreaks"),
            "n_clashing_backbone": data["metrics"].get("n_clashing.interresidue_clashes_w_backbone"),
            "n_clashing_sidechain": data["metrics"].get("n_clashing.interresidue_clashes_w_sidechain"),
            "loop_fraction": data["metrics"].get("loop_fraction"),
            "helix_fraction": data["metrics"].get("helix_fraction"),
            "sheet_fraction": data["metrics"].get("sheet_fraction"),
            "non_loop_fraction": data["metrics"].get("non_loop_fraction"),
            "radius_of_gyration": data["metrics"].get("radius_of_gyration"),
            "num_ss_elements": data["metrics"].get("num_ss_elements"),
            "alanine_content": data["metrics"].get("alanine_content"),
            "glycine_content": data["metrics"].get("glycine_content"),
            "num_hotspots": len(hotspots),
            "contig": spec.get("contig", ""),
        }
        records.append(rec)
    return records


def load_refold_scores():
    records = []
    for csv_path in sorted((RESULTS_DIR / "refold").rglob("*_ranking_scores.csv")):
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                m = re.search(r"(cond_[a-z]+)", csv_path.parts[-3] if len(csv_path.parts) >= 3 else "")
                condition = m.group(1) if m else ""
                rec = {
                    "design_id": csv_path.parent.name,
                    "condition": condition,
                    "early_stopped": row.get("early_stopped", "").lower() == "true",
                    "mean_plddt": float(row["mean_plddt"]) if row.get("mean_plddt") else np.nan,
                }
                records.append(rec)
    return records


def load_mpnn_fasta():
    records = []
    for fa_path in sorted((RESULTS_DIR / "mpnn_fasta").rglob("*.fa")):
        condition = fa_path.parts[-3] if len(fa_path.parts) >= 3 else ""
        with open(fa_path) as f:
            header = None
            seq_lines = []
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if header is not None:
                        records.append(parse_fasta_entry(header, "".join(seq_lines), condition))
                    header = line[1:]
                    seq_lines = []
                else:
                    seq_lines.append(line)
            if header is not None:
                records.append(parse_fasta_entry(header, "".join(seq_lines), condition))
    return records


def parse_fasta_entry(header, sequence, condition):
    parts = {}
    for token in header.split(","):
        token = token.strip()
        if "=" in token:
            k, v = token.split("=", 1)
            parts[k.strip()] = v.strip()
    return {
        "design_id": parts.get("", header.split(",")[0].split()[0]),
        "condition": condition,
        "confidence": float(parts.get("confidence", "nan")) if parts.get("confidence") else np.nan,
        "sequence_recovery": float(parts.get("sequence_recovery", "nan")) if parts.get("sequence_recovery") else np.nan,
        "sequence_length": len(sequence),
    }


def parse_cif_coordinates(cif_path):
    """Extract CA atom coordinates from a (possibly gzipped) CIF file."""
    opener = gzip.open if str(cif_path).endswith(".gz") else open
    coords = []
    in_atom_site = False
    atom_fields = None
    with opener(cif_path, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith("loop_"):
                in_atom_site = False
                continue
            if line.startswith("_atom_site."):
                in_atom_site = True
                atom_fields = (atom_fields or []) + [line.split(".")[1]]
                continue
            if in_atom_site and line and not line.startswith("_") and not line.startswith("#"):
                fields = line.split()
                if atom_fields and len(fields) >= len(atom_fields):
                    row = dict(zip(atom_fields, fields))
                    if row.get("label_atom_id") == "CA":
                        try:
                            coords.append([
                                float(row["Cartn_x"]),
                                float(row["Cartn_y"]),
                                float(row["Cartn_z"]),
                                row.get("label_comp_id", ""),
                                int(row.get("label_seq_id", 0)),
                            ])
                        except (ValueError, KeyError):
                            pass
    return coords


def save_fig(fig, name):
    path = OUT_DIR / name
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"saved {path}")


def plot_summary(metrics):
    by_cond = defaultdict(list)
    for r in metrics:
        by_cond[r["condition"]].append(r)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Summary Metrics by Condition", fontsize=16)

    metric_specs = [
        ("steric_clashes", "Steric Clashes", "count"),
        ("interface_contacts", "Interface Contacts", "count"),
        ("hotspot_recall", "Hotspot Recall", "fraction"),
        ("geometry_pass", "Geometry Pass Rate", "rate"),
        ("binder_length", "Binder Length", "residues"),
    ]

    for ax, (key, title, unit) in zip(axes.flat[:5], metric_specs):
        means = []
        stds = []
        conds = []
        for c in COND_ORDER:
            vals = [r[key] for r in by_cond.get(c, []) if r[key] is not None and not (isinstance(r[key], float) and np.isnan(r[key]))]
            if vals:
                conds.append(c)
                means.append(np.mean(vals))
                stds.append(np.std(vals))
        if conds:
            x = np.arange(len(conds))
            ax.bar(x, means, yerr=stds, color=[COLORS[c] for c in conds], capsize=4, alpha=0.85, edgecolor="black", linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([COND_LABELS[c].replace("\n", " ") for c in conds], rotation=15, ha="right", fontsize=8)
            ax.set_title(title)
            ax.set_ylabel(unit)
            ax.grid(axis="y", alpha=0.3)
        else:
            ax.set_visible(False)

    axes.flat[5].axis("off")
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_fig(fig, "summary_metrics_by_condition.png")


def plot_structure_metrics(records):
    by_cond = defaultdict(list)
    for r in records:
        by_cond[r["condition"]].append(r)

    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    fig.suptitle("Structure Quality Metrics by Condition", fontsize=16)

    metric_specs = [
        ("max_ca_deviation", "Max CA Deviation", "\u00c5"),
        ("n_chainbreaks", "Chain Breaks", "count"),
        ("n_clashing_backbone", "Backbone Clashes", "count"),
        ("n_clashing_sidechain", "Sidechain Clashes", "count"),
        ("helix_fraction", "Helix Fraction", "fraction"),
        ("sheet_fraction", "Sheet Fraction", "fraction"),
    ]

    for ax, (key, title, unit) in zip(axes.flat, metric_specs):
        means, stds, conds = [], [], []
        for c in COND_ORDER:
            vals = [r[key] for r in by_cond.get(c, []) if r.get(key) is not None and not (isinstance(r[key], float) and np.isnan(r[key]))]
            if vals:
                conds.append(c)
                means.append(np.mean(vals))
                stds.append(np.std(vals))
        if conds:
            x = np.arange(len(conds))
            ax.bar(x, means, yerr=stds, color=[COLORS[c] for c in conds], capsize=4, alpha=0.85, edgecolor="black", linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([COND_LABELS[c].replace("\n", " ") for c in conds], rotation=15, ha="right", fontsize=8)
            ax.set_title(title)
            ax.set_ylabel(unit)
            ax.grid(axis="y", alpha=0.3)
        else:
            ax.set_visible(False)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_fig(fig, "structure_metrics_by_condition.png")


def plot_refold_plddt(records):
    by_cond = defaultdict(list)
    for r in records:
        by_cond[r["condition"]].append(r)

    fig, ax = plt.subplots(figsize=(8, 6))
    data = []
    labels = []
    for c in COND_ORDER:
        vals = [r["mean_plddt"] for r in by_cond.get(c, []) if not np.isnan(r["mean_plddt"])]
        if vals:
            data.append(vals)
            labels.append(COND_LABELS[c].replace("\n", " "))
    if data:
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], [c for c in COND_ORDER if by_cond.get(c)]):
            patch.set_facecolor(COLORS[c])
            patch.set_alpha(0.6)
    ax.set_title("Refold pLDDT Distribution by Condition")
    ax.set_ylabel("mean pLDDT")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "refold_plddt.png")


def plot_mpnn_recovery(records):
    by_cond = defaultdict(list)
    for r in records:
        by_cond[r["condition"]].append(r)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    fig.suptitle("MPNN Sequence Design by Condition", fontsize=16)

    for ax, key, title, ylabel in [
        (axes[0], "sequence_recovery", "Sequence Recovery", "fraction"),
        (axes[1], "confidence", "Confidence", "score"),
    ]:
        means, stds, conds = [], [], []
        for c in COND_ORDER:
            vals = [r[key] for r in by_cond.get(c, []) if not np.isnan(r[key])]
            if vals:
                conds.append(c)
                means.append(np.mean(vals))
                stds.append(np.std(vals))
        if conds:
            x = np.arange(len(conds))
            ax.bar(x, means, yerr=stds, color=[COLORS[c] for c in conds], capsize=4, alpha=0.85, edgecolor="black", linewidth=0.5)
            ax.set_xticks(x)
            ax.set_xticklabels([COND_LABELS[c].replace("\n", " ") for c in conds], rotation=15, ha="right", fontsize=8)
            ax.set_title(title)
            ax.set_ylabel(ylabel)
            ax.grid(axis="y", alpha=0.3)
        else:
            ax.set_visible(False)

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_fig(fig, "mpnn_design_metrics.png")


def plot_hotspot_recall(metrics):
    by_cond = defaultdict(list)
    for r in metrics:
        by_cond[r["condition"]].append(r)

    fig, ax = plt.subplots(figsize=(8, 6))
    data, labels = [], []
    for c in COND_ORDER:
        vals = [r["hotspot_recall"] for r in by_cond.get(c, []) if not np.isnan(r["hotspot_recall"])]
        if vals:
            data.append(vals)
            labels.append(COND_LABELS[c].replace("\n", " "))
    if data:
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, medianprops=dict(color="black"))
        for patch, c in zip(bp["boxes"], [c for c in COND_ORDER if by_cond.get(c)]):
            patch.set_facecolor(COLORS[c])
            patch.set_alpha(0.6)
    ax.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    ax.set_title("Hotspot Recall Distribution by Condition")
    ax.set_ylabel("hotspot recall")
    ax.grid(axis="y", alpha=0.3)
    fig.tight_layout()
    save_fig(fig, "hotspot_recall.png")


def plot_3d_structures():
    """Plot CA trace 3D for one representative structure per condition."""
    struct_dir = RESULTS_DIR / "structures"
    reps = {}
    for cond in COND_ORDER:
        cond_dir = struct_dir / cond
        if not cond_dir.exists():
            continue
        cif_files = sorted(cond_dir.glob("*.cif.gz"))
        if cif_files:
            reps[cond] = cif_files[0]

    if not reps:
        print("No CIF files found for 3D visualization.")
        return

    fig = plt.figure(figsize=(15, 5))
    for i, (cond, path) in enumerate(reps.items()):
        ax = fig.add_subplot(1, len(reps), i + 1, projection="3d")
        coords_raw = parse_cif_coordinates(path)
        if not coords_raw:
            ax.text(0.5, 0.5, 0.5, "No CA coords", transform=ax.transAxes, ha="center")
            ax.set_title(f"{COND_LABELS[cond]}\n(no coords)")
            continue
        coords = np.array([c[:3] for c in coords_raw])
        seq_ids = np.array([c[4] for c in coords_raw])
        ax.plot(coords[:, 0], coords[:, 1], coords[:, 2], color=COLORS[cond], linewidth=2, alpha=0.8)
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=seq_ids, cmap="viridis", s=20, edgecolor="black", linewidth=0.3)
        ax.set_title(f"{COND_LABELS[cond].replace(chr(10), ' ')}\n{path.name}")
        ax.set_xlabel("X")
        ax.set_ylabel("Y")
        ax.set_zlabel("Z")
        ax.grid(alpha=0.3)

    fig.suptitle("Backbone CA Traces (Representative Structure per Condition)", fontsize=14)
    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    save_fig(fig, "structure_3d_traces.png")


def main():
    print("Loading data...")
    metrics = load_summary_metrics()
    struct_records = load_structure_jsons()
    refold_records = load_refold_scores()
    mpnn_records = load_mpnn_fasta()
    print(f"  summary metrics: {len(metrics)} rows")
    print(f"  structure jsons: {len(struct_records)} rows")
    print(f"  refold scores: {len(refold_records)} rows")
    print(f"  mpnn fasta entries: {len(mpnn_records)} rows")

    print("Generating figures...")
    plot_summary(metrics)
    plot_structure_metrics(struct_records)
    plot_refold_plddt(refold_records)
    plot_mpnn_recovery(mpnn_records)
    plot_hotspot_recall(metrics)
    plot_3d_structures()

    print(f"\nAll figures written to {OUT_DIR}")


if __name__ == "__main__":
    main()
