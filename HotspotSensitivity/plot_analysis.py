import gzip
import json
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", context="notebook")
cif_parser = MMCIFParser(QUIET=True)
pdb_parser = PDBParser(QUIET=True)

BASE_DIR = Path("HotspotSensitivity")

HOTSPOTS = {
    "cond_a_canonical": [56, 115, 123],
    "cond_b_alternative": [66, 68, 75],
    "cond_c_control": [],
}

COLORS = {
    "cond_a_canonical": "#1f77b4",
    "cond_b_alternative": "#ff7f0e",
    "cond_c_control": "#2ca02c",
}

AROMATIC_RES = {"PHE", "TYR", "TRP", "HIS"}


def load_structure(filepath: Path):
    name = filepath.name
    if name.endswith(".cif.gz"):
        with gzip.open(filepath, "rt") as handle:
            return cif_parser.get_structure("design", handle)
    if name.endswith(".cif"):
        return cif_parser.get_structure("design", str(filepath))
    return pdb_parser.get_structure("design", str(filepath))


def compute_rg(ca_coords: np.ndarray) -> float:
    """Computes Radius of Gyration (Rg) in Angstroms from CA coordinates."""
    if len(ca_coords) == 0:
        return np.nan
    center = ca_coords.mean(axis=0)
    return float(np.sqrt(np.mean(np.sum((ca_coords - center) ** 2, axis=1))))


def kabsch_rmsd(c1: np.ndarray, c2: np.ndarray) -> float:
    """Computes Kabsch-aligned CA-RMSD on length-matched N-terminal cores."""
    n = min(len(c1), len(c2))
    p = c1[:n] - c1[:n].mean(axis=0)
    q = c2[:n] - c2[:n].mean(axis=0)
    H = p.T @ q
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return float(np.sqrt(np.mean(np.sum((p - (q @ R)) ** 2, axis=1))))


def analyze_single_complex(filepath: Path, condition: str):
    structure = load_structure(filepath)
    model = structure[0]
    chains = list(model.get_chains())

    target_chain = model["A"]
    binder_chain = [c for c in chains if c.id != "A"][0]

    target_atoms = [a for a in target_chain.get_atoms() if a.element != "H"]
    binder_atoms = [a for a in binder_chain.get_atoms() if a.element != "H"]
    binder_ca = np.array([r["CA"].coord for r in binder_chain if "CA" in r])

    # Hotspot atom coordinates for both patches
    canon_coords = np.array([
        a.coord for r in target_chain if r.id[1] in HOTSPOTS["cond_a_canonical"]
        for a in r if a.element != "H"
    ])
    alt_coords = np.array([
        a.coord for r in target_chain if r.id[1] in HOTSPOTS["cond_b_alternative"]
        for a in r if a.element != "H"
    ])
    binder_xyz = np.array([a.coord for a in binder_atoms])

    min_dist_canon = float(cdist(binder_xyz, canon_coords).min()) if len(canon_coords) else np.nan
    min_dist_alt = float(cdist(binder_xyz, alt_coords).min()) if len(alt_coords) else np.nan

    # Interface geometry, H-bonds (N/O <= 3.5 A), and chemical atom breakdown
    ns = NeighborSearch(target_atoms + binder_atoms)
    clashes = 0
    total_contacts = 0
    hbonds = 0
    hydrophobic_contacts = 0
    polar_contacts = 0
    aromatic_contacts = 0

    for b_atom in binder_atoms:
        b_elem = b_atom.element.upper()
        b_resname = b_atom.get_parent().resname.upper()
        for n_atom in ns.search(b_atom.coord, radius=4.5, level="A"):
            if n_atom.get_parent().get_parent().id == "A":
                dist = np.linalg.norm(b_atom.coord - n_atom.coord)
                n_elem = n_atom.element.upper()
                if dist < 2.0:
                    clashes += 1
                if dist <= 4.5:
                    total_contacts += 1
                    if b_resname in AROMATIC_RES:
                        aromatic_contacts += 1
                    elif b_elem in {"C", "S"}:
                        hydrophobic_contacts += 1
                    elif b_elem in {"N", "O"}:
                        polar_contacts += 1

                    # Putative H-bond / salt bridge (heavy N/O pair within 3.5 A)
                    if dist <= 3.5 and b_elem in {"N", "O"} and n_elem in {"N", "O"}:
                        hbonds += 1

    # Check assigned hotspot distance
    assigned_dist = (
        min_dist_canon if condition == "cond_a_canonical"
        else min_dist_alt if condition == "cond_b_alternative"
        else min(min_dist_canon, min_dist_alt)
    )

    row = {
        "condition": condition,
        "design_id": filepath.name,
        "binder_length": len(binder_ca),
        "rg_angstrom": compute_rg(binder_ca),
        "steric_clashes": clashes,
        "interface_contacts": total_contacts,
        "putative_hbonds": hbonds,
        "hydrophobic_contacts": hydrophobic_contacts,
        "polar_contacts": polar_contacts,
        "aromatic_contacts": aromatic_contacts,
        "min_dist_canonical_A": min_dist_canon,
        "min_dist_alternative_A": min_dist_alt,
        "assigned_hotspot_dist_A": assigned_dist,
        "pass_zero_clash": clashes == 0,
        "pass_packing": (clashes == 0) and (total_contacts >= 30),
        "pass_hotspot": (clashes == 0) and (total_contacts >= 30) and (assigned_dist <= 4.5),
    }
    return row, binder_ca


def check_rf3_confidence(condition: str):
    """Optionally reads RF3 confidence JSON files if Stage 3 refolding was run."""
    refold_dir = BASE_DIR / "results" / "refold" / condition
    if not refold_dir.exists():
        return []
    scores = []
    for jf in refold_dir.rglob("*.json"):
        try:
            data = json.loads(jf.read_text())
            if isinstance(data, dict):
                plddt = data.get("mean_plddt") or data.get("plddt")
                ipae = data.get("ipae") or data.get("interface_pae")
                if plddt is not None:
                    scores.append({"condition": condition, "file": jf.name, "plddt": plddt, "ipae": ipae})
        except Exception:
            continue
    return scores


def main():
    rows = []
    ca_banks = {cond: [] for cond in HOTSPOTS}
    rf3_records = []

    for cond in HOTSPOTS:
        struct_dir = BASE_DIR / "results" / "structures" / cond
        files = sorted(
            list(struct_dir.rglob("*.cif.gz"))
            + list(struct_dir.rglob("*.cif"))
            + list(struct_dir.rglob("*.pdb"))
        )
        for f in files:
            row, ca = analyze_single_complex(f, cond)
            rows.append(row)
            ca_banks[cond].append(ca)
        rf3_records.extend(check_rf3_confidence(cond))

    if not rows:
        raise RuntimeError(f"No structure files found under {BASE_DIR / 'results' / 'structures'}/")

    df = pd.DataFrame(rows)
    df.to_csv(BASE_DIR / "results" / "advanced_per_design_metrics.csv", index=False)

    # Compute pairwise CA-RMSD distributions per condition
    rmsd_rows = []
    for cond, cas in ca_banks.items():
        for i in range(len(cas)):
            for j in range(i + 1, len(cas)):
                rmsd_rows.append({"condition": cond, "pairwise_rmsd_A": kabsch_rmsd(cas[i], cas[j])})
    rmsd_df = pd.DataFrame(rmsd_rows)

    # Export aggregated statistical summary table
    agg = df.groupby("condition").agg(
        n_designs=("design_id", "count"),
        mean_length=("binder_length", "mean"),
        mean_rg_A=("rg_angstrom", "mean"),
        mean_clashes=("steric_clashes", "mean"),
        mean_contacts=("interface_contacts", "mean"),
        mean_hbonds=("putative_hbonds", "mean"),
        dist_to_canonical_A=("min_dist_canonical_A", "mean"),
        dist_to_alternative_A=("min_dist_alternative_A", "mean"),
        funnel_zero_clash_rate=("pass_zero_clash", "mean"),
        funnel_packing_rate=("pass_packing", "mean"),
        funnel_hotspot_rate=("pass_hotspot", "mean"),
    )
    if not rmsd_df.empty:
        agg["mean_pairwise_rmsd_A"] = rmsd_df.groupby("condition")["pairwise_rmsd_A"].mean()
    agg.to_csv(BASE_DIR / "results" / "advanced_summary_stats.csv")
    print("\n=== Advanced Structural & Physicochemical Summary ===")
    print(agg.round(3))

    # --- Build 3x2 Multi-Panel Figure ---
    fig, axes = plt.subplots(3, 2, figsize=(15, 16))

    # Panel A: Minimum Distance to Canonical vs. Alternative Hotspot Patches
    ax = axes[0, 0]
    dist_melted = df.melt(
        id_vars=["condition"],
        value_vars=["min_dist_canonical_A", "min_dist_alternative_A"],
        var_name="Target_Patch",
        value_name="Min_Distance_A",
    )
    dist_melted["Target_Patch"] = dist_melted["Target_Patch"].map({
        "min_dist_canonical_A": "Canonical Patch (A56/115/123)",
        "min_dist_alternative_A": "Alternative Patch (A66/68/75)",
    })
    sns.boxplot(data=dist_melted, x="condition", y="Min_Distance_A", hue="Target_Patch", ax=ax)
    ax.axhline(4.5, color="red", linestyle="--", linewidth=1.2, label="Direct Contact Cutoff (4.5 Å)")
    ax.set_title("A. Epitope Targeting: Distance to Canonical vs. Alternative Patches", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Minimum Heavy-Atom Distance (Å)")
    ax.legend(loc="best", fontsize=8)

    # Panel B: Interface Hydrogen Bonds & Polar Complementarity
    ax = axes[0, 1]
    sns.violinplot(data=df, x="condition", y="putative_hbonds", palette=COLORS, inner="point", ax=ax)
    ax.set_title("B. Putative Inter-Chain Hydrogen Bonds (N/O ≤ 3.5 Å)", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("H-Bond / Salt-Bridge Count")

    # Panel C: Binder Compactness (Radius of Gyration vs. Length)
    ax = axes[1, 0]
    for cond, grp in df.groupby("condition"):
        ax.scatter(
            grp["binder_length"], grp["rg_angstrom"],
            c=COLORS[cond], label=cond, s=70, alpha=0.85, edgecolor="black"
        )
    # Theoretical globular protein scaling reference: Rg ~ 2.2 * N^(0.38)
    n_grid = np.linspace(55, 95, 100)
    ax.plot(n_grid, 2.2 * (n_grid ** 0.38), "k--", linewidth=1.2, label="Globular Scaling (2.2·N^0.38)")
    ax.set_title("C. Binder Fold Compactness (Radius of Gyration vs. Length)", fontweight="bold")
    ax.set_xlabel("Designed Binder Length (Residues)")
    ax.set_ylabel("Radius of Gyration Rg (Å)")
    ax.legend(loc="best", fontsize=8)

    # Panel D: Chemical Makeup of Interface Contacts
    ax = axes[1, 1]
    chem_means = df.groupby("condition")[["hydrophobic_contacts", "polar_contacts", "aromatic_contacts"]].mean()
    chem_means.plot(kind="bar", stacked=True, color=["#4c72b0", "#dd8452", "#8172b3"], ax=ax)
    ax.set_title("D. Mean Interface Contact Chemistry by Condition", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Mean Heavy-Atom Contacts (≤ 4.5 Å)")
    ax.tick_params(axis="x", rotation=0)
    ax.legend(["Aliphatic/Sulfur (C/S)", "Polar (N/O)", "Aromatic (F/Y/W/H)"], fontsize=8)

    # Panel E: Pairwise Structural Fold Diversity (CA-RMSD Distribution)
    ax = axes[2, 0]
    if not rmsd_df.empty:
        sns.kdeplot(
            data=rmsd_df, x="pairwise_rmsd_A", hue="condition",
            palette=COLORS, fill=True, common_norm=False, alpha=0.35, ax=ax
        )
        ax.set_title("E. Intra-Condition Fold Diversity (Pairwise CA-RMSD)", fontweight="bold")
        ax.set_xlabel("Pairwise Kabsch CA-RMSD (Å)")
        ax.set_ylabel("Density")
    else:
        ax.text(0.5, 0.5, "Need ≥ 2 designs per condition\nfor pairwise RMSD distribution", ha="center", va="center")
        ax.set_title("E. Intra-Condition Fold Diversity", fontweight="bold")

    # Panel F: Multi-Stage Design Funnel Pass Rate
    ax = axes[2, 1]
    funnel_df = df.groupby("condition")[["pass_zero_clash", "pass_packing", "pass_hotspot"]].mean().reset_index()
    funnel_melted = funnel_df.melt(id_vars="condition", var_name="Gate", value_name="Pass_Rate")
    funnel_melted["Gate"] = funnel_melted["Gate"].map({
        "pass_zero_clash": "1. Clash-Free",
        "pass_packing": "2. ≥30 Contacts",
        "pass_hotspot": "3. Hotspot ≤4.5Å",
    })
    sns.barplot(data=funnel_melted, x="Gate", y="Pass_Rate", hue="condition", palette=COLORS, ax=ax)
    ax.set_title("F. Sequential In Silico Design Funnel Yield", fontweight="bold")
    ax.set_xlabel("Filtering Stage")
    ax.set_ylabel("Fraction of Designs Passing")
    ax.set_ylim(0, 1.05)
    ax.legend(loc="best", fontsize=8)

    plt.tight_layout()
    out_fig = BASE_DIR / "results" / "advanced_hotspot_analysis.png"
    plt.savefig(out_fig, dpi=300)
    print(f"\nSaved 6-panel advanced analysis figure to: {out_fig}")
    print(f"Saved aggregated statistics table to: {BASE_DIR / 'results' / 'advanced_summary_stats.csv'}")


if __name__ == "__main__":
    main()
