import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, NeighborSearch, PPBuilder
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_theme(style="whitegrid", context="notebook")
cif_parser = MMCIFParser(QUIET=True)
ppb = PPBuilder()

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


def load_cif(filepath: Path):
    if filepath.name.endswith(".cif.gz"):
        with gzip.open(filepath, "rt") as handle:
            return cif_parser.get_structure("design", handle)
    return cif_parser.get_structure("design", str(filepath))


def classify_ss(phi_deg: float, psi_deg: float) -> str:
    """Classifies secondary structure from Ramachandran (phi, psi) angles."""
    if phi_deg is None or psi_deg is None:
        return "Coil"
    # Right-handed alpha-helix region
    if -100 <= phi_deg <= -30 and -67 <= psi_deg <= -7:
        return "Helix"
    # Beta-strand region
    if -170 <= phi_deg <= -70 and (psi_deg >= 90 or psi_deg <= -120):
        return "Sheet"
    return "Coil"


def extract_cif_features(filepath: Path, condition: str):
    structure = load_cif(filepath)
    model = structure[0]
    chains = list(model.get_chains())

    target_chain = model["A"]
    binder_chain = [c for c in chains if c.id != "A"][0]

    target_atoms = [a for a in target_chain.get_atoms() if a.element != "H"]
    binder_atoms = [a for a in binder_chain.get_atoms() if a.element != "H"]

    # 1. Center of Mass (COM) coordinates
    target_com = np.mean([a.coord for a in target_atoms], axis=0)
    binder_com = np.mean([a.coord for a in binder_atoms], axis=0)
    rel_com = binder_com - target_com  # Centered on target COM

    # 2. Target residue coordinates (for spatial hotspot overlay)
    res_coords = {
        r.id[1]: r["CA"].coord - target_com
        for r in target_chain
        if "CA" in r
    }

    # 3. Interface contacts (<= 4.5 A) & steric clashes (< 2.0 A)
    ns = NeighborSearch(target_atoms + binder_atoms)
    clashes = 0
    contacts = 0
    contacted_target_residues = set()

    for b_atom in binder_atoms:
        for n_atom in ns.search(b_atom.coord, radius=4.5, level="A"):
            if n_atom.get_parent().get_parent().id == "A":
                dist = np.linalg.norm(b_atom.coord - n_atom.coord)
                if dist < 2.0:
                    clashes += 1
                if dist <= 4.5:
                    contacts += 1
                    contacted_target_residues.add(n_atom.get_parent().id[1])

    # 4. Secondary structure fractions from phi/psi angles
    ss_counts = {"Helix": 0, "Sheet": 0, "Coil": 0}
    for poly in ppb.build_peptides(binder_chain):
        for phi, psi in poly.get_phi_psi_list():
            phi_d = np.degrees(phi) if phi is not None else None
            psi_d = np.degrees(psi) if psi is not None else None
            ss_counts[classify_ss(phi_d, psi_d)] += 1

    total_res = max(1, sum(ss_counts.values()))
    binder_len = sum(1 for r in binder_chain if "CA" in r)

    summary_row = {
        "condition": condition,
        "file": filepath.name,
        "binder_length": binder_len,
        "com_x": rel_com[0],
        "com_y": rel_com[1],
        "com_z": rel_com[2],
        "interface_contacts": contacts,
        "steric_clashes": clashes,
        "helix_frac": ss_counts["Helix"] / total_res,
        "sheet_frac": ss_counts["Sheet"] / total_res,
        "coil_frac": ss_counts["Coil"] / total_res,
    }
    return summary_row, contacted_target_residues, res_coords


def main():
    rows = []
    contact_matrix = {cond: {} for cond in HOTSPOTS}
    cond_counts = {cond: 0 for cond in HOTSPOTS}
    ref_res_coords = None

    for cond in HOTSPOTS:
        struct_dir = BASE_DIR / "results" / "structures" / cond
        files = sorted(
            list(struct_dir.rglob("*.cif.gz")) + list(struct_dir.rglob("*.cif"))
        )
        cond_counts[cond] = len(files)
        for f in files:
            row, contacted_res, res_coords = extract_cif_features(f, cond)
            rows.append(row)
            if ref_res_coords is None:
                ref_res_coords = res_coords
            for r_id in contacted_res:
                contact_matrix[cond][r_id] = contact_matrix[cond].get(r_id, 0) + 1

    if not rows:
        raise RuntimeError(f"No .cif or .cif.gz files found under {BASE_DIR / 'results' / 'structures'}/")

    df = pd.DataFrame(rows)
    out_csv = BASE_DIR / "results" / "cif_extracted_features.csv"
    df.to_csv(out_csv, index=False)

    # Create 2x2 diagnostic figure
    fig, axes = plt.subplots(2, 2, figsize=(15, 11))

    # --- Panel A: Spatial Docking Map (Binder COM relative to Target COM) ---
    ax = axes[0, 0]
    if ref_res_coords:
        t_xyz = np.array(list(ref_res_coords.values()))
        ax.scatter(t_xyz[:, 0], t_xyz[:, 1], c="#d0d0d0", s=25, label="Target Chain A (CA)")
        for r_id in HOTSPOTS["cond_a_canonical"]:
            if r_id in ref_res_coords:
                c = ref_res_coords[r_id]
                ax.scatter(c[0], c[1], c="#1f77b4", marker="*", s=220, edgecolor="black", zorder=5)
                ax.text(c[0] + 0.8, c[1] + 0.8, f"A{r_id}", fontsize=8, color="#1f77b4", fontweight="bold")
        for r_id in HOTSPOTS["cond_b_alternative"]:
            if r_id in ref_res_coords:
                c = ref_res_coords[r_id]
                ax.scatter(c[0], c[1], c="#ff7f0e", marker="^", s=160, edgecolor="black", zorder=5)
                ax.text(c[0] + 0.8, c[1] + 0.8, f"A{r_id}", fontsize=8, color="#d95f02", fontweight="bold")

    for cond, grp in df.groupby("condition"):
        ax.scatter(
            grp["com_x"], grp["com_y"],
            c=COLORS[cond], label=f"{cond} (n={len(grp)})",
            alpha=0.85, s=75, edgecolor="white"
        )
    ax.set_title("A. Spatial Docking Distribution (Binder COM vs. Target)", fontweight="bold")
    ax.set_xlabel("Relative X Coordinate (Å)")
    ax.set_ylabel("Relative Y Coordinate (Å)")
    ax.legend(loc="best", fontsize=8)

    # --- Panel B: Target Residue Contact Frequency Heatmap ---
    ax = axes[0, 1]
    all_res = sorted(ref_res_coords.keys()) if ref_res_coords else list(range(18, 133))
    heatmap_data = []
    for cond in HOTSPOTS:
        n = max(1, cond_counts[cond])
        freqs = [contact_matrix[cond].get(r, 0) / n for r in all_res]
        heatmap_data.append(freqs)
    hm_df = pd.DataFrame(heatmap_data, index=list(HOTSPOTS.keys()), columns=all_res)
    sns.heatmap(hm_df, cmap="YlGnBu", cbar_kws={"label": "Contact Fraction"}, ax=ax)
    ax.set_title("B. Target Interface Footprint (Residue Contact Frequency)", fontweight="bold")
    ax.set_xlabel("Target Residue Number (Chain A)")
    ax.set_ylabel("")

    # --- Panel C: Secondary Structure Composition per Condition ---
    ax = axes[1, 0]
    ss_melted = df.melt(
        id_vars=["condition"],
        value_vars=["helix_frac", "sheet_frac", "coil_frac"],
        var_name="SS_Type",
        value_name="Fraction",
    )
    ss_melted["SS_Type"] = ss_melted["SS_Type"].map({
        "helix_frac": "α-Helix",
        "sheet_frac": "β-Sheet",
        "coil_frac": "Loop / Coil",
    })
    sns.barplot(data=ss_melted, x="condition", y="Fraction", hue="SS_Type", ax=ax)
    ax.set_title("C. Designed Binder Secondary Structure Composition", fontweight="bold")
    ax.set_xlabel("")
    ax.set_ylabel("Mean Residue Fraction")
    ax.set_ylim(0, 1)

    # --- Panel D: Interface Packing vs. Steric Clashes & Length ---
    ax = axes[1, 1]
    for cond, grp in df.groupby("condition"):
        ax.scatter(
            grp["interface_contacts"],
            grp["steric_clashes"],
            s=grp["binder_length"] * 1.5,
            c=COLORS[cond],
            alpha=0.8,
            edgecolor="black",
            label=cond,
        )
    ax.axvline(30, color="gray", linestyle="--", linewidth=1, label="Min Contact Threshold (30)")
    ax.axhline(0.5, color="red", linestyle=":", linewidth=1, label="Zero-Clash Cutoff")
    ax.set_title("D. Interface Packing vs. Clashes (Bubble Size = Binder Length)", fontweight="bold")
    ax.set_xlabel("Heavy-Atom Interface Contacts (≤ 4.5 Å)")
    ax.set_ylabel("Inter-Chain Steric Clashes (< 2.0 Å)")
    ax.legend(loc="upper right", fontsize=8)

    plt.tight_layout()
    out_img = BASE_DIR / "results" / "cif_structural_diagnostics.png"
    plt.savefig(out_img, dpi=300)
    print(f"Saved 4-panel CIF diagnostic plot to: {out_img}")
    print(f"Saved extracted per-structure table to: {out_csv}")


if __name__ == "__main__":
    main()
