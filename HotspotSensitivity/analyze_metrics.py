import gzip
from pathlib import Path
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch
import matplotlib.pyplot as plt
import seaborn as sns

cif_parser = MMCIFParser(QUIET=True)
pdb_parser = PDBParser(QUIET=True)

def load_structure(filepath: Path):
    name = filepath.name
    if name.endswith(".cif.gz"):
        with gzip.open(filepath, "rt") as handle:
            return cif_parser.get_structure("design", handle)
    if name.endswith(".cif"):
        return cif_parser.get_structure("design", str(filepath))
    return pdb_parser.get_structure("design", str(filepath))

def evaluate_structure(filepath: Path, hotspot_res_ids: list[str]) -> dict:
    structure = load_structure(filepath)
    model = structure[0]
    chains = list(model.get_chains())

    target_chain = model["A"]
    binder_chain = [c for c in chains if c.id != "A"][0]

    target_atoms = [a for a in target_chain.get_atoms() if a.element != "H"]
    binder_atoms = [a for a in binder_chain.get_atoms() if a.element != "H"]

    ns = NeighborSearch(target_atoms + binder_atoms)
    clashes = 0
    interface_pairs = 0
    satisfied_hotspots = set()

    for b_atom in binder_atoms:
        for n_atom in ns.search(b_atom.coord, radius=4.5, level="A"):
            if n_atom.get_parent().get_parent().id == "A":
                dist = np.linalg.norm(b_atom.coord - n_atom.coord)
                if dist < 2.0:
                    clashes += 1
                if dist <= 4.5:
                    interface_pairs += 1
                    res_num = f"A{n_atom.get_parent().id[1]}"
                    if res_num in hotspot_res_ids:
                        satisfied_hotspots.add(res_num)

    hotspot_recall = (
        len(satisfied_hotspots) / len(hotspot_res_ids)
        if hotspot_res_ids else np.nan
    )
    ca_coords = np.array([r["CA"].coord for r in binder_chain if "CA" in r])

    return {
        "binder_length": len(ca_coords),
        "steric_clashes": clashes,
        "interface_contacts": interface_pairs,
        "hotspot_recall": hotspot_recall,
        "ca_coords": ca_coords,
    }

def compute_pairwise_diversity(ca_list: list[np.ndarray]) -> float:
    rmsds = []
    for i in range(len(ca_list)):
        for j in range(i + 1, len(ca_list)):
            min_len = min(len(ca_list[i]), len(ca_list[j]))
            c1 = ca_list[i][:min_len] - ca_list[i][:min_len].mean(axis=0)
            c2 = ca_list[j][:min_len] - ca_list[j][:min_len].mean(axis=0)
            H = c1.T @ c2
            U, _, Vt = np.linalg.svd(H)
            R = Vt.T @ U.T
            if np.linalg.det(R) < 0:
                Vt[-1, :] *= -1
                R = Vt.T @ U.T
            rmsd = np.sqrt(np.mean(np.sum((c1 - (c2 @ R)) ** 2, axis=1)))
            rmsds.append(rmsd)
    return float(np.mean(rmsds)) if rmsds else 0.0

def main():
    conditions = {
        "cond_a_canonical": ["A56", "A115", "A123"],
        "cond_b_alternative": ["A66", "A68", "A75"],
        "cond_c_control": [],
    }
    records = []
    diversity_summary = {}

    for cond, hotspots in conditions.items():
        base_dir = Path(f"test/results/structures/{cond}")
        files = sorted(
            list(base_dir.rglob("*.cif.gz"))
            + list(base_dir.rglob("*.cif"))
            + list(base_dir.rglob("*.pdb"))
        )
        print(f"[{cond}] Found {len(files)} structure files")
        ca_bank = []
        for f in files:
            m = evaluate_structure(f, hotspots)
            ca_bank.append(m.pop("ca_coords"))
            m["condition"] = cond
            m["design_id"] = f.name
            m["geometry_pass"] = (m["steric_clashes"] == 0) and (m["interface_contacts"] >= 30)
            records.append(m)
        diversity_summary[cond] = compute_pairwise_diversity(ca_bank)

    if not records:
        raise RuntimeError("No structure files found under test/results/structures/")

    df = pd.DataFrame(records)
    df.to_csv("test/results/summary_metrics.csv", index=False)

    summary = df.groupby("condition").agg(
        n_designs=("design_id", "count"),
        mean_clashes=("steric_clashes", "mean"),
        mean_contacts=("interface_contacts", "mean"),
        hotspot_recall=("hotspot_recall", "mean"),
        geometry_pass_rate=("geometry_pass", "mean"),
    )
    summary["mean_pairwise_rmsd_A"] = summary.index.map(diversity_summary)
    print("\n=== Hotspot Sensitivity Comparison ===")
    print(summary.round(3))

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.5))
    sns.boxplot(data=df, x="condition", y="interface_contacts", ax=axes[0])
    axes[0].set_title("Interface Atom Contacts (<= 4.5 Å)")

    sns.boxplot(data=df, x="condition", y="steric_clashes", ax=axes[1])
    axes[1].set_title("Inter-Chain Steric Clashes (< 2.0 Å)")

    sns.barplot(
        data=df,
        x="condition",
        y="geometry_pass",
        estimator=np.mean,
        errorbar=("ci", 95),
        ax=axes[2],
    )
    axes[2].set_title("Geometric Pass Rate (95% CI)")
    axes[2].set_ylim(0, 1)

    plt.tight_layout()
    plt.savefig("test/results/hotspot_sensitivity_comparison.png", dpi=300)
    print("\nSaved figure to test/results/hotspot_sensitivity_comparison.png")

if __name__ == "__main__":
    main()
