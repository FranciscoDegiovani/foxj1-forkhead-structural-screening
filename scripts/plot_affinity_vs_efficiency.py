# Scatter plot of all successful docking results: raw affinity (X) vs.
# ligand efficiency (Y). Highlights the top-30 combined-rank hits, and the
# 3 compounds selected for pose inspection / molecular dynamics follow-up.

import csv
import matplotlib.pyplot as plt
from rdkit import Chem

RESULTS_CSV = "data/processed/docking_results.csv"
LIBRARY_CSV = "data/processed/approved_drugs_library.csv"
COMBINED_CSV = "data/processed/top_hits_combined.csv"
OUTPUT_PNG = "results/figures/affinity_vs_efficiency_scatter.png"

MD_CANDIDATES = {
    "CHEMBL386630": "Testosterone",
    "CHEMBL108": "Carbamazepine",
    "CHEMBL1387": "Norethynodrel",
}

top30_ids = set()
with open(COMBINED_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        top30_ids.add(row["chembl_id"])

library = {}
with open(LIBRARY_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        library[row["chembl_id"]] = row

all_pts, top30_pts, md_pts = [], [], []

with open(RESULTS_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["status"] != "ok" or not row["best_affinity_kcal_mol"]:
            continue
        chembl_id = row["chembl_id"]
        lib_entry = library.get(chembl_id)
        if lib_entry is None:
            continue
        mol = Chem.MolFromSmiles(lib_entry["smiles"])
        if mol is None:
            continue
        heavy_atoms = mol.GetNumHeavyAtoms()
        if heavy_atoms == 0:
            continue
        affinity = float(row["best_affinity_kcal_mol"])
        le = -affinity / heavy_atoms

        all_pts.append((affinity, le))
        if chembl_id in top30_ids:
            top30_pts.append((affinity, le))
        if chembl_id in MD_CANDIDATES:
            md_pts.append((affinity, le, MD_CANDIDATES[chembl_id]))

fig, ax = plt.subplots(figsize=(10, 7))

ax.scatter(*zip(*all_pts), s=8, alpha=0.25, color="grey",
           label=f"Full library (n={len(all_pts)})")
ax.scatter(*zip(*top30_pts), s=35, alpha=0.85, color="orange", edgecolor="none",
           label=f"Top 30 combined-rank hits (n={len(top30_pts)})")

for affinity, le, name in md_pts:
    ax.scatter(affinity, le, s=150, color="crimson", zorder=5, edgecolor="black", linewidth=1.2)
    ax.annotate(name, (affinity, le), textcoords="offset points", xytext=(9, 9),
                fontsize=11, fontweight="bold")

ax.set_xlabel("Predicted binding affinity (kcal/mol)", fontsize=12)
ax.set_ylabel("Ligand Efficiency (LE = -affinity / heavy atoms)", fontsize=12)
ax.set_title("Virtual screening against FOXJ1 (Pocket 3): raw affinity vs. ligand efficiency", fontsize=13)
ax.invert_xaxis()
ax.legend(loc="upper left", fontsize=9)
ax.grid(alpha=0.2)

plt.tight_layout()
plt.savefig(OUTPUT_PNG, dpi=300)
print(f"Salvo em {OUTPUT_PNG}")