# Cross-reference the two rankings (raw affinity vs. ligand efficiency) to
# find compounds that score well under BOTH criteria - a stronger signal of
# a genuine hit than either metric alone, since each has its own known bias
# (raw affinity favors large molecules; LE favors tiny ones).

import csv
from rdkit import Chem

RESULTS_CSV = "data/processed/docking_results.csv"
LIBRARY_CSV = "data/processed/approved_drugs_library.csv"
OUTPUT_CSV = "data/processed/robust_hits.csv"
TOP_N = 50
MIN_HEAVY_ATOMS = 15

library = {}
with open(LIBRARY_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        library[row["chembl_id"]] = row

results = []
with open(RESULTS_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["status"] != "ok" or not row["best_affinity_kcal_mol"]:
            continue
        lib_entry = library.get(row["chembl_id"])
        if lib_entry is None:
            continue
        mol = Chem.MolFromSmiles(lib_entry["smiles"])
        if mol is None:
            continue
        heavy_atoms = mol.GetNumHeavyAtoms()
        affinity = float(row["best_affinity_kcal_mol"])
        results.append({
            "chembl_id": row["chembl_id"], "name": lib_entry["name"],
            "affinity": affinity, "heavy_atoms": heavy_atoms,
            "le": -affinity / heavy_atoms,
        })

by_affinity = sorted(results, key=lambda r: r["affinity"])[:TOP_N]
by_le = sorted(
    [r for r in results if r["heavy_atoms"] >= MIN_HEAVY_ATOMS],
    key=lambda r: r["le"], reverse=True
)[:TOP_N]

ids_affinity = {r["chembl_id"] for r in by_affinity}
ids_le = {r["chembl_id"] for r in by_le}
overlap_ids = ids_affinity & ids_le

overlap = [r for r in results if r["chembl_id"] in overlap_ids]
overlap.sort(key=lambda r: r["affinity"])

print(f"Compostos no top {TOP_N} de afinidade bruta: {len(ids_affinity)}")
print(f"Compostos no top {TOP_N} de eficiencia (LE): {len(ids_le)}")
print(f"Compostos presentes nos DOIS rankings: {len(overlap)}\n")

for r in overlap:
    print(f"  {r['chembl_id']:15s} {r['name']:30s} "
          f"{r['affinity']:.2f} kcal/mol  LE={r['le']:.3f}")

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["chembl_id", "name", "affinity_kcal_mol", "heavy_atoms", "ligand_efficiency"])
    writer.writeheader()
    for r in overlap:
        writer.writerow({
            "chembl_id": r["chembl_id"], "name": r["name"],
            "affinity_kcal_mol": r["affinity"], "heavy_atoms": r["heavy_atoms"],
            "ligand_efficiency": round(r["le"], 3),
        })

print(f"\nSalvo em {OUTPUT_CSV}")