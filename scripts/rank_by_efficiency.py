# Re-rank docking hits by Ligand Efficiency (LE = -affinity / heavy atom count)
# instead of raw affinity, to correct for the well-known size bias in
# empirical scoring functions like Vina's (larger molecules accumulate more
# favorable contacts regardless of true binding specificity).

import csv
from rdkit import Chem

RESULTS_CSV = "data/processed/docking_results.csv"
LIBRARY_CSV = "data/processed/approved_drugs_library.csv"
OUTPUT_CSV = "data/processed/top_hits_by_efficiency.csv"
TOP_N = 30

library = {}
with open(LIBRARY_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        library[row["chembl_id"]] = row

results = []
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
        affinity = float(row["best_affinity_kcal_mol"])
        le = -affinity / heavy_atoms
        results.append({
            "chembl_id": chembl_id,
            "name": lib_entry["name"],
            "affinity": affinity,
            "heavy_atoms": heavy_atoms,
            "le": le,
        })

# Exclude very small molecules (heavy_atoms < 15): LE becomes mathematically
# inflated and unreliable for tiny fragments/solvents/anesthetics, which tend
# to bind promiscuously to any hydrophobic pocket regardless of specificity.
MIN_HEAVY_ATOMS = 15
results = [r for r in results if r["heavy_atoms"] >= MIN_HEAVY_ATOMS]

results.sort(key=lambda r: r["le"], reverse=True)

print(f"Total avaliado: {len(results)}")
print(f"\nTop {TOP_N} por Ligand Efficiency (afinidade normalizada pelo tamanho):\n")
for i, r in enumerate(results[:TOP_N], start=1):
    print(f"{i:2d}. {r['chembl_id']:15s} {r['name']:30s} "
          f"LE={r['le']:.3f}  ({r['affinity']:.2f} kcal/mol, {r['heavy_atoms']} atomos pesados)")

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["rank", "chembl_id", "name", "affinity_kcal_mol", "heavy_atoms", "ligand_efficiency"])
    writer.writeheader()
    for i, r in enumerate(results[:TOP_N], start=1):
        writer.writerow({
            "rank": i, "chembl_id": r["chembl_id"], "name": r["name"],
            "affinity_kcal_mol": r["affinity"], "heavy_atoms": r["heavy_atoms"],
            "ligand_efficiency": round(r["le"], 3),
        })

print(f"\nSalvo em {OUTPUT_CSV}")