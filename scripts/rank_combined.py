# Combine the two ranking criteria (raw affinity and ligand efficiency) into
# a single balanced score, instead of a brittle fixed top-N intersection
# (which returned zero overlap - the two metrics pull toward opposite ends
# of the molecular size spectrum). Each compound gets a rank position under
# each criterion; the combined score is the sum of both rank positions
# (lower = better under both criteria simultaneously).

import csv
from rdkit import Chem

RESULTS_CSV = "data/processed/docking_results.csv"
LIBRARY_CSV = "data/processed/approved_drugs_library.csv"
OUTPUT_CSV = "data/processed/top_hits_combined.csv"
TOP_N = 30
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
        if heavy_atoms < MIN_HEAVY_ATOMS:
            continue  # exclude tiny/promiscuous-binder-prone fragments, as before
        affinity = float(row["best_affinity_kcal_mol"])
        results.append({
            "chembl_id": row["chembl_id"], "name": lib_entry["name"],
            "affinity": affinity, "heavy_atoms": heavy_atoms,
            "le": -affinity / heavy_atoms,
        })

# Rank by affinity (most negative = rank 1)
for rank, r in enumerate(sorted(results, key=lambda r: r["affinity"]), start=1):
    r["rank_affinity"] = rank

# Rank by ligand efficiency (highest = rank 1)
for rank, r in enumerate(sorted(results, key=lambda r: r["le"], reverse=True), start=1):
    r["rank_le"] = rank

for r in results:
    r["combined_rank_score"] = r["rank_affinity"] + r["rank_le"]

results.sort(key=lambda r: r["combined_rank_score"])

print(f"Total avaliado (>= {MIN_HEAVY_ATOMS} atomos pesados): {len(results)}\n")
print(f"Top {TOP_N} por ranking combinado (afinidade + eficiencia):\n")
for i, r in enumerate(results[:TOP_N], start=1):
    print(f"{i:2d}. {r['chembl_id']:15s} {r['name']:25s} "
          f"{r['affinity']:.2f} kcal/mol (rank {r['rank_affinity']:4d})  "
          f"LE={r['le']:.3f} (rank {r['rank_le']:4d})")

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    fieldnames = ["combined_rank", "chembl_id", "name", "affinity_kcal_mol",
                  "rank_affinity", "ligand_efficiency", "rank_le", "combined_rank_score"]
    writer = csv.DictWriter(f, fieldnames=fieldnames)
    writer.writeheader()
    for i, r in enumerate(results[:TOP_N], start=1):
        writer.writerow({
            "combined_rank": i, "chembl_id": r["chembl_id"], "name": r["name"],
            "affinity_kcal_mol": r["affinity"], "rank_affinity": r["rank_affinity"],
            "ligand_efficiency": round(r["le"], 3), "rank_le": r["rank_le"],
            "combined_rank_score": r["combined_rank_score"],
        })

print(f"\nSalvo em {OUTPUT_CSV}")