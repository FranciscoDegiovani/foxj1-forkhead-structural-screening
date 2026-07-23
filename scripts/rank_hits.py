# Rank docking results by predicted binding affinity (most negative = best),
# joining with the drug library manifest to show compound names instead of
# just ChEMBL IDs, and export a clean top-N hit list for the manuscript.

import csv

RESULTS_CSV = "data/processed/docking_results.csv"
LIBRARY_CSV = "data/processed/approved_drugs_library.csv"
OUTPUT_CSV = "data/processed/top_hits.csv"
TOP_N = 30

# Load compound names from the curated library manifest
names = {}
with open(LIBRARY_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        names[row["chembl_id"]] = row["name"]

# Load docking results, keep only successful runs with a valid score
results = []
with open(RESULTS_CSV, newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["status"] == "ok" and row["best_affinity_kcal_mol"]:
            results.append({
                "chembl_id": row["chembl_id"],
                "name": names.get(row["chembl_id"], "desconhecido"),
                "affinity": float(row["best_affinity_kcal_mol"]),
            })

# Sort by affinity ascending (more negative kcal/mol = stronger predicted binding)
results.sort(key=lambda r: r["affinity"])

print(f"Total de resultados validos: {len(results)}")
print(f"\nTop {TOP_N} candidatos por afinidade predita:\n")
for i, r in enumerate(results[:TOP_N], start=1):
    print(f"{i:2d}. {r['chembl_id']:15s} {r['name']:35s} {r['affinity']:.2f} kcal/mol")

with open(OUTPUT_CSV, "w", newline="", encoding="utf-8") as f:
    writer = csv.DictWriter(f, fieldnames=["rank", "chembl_id", "name", "affinity_kcal_mol"])
    writer.writeheader()
    for i, r in enumerate(results[:TOP_N], start=1):
        writer.writerow({
            "rank": i,
            "chembl_id": r["chembl_id"],
            "name": r["name"],
            "affinity_kcal_mol": r["affinity"],
        })

print(f"\nSalvo em {OUTPUT_CSV}")