import csv
from rdkit import Chem

CHEMBL_ID = "CHEMBL1201165"  # Quinestrol

library = {}
with open("data/processed/approved_drugs_library.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        library[row["chembl_id"]] = row

with open("data/processed/docking_results.csv", newline="", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if row["chembl_id"] == CHEMBL_ID:
            affinity = float(row["best_affinity_kcal_mol"])
            mol = Chem.MolFromSmiles(library[CHEMBL_ID]["smiles"])
            heavy = mol.GetNumHeavyAtoms()
            print(f"{CHEMBL_ID} ({library[CHEMBL_ID]['name']})")
            print(f"Afinidade: {affinity:.2f} kcal/mol")
            print(f"Atomos pesados: {heavy}")
            print(f"LE: {-affinity/heavy:.3f}")