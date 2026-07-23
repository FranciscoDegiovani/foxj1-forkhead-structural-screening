from Bio.PDB import PDBParser, PDBIO, Select

INPUT_PDB = "data/raw/AF-Q92949-F1-model_v6.pdb"
OUTPUT_PDB = "data/processed/FOXJ1_forkhead_domain.pdb"
DOMAIN_START = 105
DOMAIN_END = 215

class DomainSelect(Select):
    def accept_residue(self, residue):
        resnum = residue.get_id()[1]
        return DOMAIN_START <= resnum <= DOMAIN_END

parser = PDBParser(QUIET=True)
structure = parser.get_structure("FOXJ1", INPUT_PDB)

plddts = []
for residue in structure[0]["A"]:
    resnum = residue.get_id()[1]
    if DOMAIN_START <= resnum <= DOMAIN_END:
        for atom in residue:
            plddts.append(atom.get_bfactor())
            break

mean_plddt = sum(plddts) / len(plddts)
print(f"Residuos no dominio: {len(plddts)}")
print(f"pLDDT medio do dominio forkhead (aa {DOMAIN_START}-{DOMAIN_END}): {mean_plddt:.2f}")

io = PDBIO()
io.set_structure(structure)
io.save(OUTPUT_PDB, DomainSelect())
print(f"Dominio truncado salvo em {OUTPUT_PDB}")
