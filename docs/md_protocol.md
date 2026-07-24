# Protocolo de Dinâmica Molecular — Diário de Bordo

Ambiente: WSL2 Ubuntu 24.04 LTS, GROMACS 2023.3, AmberTools/acpype via conda.
Campo de força: AMBER99SB-ILDN (proteína) + GAFF2 (ligantes, cargas AM1-BCC).

## Candidato 1: Testosterona

| Etapa | Status | Resultado-chave |
|---|---|---|
| Preparo ligante (correção de valência + acpype) | OK | Fórmula confirmada C19H28O2 |
| Preparo receptor (pdb2gmx, AMBER99SB-ILDN) | OK | 1806 átomos |
| Montagem do complexo | OK | 1855 átomos totais |
| Solvatação (TIP3P) + íons (0,15 M NaCl) | OK | 62143 átomos, neutralizado (57 Na+, 64 Cl-) |
| Minimização de energia | OK | Convergiu, Epot = -1.027.263 kJ/mol (2767 passos) |
| Equilíbrio NVT (100 ps, 300K) | OK | Temp média 299,79 K, RMSD 2,63 |
| Equilíbrio NPT (100 ps) | OK | [PENDENTE conferir pressão/densidade] |
| Produção MD | PENDENTE | |

