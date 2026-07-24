# Diário de Bordo Completo do Pipeline — FOXJ1 Forkhead Structural Screening

Registro cronológico de todas as etapas técnicas do projeto, incluindo decisões de design, obstáculos encontrados e como foram resolvidos. Complementa o `manuscript_draft.md` (que documenta metodologia formal) com o processo real de execução.

---

## Fase 1 — Estrutura e validação (Google Colab)

- Download da estrutura predita do FOXJ1 humano via AlphaFold DB (UniProt Q92949, modelo v6).
- Truncagem do domínio forkhead (resíduos 105–215) com Biopython.
- Validação de confiança: pLDDT médio de 87,57 na região truncada — "confiante" na escala do AlphaFold.
- **Obstáculo:** primeira tentativa de download via `curl`/URL fixa falhou (AlphaFold DB havia migrado de v4 para v6 recentemente) — resolvido baixando manualmente pelo navegador.

## Fase 2 — Detecção de cavidades e mapeamento funcional (Colab)

- fpocket identificou 6 cavidades candidatas; Pocket 3 selecionada por maior druggability score (0,287), apesar de não ter o maior score geométrico.
- Alinhamento de sequência (Biopython PairwiseAligner) entre FOXJ1 e FOXN1 (PDB 6EL8) mapeou a hélice de reconhecimento de DNA para os resíduos N167/S168/R170/H171.
- Confirmado: nenhum resíduo da Pocket 3 sobrepõe a hélice de reconhecimento — consistente com hipótese de modulação alostérica.
- **Obstáculos técnicos recorrentes no Colab:** o ambiente tinha dois interpretadores Python distintos (shell `pip` vs. kernel), exigindo `{sys.executable} -m pip install` repetidamente; sessões reiniciavam sem aviso, perdendo variáveis em memória (mitigado com uma célula de "bootstrap" reexecutável).

## Fase 3 — Curadoria da biblioteca de ligantes (Colab)

- Consulta à API do ChEMBL (max_phase=4, moléculas pequenas): 3.475 → 3.311 com SMILES válido → 2.860 após filtro de massa molecular (100–600 Da).
- Conversão para 3D (RDKit ETKDG + MMFF94) e PDBQT (Meeko): 1.924 sucessos na primeira passada; 936 falhas por fragmentação de sal (ex. "Valproate Sodium"), das quais 911 recuperadas retendo o maior fragmento molecular.
- **Nota de reprodutibilidade:** os 911 recuperados não chegaram a ser sincronizados para o GitHub antes da sessão do Colab expirar — biblioteca de trabalho usada no docking ficou em 1.924 compostos; reprocessamento dos 911 pendente.

## Fase 4 — Virtual screening por docking molecular

- Preparo do receptor (mk_prepare_receptor/Meeko) e caixa de busca (22×22×22 Å, centrada na Pocket 3).
- **Migração de ambiente:** Colab → execução local (Windows/PowerShell) por instabilidade de sessão em execuções longas.
- AutoDock Vina (instalado via binário standalone, já que `pip install vina` falha no Windows por dependência de Boost não resolvida) rodado em lote paralelizado (8 processos, 2 CPUs cada) com checkpoint incremental em CSV.
- **Obstáculo relevante:** primeira rodada completa "silenciosamente" reportou 100% de sucesso falso — bug no script que não capturava stderr; corrigido, e a segunda rodada revelou a causa real de falhas (timeout por suspensão do PC durante execução, resolvido desabilitando hibernação).
- Resultado final: 1.910/1.924 sucesso (99,3%); 8 timeouts + 6 falhas de parametrização de boro.

## Fase 5 — Correção de vieses e seleção de candidatos

- Ranking por afinidade bruta: dominado por moléculas grandes (viés de tamanho conhecido de funções de pontuação empíricas).
- Ranking por eficiência de ligante (LE) sem filtro: dominado por fragmentos/anestésicos voláteis (viés oposto).
- Ranking combinado (posição em afinidade + posição em LE, com piso de 15 átomos pesados) selecionado como critério final.
- Candidatos prioritários: testosterona (andrógeno), norethynodrel (progestágeno), quinestrol (estrogênio, adicionado após checagem de literatura sobre regulação estrogênica de FOXJ1) e carbamazepina (controle não esteroide).
- Validação visual das 4 poses no PyMOL — todas confirmadas dentro do volume da Pocket 3.

## Fase 6 — Manuscrito (versão inicial)

- Rascunho estruturado em formato ABNT (Introdução, Objetivos, Métodos, Resultados parciais, Discussão, Referências) gerado em `.docx`, com 9 figuras embutidas e 4 tabelas de resultados reais.
- Referências verificadas individualmente (Weir 2023, Walentek 2016, Newman 2018, Okada 2004, Yang 2025, ChEMBL, fpocket, AutoDock Vina); uma citação (Cochrane et al., 2020) não pôde ser verificada e foi sinalizada como pendente.

## Fase 7 — Dinâmica molecular (em andamento)

Ambiente: WSL2 (Ubuntu 24.04 LTS), GROMACS 2023.3, AmberTools/acpype via conda (necessário após falhas em cascata de bibliotecas do sistema na instalação via pip: `libarpack2t64`, `libhdf4-0`, `libhdf5_hl` — resolvido migrando para conda, que empacota dependências isoladas).

**Obstáculo relevante de química:** a conversão direta do PDBQT de saída do Vina para mol2 (via OpenBabel) gerou valência incorreta em 3 dos 4 ligantes esteroides (contagem de elétrons ímpar) — o formato PDBQT usa tipos de átomo "unidos" (hidrogênios implícitos) que o OpenBabel nem sempre reconstrói corretamente em esqueletos complexos. Corrigido com RDKit (`AssignBondOrdersFromTemplate`), usando o SMILES correto do ChEMBL como molde de conectividade aplicado às coordenadas 3D da pose já dockada — fórmulas moleculares confirmadas (testosterona C19H28O2, norethynodrel C20H26O2, quinestrol C25H32O2).

### Candidato 1: Testosterona

| Etapa | Status | Resultado-chave |
|---|---|---|
| Preparo ligante (correção de valência + acpype) | OK | C19H28O2 confirmada |
| Preparo receptor (pdb2gmx, AMBER99SB-ILDN) | OK | 1806 átomos |
| Montagem do complexo | OK | 1855 átomos |
| Solvatação (TIP3P) + íons (0,15 M NaCl) | OK | 62143 átomos, neutralizado |
| Minimização de energia | OK | Convergiu, Epot = -1.027.263 kJ/mol |
| Equilíbrio NVT (100 ps, 300K) | OK | Temp média 299,79 K |
| Equilíbrio NPT (100 ps) | OK | Densidade 989,22 kg/m³, pressão média -1,33 bar |
| Produção MD | PENDENTE | |

### Candidatos 2–4 (carbamazepina, norethynodrel, quinestrol)

PENDENTE — repetir fases de montagem de sistema em diante.

