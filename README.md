# Prospecção Computacional de Moduladores Estruturais do Domínio Forkhead do FOXJ1

> Projeto exploratório de bioinformática estrutural — modelagem, virtual screening e dinâmica molecular.
> Não utiliza dados da dissertação de mestrado do autor (confidenciais, sob aprovação de orientador).

## Motivação

A perda de expressão do fator de transcrição FOXJ1 tem sido associada a desdiferenciação
epitelial e pior prognóstico em neoplasias ginecológicas (WEIR et al., 2023; COCHRANE et al.,
2020; WALENTEK, 2016). FOXJ1 é um fator de transcrição da família forkhead, sem bolsão
hidrofóbico clássico e sem estrutura cristalográfica depositada — o que o torna um alvo
estruturalmente pouco explorado. Este projeto investiga, de forma exploratória, se existem
moduladores de pequenas moléculas capazes de interagir com o domínio forkhead do FOXJ1
de forma a estabilizar sua conformação funcional.

**Nota sobre escopo:** dado que a literatura aponta a perda de FOXJ1 nesses tumores como um
fenômeno predominantemente de expressão/epigenética (e não de mutação desestabilizante da
proteína, ao contrário do que se observa para TP53), este projeto é uma prova de conceito
metodológica em bioinformática estrutural — não uma alegação terapêutica.

## Alvo

- **Proteína:** FOXJ1 humano (UniProt Q92949)
- **Domínio de interesse:** domínio forkhead, aa ~105–215 (IPR047512)
- **Estrutura de partida:** AlphaFold DB (AF-Q92949-F1)
- **Motivo de DNA reconhecido (referência):** 5'-HWDTGTTTGTTTA-3' / 5'-KTTTGTTGTTKTW-3'

## Pipeline

1. Preparo da estrutura (truncagem do domínio, validação de pLDDT)
2. Detecção de cavidades/bolsões (fpocket / P2Rank)
3. Curadoria da biblioteca de compostos (reposicionamento/produtos naturais)
4. Virtual screening — docking com AutoDock Vina
5. Reclassificação de poses via modelo de ML
6. Dinâmica molecular (GROMACS) dos candidatos priorizados
7. Análise e relatório final

## Estrutura do repositório
data/raw/ # dados originais, sem alteração
data/processed/ # estruturas preparadas
notebooks/ # notebooks de análise
scripts/ # scripts reutilizáveis
results/ # saídas de cada etapa
docs/ # notas metodológicas e referências


## Status

 Em desenvolvimento — projeto de portfólio pessoal.

## Referências

- Weir et al., 2023
- Cochrane et al., 2020
- Walentek, 2016

## Licença

MIT
