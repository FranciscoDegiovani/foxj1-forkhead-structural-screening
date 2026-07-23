> **RASCUNHO — VERSÃO PRELIMINAR (documento vivo)**
> Este manuscrito documenta um projeto de portfólio em andamento. Seções ou dados marcados como **[PENDENTE]** ainda não foram executados. A lista de referências deve ser conferida contra as fontes originais antes de qualquer submissão formal — uma citação específica (Cochrane et al., 2020) não pôde ser verificada de forma independente e está sinalizada abaixo.

---

# Prospecção Computacional de Moduladores Estruturais do Domínio Forkhead do FOXJ1: Modelagem, Detecção de Sítios e Triagem Virtual como Estudo Exploratório em Neoplasias Ginecológicas

**Autor:** Francisco Degiovani
**Repositório do projeto:** https://github.com/FranciscoDegiovani/foxj1-forkhead-structural-screening

---

## Resumo

A perda de expressão do fator de transcrição FOXJ1 tem sido associada a desdiferenciação epitelial e a desfechos clínicos distintos em neoplasias ginecológicas, com evidência sugerindo um papel protetor em carcinoma seroso de alto grau de ovário. FOXJ1 é um fator de transcrição da família forkhead sem estrutura cristalográfica depositada e sem sítio de ligação a pequenas moléculas caracterizado, o que o torna um alvo estruturalmente pouco explorado do ponto de vista farmacológico. Este trabalho descreve um pipeline exploratório de bioinformática estrutural para identificação computacional de candidatos a moduladores de pequenas moléculas do domínio forkhead do FOXJ1, compreendendo: obtenção e validação de um modelo estrutural via AlphaFold Protein Structure Database; truncagem do domínio funcional; detecção de cavidades passíveis de ligação por fpocket; localização da hélice de reconhecimento de DNA por alinhamento de sequência com a estrutura cristalográfica do FOXN1 (PDB 6EL8); curadoria de uma biblioteca de fármacos aprovados a partir do ChEMBL; e triagem virtual por docking molecular com AutoDock Vina. O domínio truncado (resíduos 105–215) apresentou confiança estrutural elevada (pLDDT médio de 87,57). Seis cavidades foram identificadas por fpocket; a cavidade de maior escore de druggability (Pocket 3) não se sobrepõe aos resíduos de contato direto com bases de DNA mapeados por homologia (N167, S168, R170, H171), sendo compatível com um mecanismo hipotético de modulação alostérica em vez de bloqueio direto da interface proteína-DNA. Uma biblioteca de 2.835 compostos aprovados foi curada e preparada em formato PDBQT. A etapa de triagem virtual por docking está **em andamento** [PENDENTE — conclusão]; dinâmica molecular dos candidatos prioritários e reclassificação por aprendizado de máquina permanecem como próximas etapas [PENDENTE]. Os resultados aqui apresentados têm caráter hipótese-gerador e não constituem evidência de eficácia farmacológica, dada a ausência de validação experimental.

**Palavras-chave:** FOXJ1. Domínio forkhead. AlphaFold. Docking molecular. Triagem virtual. Bioinformática estrutural. Câncer ginecológico.

---

## 1 Introdução

Os fatores de transcrição da família forkhead box (FOX) constituem uma das maiores famílias de fatores de transcrição em humanos, compartilhando um domínio de ligação a DNA de aproximadamente 100 aminoácidos — o domínio forkhead (FH), uma subclasse da superfamília winged-helix (NEWMAN et al., 2018). A especificidade de reconhecimento de DNA nessa família é mediada predominantemente pela terceira α-hélice do domínio (hélice de reconhecimento, α3), que se insere no sulco maior da dupla-fita de DNA (NEWMAN et al., 2018).

FOXJ1 (Forkhead Box J1) é o regulador transcricional mestre do programa de ciliogênese motora, ativando a transcrição de genes necessários à montagem e função de cílios móveis em diversos tecidos epiteliais (STUBBS et al., 2008 apud NEWMAN et al., 2018). Para além de seu papel no desenvolvimento ciliar, a expressão de FOXJ1 tem sido investigada como biomarcador em contextos oncológicos, com achados que sugerem um papel context-dependente: em carcinoma seroso de alto grau de ovário, maior expressão proteica de FOXJ1 associou-se a melhor sobrevida global, num estudo de grande coorte do Ovarian Tumor Tissue Analysis Consortium (WEIR et al., 2023). Revisões sobre fatores de transcrição ciliares em câncer também descrevem a perda de expressão de FOXJ1 como característica de neoplasias mais agressivas em outros contextos, como ependimomas e tumores de plexo coroide (WALENTEK, 2016).

É relevante notar que a literatura recente também documenta um papel oposto de FOXJ1 em determinados contextos tumorais — por exemplo, sua participação na resistência a taxanos em linhagens de câncer resistentes à quimioterapia, onde a supressão farmacológica de FOXJ1 restaura a sensibilidade ao paclitaxel. Isso reforça que a farmacologia de FOXJ1 é dependente de tecido e contexto tumoral, e que a estratégia de restauração/estabilização de sua função — em vez de inibição — é justificada especificamente pelo cenário de desdiferenciação epitelial ginecológica, não como generalização para todos os tumores.

Do ponto de vista estrutural, FOXJ1 apresenta uma lacuna relevante: não existe estrutura cristalográfica depositada no Protein Data Bank para o domínio forkhead humano, ao contrário de diversos parálogos da família (FOXN1, FOXO1, FOXO3, FOXO4, FOXA2, FOXK1, FOXP2, FOXM1), todos com estruturas resolvidas em complexo com DNA. Adicionalmente, por se tratar de um fator de transcrição, FOXJ1 não possui um bolsão hidrofóbico catalítico clássico como alvos enzimáticos tradicionais, tornando a prospecção de moduladores de pequenas moléculas uma tarefa estruturalmente mais especulativa e exploratória do que em alvos já validados farmacologicamente.

Ainda assim, há precedente metodológico direto para essa abordagem dentro da própria família FOX: um estudo de triagem virtual estrutural identificou candidatos a inibidores do FOXM1 — outro membro da família forkhead sem bolsão catalítico clássico — como possíveis agentes terapêuticos em câncer de ovário, empregando docking com AutoDock Vina contra a estrutura da proteína isolada, sem necessidade de modelar o complexo com DNA. De forma semelhante, trabalhos recentes de modelagem estrutural via AlphaFold têm sido aplicados com sucesso a outros membros da família FOX com arquitetura semelhante (domínio bem enovelado seguido por regiões intrinsecamente desordenadas), como no caso do FOXP2.

O advento de métodos de predição de estrutura proteica baseados em aprendizado profundo, notadamente o AlphaFold (JUMPER et al., 2021), e sua disponibilização em larga escala através do AlphaFold Protein Structure Database (BERTONI et al., 2026), tornou viável a obtenção de um modelo estrutural de alta confiança para o domínio forkhead do FOXJ1 mesmo na ausência de dados experimentais, desde que as limitações inerentes a esse tipo de predição — sobretudo em regiões flexíveis e na ausência de conformação induzida por ligante — sejam adequadamente consideradas na interpretação dos resultados.

Este trabalho não deriva de nem reutiliza dados da dissertação de mestrado do autor (de natureza confidencial, sujeitos à aprovação do orientador), constituindo um estudo estrutural exploratório independente, ancorado exclusivamente em dados públicos.

---

## 2 Objetivos

### 2.1 Objetivo geral

Desenvolver e documentar um pipeline computacional de bioinformática estrutural para a identificação exploratória de candidatos a moduladores de pequenas moléculas do domínio forkhead do FOXJ1 humano, como estudo de prova de conceito metodológica no contexto de neoplasias ginecológicas.

### 2.2 Objetivos específicos

a) Obter e validar computacionalmente um modelo estrutural do domínio forkhead do FOXJ1 humano;
b) Identificar e caracterizar cavidades candidatas a sítios de ligação de pequenas moléculas na estrutura do domínio;
c) Localizar, por homologia estrutural/sequencial com um membro da família FOX de estrutura resolvida, a posição da hélice de reconhecimento de DNA, permitindo a interpretação mecanística das cavidades identificadas;
d) Curar uma biblioteca de compostos com estratégia de reposicionamento farmacológico (fármacos já aprovados clinicamente);
e) Executar triagem virtual por docking molecular dos compostos curados contra a cavidade selecionada;
f) [PENDENTE] Refinar os candidatos prioritários por dinâmica molecular;
g) [PENDENTE] Aplicar modelo de aprendizado de máquina para reclassificação das poses de docking.

---

## 3 Material e Métodos

### 3.1 Obtenção e validação da estrutura do domínio forkhead

A sequência e o modelo estrutural preditos do FOXJ1 humano (UniProt Q92949, isoforma canônica, 421 resíduos) foram obtidos do AlphaFold Protein Structure Database, entrada AF-Q92949-F1 (JUMPER et al., 2021; BERTONI et al., 2026), em formato PDB (versão de modelo v6). O domínio forkhead foi delimitado computacionalmente entre os resíduos 105 e 215, com margem de segurança sobre o intervalo funcional anotado (aproximadamente 120–210) descrito em bancos de anotação de domínio (InterPro IPR047512) e em material técnico de anticorpos comerciais validados contra a região. A confiança da predição foi avaliada por meio do escore pLDDT (predicted Local Distance Difference Test), extraído da coluna de fator-B do arquivo PDB, seguindo a convenção do AlphaFold DB.

O domínio foi truncado computacionalmente com a biblioteca Biopython (`Bio.PDB`), preservando apenas os resíduos no intervalo definido, e o pLDDT médio da região truncada foi calculado como medida-resumo de confiança estrutural.

### 3.2 Detecção de cavidades

A estrutura truncada foi submetida à ferramenta fpocket para detecção geométrica de cavidades passíveis de ligação de pequenas moléculas, com parâmetros padrão. Foram registradas, para cada cavidade identificada, as métricas nativas da ferramenta: escore geométrico (Score), escore de druggability (Druggability Score), volume, número de esferas alfa, proporção apolar, e escores de hidrofobicidade e polaridade.

### 3.3 Mapeamento da hélice de reconhecimento por homologia

Dado que FOXJ1 não possui estrutura resolvida em complexo com DNA, a localização de sua hélice de reconhecimento (α3) foi inferida por alinhamento de sequência par-a-par (algoritmo global, `Bio.Align.PairwiseAligner`, Biopython) entre o domínio forkhead truncado do FOXJ1 e o domínio forkhead do FOXN1 humano (UniProt O15353, resíduos 270–366), para o qual há estrutura cristalográfica resolvida em complexo com DNA (PDB 6EL8; NEWMAN et al., 2018). Os resíduos de contato direto com bases nitrogenadas descritos para o FOXN1 (N317, S318, R320, H321) foram mapeados, via alinhamento, para as posições equivalentes na sequência do FOXJ1, permitindo comparação posicional com as cavidades identificadas na etapa anterior.

### 3.4 Curadoria da biblioteca de ligantes

Compostos aprovados clinicamente foram obtidos programaticamente via API do ChEMBL (`chembl_webresource_client`), filtrando-se entradas classificadas como `molecule_type = "Small molecule"` e `max_phase = 4` (fase máxima de desenvolvimento clínico, correspondente à aprovação para uso humano). Entradas sem estrutura SMILES válida foram descartadas. A biblioteca foi adicionalmente filtrada por massa molecular (100–600 Da), critério padrão de similaridade a fármacos (drug-likeness), e por remoção de entradas contendo átomos metálicos (sais e agentes de contraste sem relevância como ligantes orgânicos de sítio proteico). Estruturas em forma salina, cujo SMILES continha múltiplos fragmentos desconectados (ex.: fármaco + contraíon sódio/citrato/cloridrato), foram recuperadas mediante retenção exclusiva do maior fragmento molecular, seguida de nova tentativa de conversão.

### 3.5 Preparo dos ligantes

Para cada composto curado, uma conformação tridimensional inicial foi gerada com o algoritmo ETKDG (RDKit), seguida de otimização de geometria pelo campo de força MMFF94. As estruturas 3D foram convertidas para o formato PDBQT — necessário para docking com AutoDock Vina — utilizando o Meeko, ferramenta oficial de preparo de ligantes do projeto AutoDock.

### 3.6 Preparo do receptor e definição da caixa de busca

O domínio truncado do FOXJ1 foi convertido para formato PDBQT com a ferramenta `mk_prepare_receptor` (Meeko), tratado como receptor rígido (docking rígido). O centro geométrico da cavidade selecionada (Pocket 3, ver Resultados) foi calculado como a média das coordenadas cartesianas de seus átomos de esfera alfa, definindo o centro de uma caixa de busca cúbica de 22 × 22 × 22 Å.

### 3.7 Triagem virtual por docking molecular

O docking molecular foi conduzido com AutoDock Vina (TROTT; OLSON, 2010; EBERHARDT et al., 2021), versão 1.2.7, função de pontuação padrão (`vina`), exaustividade de busca 8 e execução em lote paralelizada por processo (`ProcessPoolExecutor`, Python), com checkpoint incremental de resultados em arquivo CSV para permitir interrupção e retomada da triagem sem perda de progresso. [PENDENTE — execução em andamento; resultados completos a serem reportados em versão futura deste manuscrito]

### 3.8 [PENDENTE] Dinâmica molecular

Está prevista, como etapa subsequente, a simulação de dinâmica molecular (GROMACS) dos complexos formados pelos candidatos mais bem pontuados no docking, com o objetivo de avaliar a estabilidade temporal da interação proteína-ligante além da pontuação estática fornecida pelo Vina.

### 3.9 [PENDENTE] Reclassificação por aprendizado de máquina

Está prevista a aplicação de um modelo de aprendizado de máquina para reclassificação (rescoring) das poses de docking, como camada adicional de priorização dos candidatos, complementar à função de pontuação nativa do Vina.

### 3.10 Ambiente computacional e disponibilidade de dados

O pipeline foi desenvolvido em ambiente híbrido: prototipagem e etapas iniciais (obtenção de estrutura, truncagem, detecção de cavidades, alinhamento de sequência, curadoria e preparo da biblioteca de ligantes) em Google Colaboratory; a etapa de triagem virtual em lote foi migrada para execução local, por questões de estabilidade de sessão e de tempo total de processamento incompatível com os limites de uso gratuito da plataforma em nuvem. Todo o código, dados intermediários e resultados estão disponíveis publicamente em repositório versionado no GitHub: https://github.com/FranciscoDegiovani/foxj1-forkhead-structural-screening, sob licença MIT.

---

## 4 Resultados (parciais)

### 4.1 Estrutura do domínio forkhead do FOXJ1

O modelo estrutural do domínio forkhead do FOXJ1 (resíduos 105–215, 111 resíduos após truncagem) apresentou pLDDT médio de **87,57**, valor classificado como "confiante" na escala do próprio AlphaFold (70–90) e substancialmente superior à confiança observada nas regiões N- e C-terminais flanqueadoras da proteína completa — um padrão consistente com a arquitetura esperada de domínio bem enovelado seguido por regiões intrinsecamente desordenadas, típica de fatores de transcrição da família FOX.

### 4.2 Cavidades identificadas

A análise por fpocket identificou seis cavidades candidatas na estrutura do domínio (Tabela 1). A cavidade de maior escore de druggability (Pocket 3, Druggability Score = 0,287) não correspondeu à cavidade de maior escore geométrico (Pocket 1, Score = 0,377), evidenciando a dissociação, já esperada para um alvo sem bolsão catalítico clássico, entre "cavidade bem definida geometricamente" e "cavidade quimicamente compatível com ligante do tipo fármaco".

**Tabela 1** — Cavidades identificadas por fpocket na estrutura do domínio forkhead do FOXJ1.

| Cavidade | Score | Druggability Score | Volume (Å³) | Nº esferas α |
|---|---|---|---|---|
| Pocket 1 | 0,377 | 0,010 | 284,12 | 24 |
| Pocket 2 | 0,356 | 0,020 | 326,51 | 22 |
| Pocket 3 | 0,266 | **0,287** | 210,35 | 19 |
| Pocket 4 | 0,206 | 0,000 | 327,74 | 18 |
| Pocket 5 | 0,170 | 0,008 | 274,31 | 16 |
| Pocket 6 | 0,165 | 0,001 | 375,42 | 15 |

*Fonte: elaborado pelo autor a partir da saída do fpocket.*

### 4.3 Relação espacial com a hélice de reconhecimento de DNA

O alinhamento de sequência entre o domínio forkhead do FOXJ1 e o domínio forkhead do FOXN1 (identidade elevada na região central do alinhamento, compreendendo o motivo `WxNSIRHNLSLNKCF`, altamente conservado na família FOX) permitiu mapear os resíduos de contato direto com bases de DNA descritos para o FOXN1 (N317, S318, R320, H321) para as posições equivalentes **N167, S168, R170 e H171** no FOXJ1.

Nenhum resíduo das cavidades Pocket 1 (134, 135, 139–142, 195–197, 202, 206) ou Pocket 3 (133, 136, 137, 142, 146, 149, 150, 153) sobrepõe-se aos resíduos de leitura direta de base mapeados (167–171). A Pocket 3 localiza-se imediatamente anterior, na sequência primária, à hélice de reconhecimento — posição estruturalmente compatível com a face oposta ou adjacente ao feixe helicoidal α1/α2 que sustenta a α3, e não com a própria interface de leitura do DNA. Essa observação é consistente com a hipótese de que um ligante nessa cavidade atuaria por modulação alostérica da conformação/estabilidade da hélice de reconhecimento, e não por bloqueio competitivo direto da interação proteína-DNA — mecanismo estruturalmente mais plausível para um fator de transcrição do que a inibição competitiva clássica.

### 4.4 Biblioteca de ligantes curada

A consulta ao ChEMBL (`max_phase = 4`, moléculas pequenas) retornou 3.475 entradas, das quais 3.311 continham SMILES válido. Após filtro de massa molecular (100–600 Da), restaram 2.860 compostos. A conversão inicial para PDBQT (Meeko) obteve sucesso em 1.924 compostos; a recuperação de estruturas em forma salina (retenção do maior fragmento molecular) elevou o total para **2.835 compostos** processados com sucesso, restando 25 falhas residuais atribuídas a limitações de parametrização para átomos metálicos e estados de carga incomuns.

**Nota de reprodutibilidade:** em razão de uma interrupção de sessão computacional durante o desenvolvimento, a etapa de triagem virtual em lote reportada neste manuscrito foi executada, até o momento, sobre o subconjunto de 1.924 compostos originalmente convertidos, sincronizado com o ambiente de execução local. A execução de um segundo lote com os 911 compostos recuperados por remoção de contraíon está pendente [PENDENTE] e será incorporada em versão futura deste documento.

### 4.5 Triagem virtual por docking molecular

A execução da triagem virtual está em andamento no momento da elaboração desta versão do manuscrito, com taxa de processamento observada de aproximadamente 14 ligantes/minuto em ambiente local (16 núcleos de CPU, 8 processos paralelos, 2 núcleos por processo), projetando um tempo total estimado da ordem de poucas horas para o lote de 1.924 compostos. [PENDENTE — resultados quantitativos de afinidade preditos, ranking de candidatos e seleção de hits prioritários serão reportados em versão subsequente]

---

## 5 Discussão (parcial)

Os resultados estruturais obtidos até o momento sustentam a viabilidade metodológica da abordagem proposta. A confiança elevada do modelo AlphaFold especificamente na região do domínio forkhead — em contraste com a baixa confiança nas regiões flanqueadoras desordenadas — está alinhada ao comportamento esperado e documentado na literatura para predições de fatores de transcrição com arquitetura de domínio único bem enovelado. A ausência de sobreposição entre as cavidades detectadas e os resíduos de leitura direta de DNA é um achado que fortalece — sem comprovar — a plausibilidade mecanística de uma estratégia de modulação alostérica, e distingue este trabalho de uma tentativa ingênua de bloqueio competitivo do sítio de ligação ao DNA, que seria mecanisticamente menos defensável.

É necessário, contudo, enquadrar adequadamente o escopo interpretativo destes resultados. Em primeiro lugar, o modelo estrutural utilizado corresponde à conformação apo (não ligada a DNA) predita computacionalmente; não há garantia de que a cavidade identificada permaneça topologicamente acessível, ou mantenha a mesma geometria, na conformação biologicamente ativa da proteína. Em segundo lugar, a literatura recente de avaliação de modelos AlphaFold para docking de pequenas moléculas indica desempenho heterogêneo, particularmente em regiões flexíveis ou não representadas de forma robusta no conjunto de treinamento do modelo — uma ressalva diretamente aplicável a um alvo, como FOXJ1, desprovido de bolsão catalítico canônico. Em terceiro lugar, a ausência de mutações desestabilizantes recorrentes documentadas para FOXJ1 nos tumores de interesse (ao contrário do que se observa para TP53, com hotspots bem caracterizados e reativadores estruturais validados como controle positivo, a exemplo do APR-246) torna a narrativa de "correção estrutural" mais especulativa para FOXJ1 do que para alvos com paradigma de reativação já estabelecido.

Este trabalho deve, portanto, ser interpretado como um estudo de prova de conceito metodológica em bioinformática estrutural — hipótese-gerador, exploratório e não confirmatório — e não como evidência de eficácia terapêutica ou de existência de um modulador farmacológico real do FOXJ1.

---

## 6 Considerações Finais e Próximas Etapas

O pipeline desenvolvido até o presente momento demonstra a aplicabilidade de um fluxo de trabalho integrado — modelagem estrutural por AlphaFold, detecção de cavidades, mapeamento funcional por homologia e triagem virtual por docking — a um alvo estruturalmente não caracterizado da família forkhead. As seguintes etapas permanecem pendentes para a conclusão do estudo:

a) Conclusão da triagem virtual do lote de 1.924 compostos em execução;
b) Execução de segundo lote de docking para os 911 compostos recuperados por remoção de contraíon salino, totalizando a biblioteca completa de 2.835 candidatos;
c) Seleção dos candidatos com melhor pontuação preditiva de afinidade (hits prioritários) para refinamento;
d) Simulação de dinâmica molecular dos complexos formados pelos hits prioritários, para avaliação de estabilidade temporal da interação;
e) Reclassificação das poses de docking por modelo de aprendizado de máquina;
f) Recomenda-se, como validação metodológica adicional não prevista no escopo original, um ensaio de redocking contra um sistema de controle positivo com afinidade experimentalmente conhecida (por exemplo, no âmbito da Rota B do projeto original — reativadores estruturais de p53 mutante, com o composto APR-246 como referência), a fim de calibrar e validar a robustez do protocolo de docking empregado;
g) Validação experimental in vitro está fora do escopo deste projeto de portfólio e não é reivindicada em nenhuma etapa.

---

## Disponibilidade de Dados e Código

Todo o código-fonte, dados intermediários, resultados e o presente manuscrito em desenvolvimento estão disponíveis publicamente em: https://github.com/FranciscoDegiovani/foxj1-forkhead-structural-screening (licença MIT).

---

## Referências

BERTONI, D. et al. AlphaFold Protein Structure Database 2025: a redesigned interface and updated structural coverage. **Nucleic Acids Research**, v. 54, n. D1, p. D358-D362, 2026. DOI: 10.1093/nar/gkaf1226.

COCHRANE, D. R. et al. [Citação não verificada de forma independente pelo autor deste rascunho — conferir referência completa (periódico, volume, páginas, DOI) diretamente na lista de referências da dissertação de mestrado do autor, fonte original desta citação]. 2020.

EBERHARDT, J. et al. AutoDock Vina 1.2.0: New Docking Methods, Expanded Force Field, and Python Bindings. **Journal of Chemical Information and Modeling**, v. 61, n. 8, p. 3891-3898, 2021. DOI: 10.1021/acs.jcim.1c00203.

JUMPER, J. et al. Highly accurate protein structure prediction with AlphaFold. **Nature**, v. 596, p. 583-589, 2021. DOI: 10.1038/s41586-021-03819-2.

NEWMAN, J. A. et al. The structural basis for forkhead box family specificity revealed by the crystal structure of human FOXN1 in complex with DNA. **bioRxiv**, 2018. DOI: 10.1101/428011.

TROTT, O.; OLSON, A. J. AutoDock Vina: improving the speed and accuracy of docking with a new scoring function, efficient optimization, and multithreading. **Journal of Computational Chemistry**, v. 31, n. 2, p. 455-461, 2010. DOI: 10.1002/jcc.21334.

WALENTEK, P. Ciliary transcription factors in cancer – how understanding ciliogenesis can promote the detection and prognosis of cancer types. **The Journal of Pathology**, 2016. DOI: 10.1002/path.4703.

WEIR, A. et al. Increased FOXJ1 protein expression is associated with improved overall survival in high-grade serous ovarian carcinoma: an Ovarian Tumor Tissue Analysis Consortium Study. **British Journal of Cancer**, v. 128, p. 137-147, 2023. DOI: 10.1038/s41416-022-02014-y.

*Nota: referências relativas às ferramentas fpocket, ChEMBL, RDKit, Meeko, Biopython e InterPro devem ser adicionadas em revisão futura deste rascunho, com verificação bibliográfica completa antes de qualquer submissão.*
