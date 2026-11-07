# Simulador do Pipeline DNA Fountain

[![Python Version](https://img.shields.io/badge/python-3.12%2B-blue.svg)](https://www.python.org/)
[![Test Suite](https://img.shields.io/badge/pytest-18%20passed-brightgreen.svg)](tests/)
[![Release](https://img.shields.io/badge/release-v1.0--defesa-orange.svg)](https://github.com/GabrielSMonteiro/DNA_Fountain/releases/tag/v1.0-defesa)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Institution](https://img.shields.io/badge/UnB-FCTE-darkred.svg)](https://fga.unb.br/)

Este repositório contém a implementação de referência, a suíte de testes determinísticos automatizados, os dados experimentais brutos e os scripts de reprodução do Trabalho de Conclusão de Curso (TCC) em Engenharia de Software na Universidade de Brasília (UnB / FCTE):

> **Título:** Implementação e Avaliação Experimental de um Pipeline DNA Fountain para Armazenamento de Dados em DNA Sintético em Regime de Blocos Finitos  
> **Autor:** Gabriel Monteiro  
> **Orientador:** Prof. Dr. Fábio Macedo Mendes  
> **Versão Oficial de Defesa:** Release [`v1.0-defesa`](https://github.com/GabrielSMonteiro/DNA_Fountain/releases/tag/v1.0-defesa)

---

## Sumário

1. [Visão Geral e Escopo](#1-visão-geral-e-escopo)
2. [Arquitetura de Software e Módulos](#2-arquitetura-de-software-e-módulos)
3. [Instalação e Dependências](#3-instalação-e-dependências)
4. [Reprodutibilidade Estrita e Controle de PRNGs](#4-reprodutibilidade-estrita-e-controle-de-prngs)
5. [Principais Resultados Experimentais](#5-principais-resultados-experimentais)
6. [Suíte de Testes Automatizados (18 testes)](#6-suíte-de-testes-automatizados-18-testes)
7. [Guia de Execução e Reprodução dos Resultados](#7-guia-de-execução-e-reprodução-dos-resultados)
8. [Estrutura do Repositório e Artefatos de Dados](#8-estrutura-do-repositório-e-artefatos-de-dados)
9. [Como Citar](#9-como-citar)
10. [Licença](#10-licença)

---

## 1. Visão Geral e Escopo

O projeto modela computacionalmente, em nível lógico, o pipeline completo da arquitetura **DNA Fountain** ([Erlich & Zielinski, *Science*, 2017](https://doi.org/10.1126/science.aaj2038)), com foco específico na caracterização do comportamento em regime de blocos finitos ($K \le 64$).

O pipeline cobre as dez etapas fundamentais de armazenamento:
1. **Ingestão e Segmentação**: Divisão da carga binária em $K$ blocos de tamanho fixo ($L = 4$ bytes) com preenchimento determinístico (*zero-padding*);
2. **Sorteio Estocástico via DRS**: Amostragem de grau $d$ e vizinhança a partir da Distribuição Robusta de Soliton (DRS);
3. **Combinação Linear Booleana**: Geração do payload da gota (*droplet*) via operações bitwise XOR sobre $\mathrm{GF}(2)$;
4. **Embaralhamento Aditivo (*Scrambling*)**: Desestruturação de padrões periódicos com máscara pseudoaleatória auto-inversível $(\text{seed} \oplus \texttt{0x5A5A})$;
5. **Ofuscação de Cabeçalho**: Aplicação da máscara XOR estática $\texttt{0x1B1B}$ ao identificador de 16 bits;
6. **Camada de Integridade**: Cálculo do CRC-8 (polinômio $x^8 + x^2 + x + 1$, $\texttt{0x07}$);
7. **Conversão Binário $\to$ DNA**: Mapeamento biunívoco a 2 bits por base (`00` $\to$ A, `01` $\to$ C, `10` $\to$ G, `11` $\to$ T);
8. **Filtragem Físico-Química (*Screening*)**: Descarte de sequências fora da janela de GC ($40\% \le \text{GC} \le 60\%$) ou com homopolímeros longos ($> 3$ bases idênticas consecutivas);
9. **Canal de Apagamento Estocástico**: Simulação de perdas de síntese e sequenciamento sob probabilidades $p_{\text{loss}} \in \{0{,}0;\; 0{,}10;\; 0{,}20\}$;
10. **Decodificação Iterativa (*Peeling*)**: Reconstrução dos blocos-fonte originais via propagação de equações unitárias em grafo bipartido com indexação invertida.

### Delimitação do Escopo
- **Canal de Comunicação**: O canal biológico é modelado como um Canal de Apagamento Binário (BEC) sem ruído de fase ou substituição pontual de bases (apagamento total ou recepção íntegra).
- **Sem Síntese Úmida**: O estudo é uma investigação estrita de Engenharia de Software e teoria da informação, sem manipulação de reagentes químicos ou sequenciamento físico em bancada (*wet lab*).

---

## 2. Arquitetura de Software e Módulos

O simulador é estruturado sob alta coesão e baixo acoplamento no pacote `dna_fountain/`:

| Módulo | Responsabilidade |
|---|---|
| `dna_fountain/config.py` | Configuração centralizada (`DNAFountainConfig`) com parâmetros da DRS, enquadramento e tolerâncias de *screening*. |
| `dna_fountain/preprocessing.py` | Segmentação da mensagem em $K$ blocos, cálculo de *zero-padding* determinístico e truncamento pós-decodificação. |
| `dna_fountain/soliton.py` | Implementação da Distribuição Robusta de Soliton: cálculo de $\rho(d)$, $\tau(d)$, $\mu(d)$, salvaguarda $m = \max(1, \lfloor K/R \rfloor)$ e sorteio determinístico via CDF. |
| `dna_fountain/lt_codec.py` | Codificador (`DNAFountainEncoder`) com *screening* com reposição e decodificador *peeling* (`DNAFountainDecoder`) com indexação reversa e rastreamento da fila de *ripple*. |
| `dna_fountain/gaussian.py` | Decodificador algébrico via Eliminação Gaussiana em $\mathrm{GF}(2)$ (`GaussianGF2Decoder`) para certificação de posto de matriz e auditoria de *stopping sets*. |
| `dna_fountain/screening.py` | Validador de restrições biológicas: janela de conteúdo GC ($40\% \le \text{GC} \le 60\%$) e limite estrito de homopolímeros ($H \le 3$). |
| `dna_fountain/channel.py` | Simulação do canal de apagamento estocástico (`ErasureChannel`) com gerador PRNG dedicado. |
| `dna_fountain/dna_conversion.py` | Mapeamento biunívoco bidirecional entre bytes e sequências de nucleotídeos (2 bits/base). |
| `dna_fountain/utils.py` | Rotinas de baixo nível: XOR vetorizado de blocos, CRC-8 padrão SMBus (`0x07`) e *scrambler* aditivo auto-inversível. |
| `dna_fountain/experiments.py` | Orquestração da execução em lote dos ensaios da grade fatorial e exportação dos resultados agregados. |

---

## 3. Instalação e Dependências

### Núcleo de Execução
O núcleo do simulador (`dna_fountain/`) foi projetado com **zero dependências externas**:
- Requer apenas **Python 3.12+ (CPython)** e utiliza exclusivamente a biblioteca padrão (`math`, `random`, `struct`, `collections`, `typing`, `csv`, `json`, `dataclasses`, `os`, `sys`).

### Ambiente de Testes e Visualização
Para rodar a suíte de testes unitários automatizados e regenerar os gráficos vetoriais da monografia:

```bash
# 1. Clonar o repositório
git clone https://github.com/GabrielSMonteiro/DNA_Fountain.git
cd DNA_Fountain

# 2. Criar e ativar ambiente virtual (recomendado)
python -m venv .venv
# No Windows:
.venv\Scripts\Activate.ps1
# No Linux/macOS:
source .venv/bin/activate

# 3. Instalar dependências de desenvolvimento
pip install -r requirements-dev.txt
```

O arquivo `requirements-dev.txt` inclui:
- `pytest>=8.0.0`: Suíte de testes unitários e de integração;
- `matplotlib>=3.8.0`: Geração das figuras vetoriais da monografia.

---

## 4. Reprodutibilidade Estrita e Controle de PRNGs

A reprodutibilidade exata em ambiente CPython 3.12 é assegurada pelo desacoplamento de três fluxos pseudoaleatórios independentes a partir de uma chave mestre composta:

### 4.1 Derivação da Semente Mestre (`seed_base`)
Em ensaios fatoriais (`dna_fountain/experiments.py`), cada ponto da grade calcula:
$$\text{seed}_{\text{base}} = \lfloor K \times 10^5 \rfloor + \lfloor c \times 10^4 \rfloor + \lfloor \delta \times 10^3 \rfloor + N_{\text{tx}} + \lfloor 100 \times p_{\text{loss}} \rfloor$$

A semente `seed_base` inicializa um gerador mestre que deriva três sementes de 31 bits dedicadas:
1. `rng_payload`: Gera os blocos sintéticos de carga útil (dados a armazenar);
2. `rng_seeds`: Sorteia uma semente inicial independente ($s_0 \in [0, 65535]$) para o codificador a cada repetição;
3. `rng_channel`: Controla os sorteios de apagamento no canal estocástico (parametrizado por $p_{\text{loss}}$).

### 4.2 Amostragem do Ensamble e Espaço de Sementes
- O sorteio de uma nova semente inicial $s_0$ a cada repetição assegura que as 1.000 repetições explorem grafos bipartidos distintos, medindo a taxa de sucesso do **ensamble do código LT** e não de uma trajetória pontual fixa de sementes.
- A busca interna por candidatos avança com incremento unitário ($s = s_0, s_0 + 1, \dots$).
- Sob o espaço de endereçamento de 16 bits ($2^{16} = 65.536$ sementes) e rejeição biológica observada de $42{,}6\%$, restam aproximadamente **37.600 sequências válidas**, capacidade mais que suficiente para os blocos investigados ($K \le 64$).

---

## 5. Principais Resultados Experimentais

A avaliação experimental compreendeu a grade fatorial completa de 216 combinações de parâmetros com **1.000 repetições independentes por ponto (216.000 simulações completas)**, consolidada em `resultados_1000.csv`:

1. **Validação Funcional sob Canal Ideal ($p_{\text{loss}} = 0{,}0$, $\epsilon = 1{,}00$)**:
   - Para $K = 32$, a configuração de referência histórica ($c = 0{,}1$, $\delta = 0{,}10$) atingiu **$92{,}0\%$ de sucesso** (IC 95\% Clopper--Pearson: $[90{,}1\%, 93{,}6\%]$).
   - Para $K = 64$, atingiu **$98{,}5\%$ de sucesso** (IC 95\%: $[97{,}5\%, 99{,}2\%]$).
   - Todos os blocos recuperados foram verificados byte a byte contra os dados originais, com zero erros lógicos.
2. **Comparativo Paramétrico: $c = 0{,}1$ vs. $c = 0{,}3$**:
   - No regime de blocos finitos ($K \le 64$), $c = 0{,}3$ apresentou desempenho igual ou superior em **8 das 9 células paramétricas** ($K \times \epsilon$).
   - Em orçamentos de redundância restritos ($\epsilon = 0{,}25$), $c = 0{,}3$ superou significativamente $c = 0{,}1$ tanto em $K = 16$ ($27{,}2\%$ vs. $16{,}9\%$, $p < 10^{-7}$) quanto em $K = 32$ ($23{,}9\%$ vs. $11{,}8\%$, $p < 10^{-10}$).
   - Sob $\epsilon = 1{,}00$, $c = 0{,}3$ atingiu $93{,}9\%$ vs. $92{,}0\%$ em $K = 32$ ($p = 0{,}12$, empate estatístico) e $89{,}3\%$ vs. $84{,}7\%$ em $K = 16$ ($p = 0{,}0027$).
   - Apenas em $K = 64$ sob $\epsilon = 1{,}00$ a configuração $c = 0{,}1$ foi superior ($98{,}5\%$ vs. $96{,}6\%$, $p = 0{,}0086$).
3. **Mecânica da Configuração Inflada ($c = 0{,}5$)**:
   - Em $K = 32$, a escala de dispersão $R \approx 16{,}32 > K/2$ induz o colapso do modo analítico para $m = \lfloor K/R \rfloor = 1$, concentrando cerca de $73\%$ da massa em grau 1. Isso reduz as gotas interligadas e limita o sucesso a $38{,}5\%$ sob $\epsilon = 1{,}00$ (e $3{,}7\%$ sob $\epsilon = 0{,}50$).
   - Em $K = 64$, a escala expande-se ($R \approx 25{,}85 < K/2$), promovendo o deslocamento modal para $m = 2$ ($\mu(2) \approx 0{,}752$), o que sustenta $69{,}2\%$ de sucesso sob $\epsilon = 1{,}00$.
4. **Mecânica do Screening e do CRC-8**:
   - Avaliação sobre 100.000 candidatos enquadrados revelou taxa de rejeição de $42{,}6\%$, contra $49{,}8\%$ do modelo analítico uniforme não restrito.
   - O checksum CRC-8 atua como restrição linear nos 6 bytes anteriores, contraindo a variância de GC e reduzindo violações de GC de $34{,}7\%$ para $24{,}8\%$ (violações de homopolímeros permanecem em $26{,}5\%$).
   - A máscara XOR `0x1B1B` é uma involução bijetora que mantém a uniformidade estrita do espaço de 16 bits e protege os prefixos de sementes baixas (rejeição de $39{,}2\%$ em $s \in [0, 2000]$).
   - O teste qui-quadrado de homogeneidade de Pearson entre aceitos e rejeitados não detectou diferença estatisticamente significativa entre as distribuições de graus ($\chi^2 = 23{,}66$, $\text{df} = 17$, $p = 0{,}1289$).
5. **Auditoria de Falhas em $K = 16$ com Eliminação Gaussiana em $\mathrm{GF}(2)$**:
   - Nas 153 falhas observadas em 1.000 repetições ($c = 0{,}1, \delta = 0{,}10, \epsilon = 1{,}00$), a validação cruzada algébrica provou que **$98{,}7\%$ (151 falhas)** foram causadas por *stopping sets* cíclicos com posto completo de matriz ($99{,}8\%$ de resolução pela Gaussiana), e apenas 2 falhas ($0{,}2\%$) decorreram de falta inicial de cobertura no grafo.

---

## 6. Suíte de Testes Automatizados (18 testes)

A corretude algorítmica e a consistência matemática são certificadas por **18 testes automatizados**:

```bash
python -m pytest tests/ -v
```

| # | Teste Implementado | Tipo | Propriedade Validada |
|---|---|---|---|
| 1 | `test_conversao_bytes_dna_reversibilidade` | Unitário | Bijeção estrita nos 256 valores de bytes e fluxos aleatórios de 1 MB. |
| 2 | `test_somatoria_drs` | Unitário | Normalização estrita da DRS ($\sum_{d=1}^K \mu(d) = 1{,}0$) em grade de $K, c, \delta$. |
| 3 | `test_somatoria_drs_salvaguarda_m_menor_1` | Unitário | Salvaguarda de limite $m_{\text{raw}} < 1 \implies m = 1$ para $c \ge 0{,}5$. |
| 4 | `test_screening_gc_limits` | Unitário | Rejeição correta de limites estritos: GC $< 40\%$ e GC $> 60\%$. |
| 5 | `test_screening_homopolimeros` | Unitário | Aceitação de homopolímeros $\le 3$ e rejeição estrita de corridas $\ge 4$ nt. |
| 6 | `test_segmentacao_padding_truncamento` | Unitário | Reversibilidade exata do preenchimento por zero-padding em mensagens não múltiplas de $L$. |
| 7 | `test_amostragem_deterministica_seed` | Unitário | Determinismo estrito: a mesma semente reproduz exatamente o mesmo grau e vizinhos. |
| 8 | `test_consistencia_vizinhancas_encoder_decoder` | Integração | Concordância exata dos grafos bipartidos regenerados entre codificador e decodificador. |
| 9 | `test_canal_apagamento_extremos` | Unitário | Preservação total sob $p_{\text{loss}} = 0{,}0$ e perda total de pacotes sob $p_{\text{loss}} = 1{,}0$. |
| 10 | `test_exportacao_metricas_csv` | Integração | Serialização e desserialização íntegras de métricas estruturadas em CSV. |
| 11 | `test_peeling_canal_ideal` | Integração | Pipeline completo ($K=32, \epsilon=1{,}00, p_{\text{loss}}=0{,}0$): recuperação de 100% dos blocos. |
| 12 | `test_grafo_canonico_gotas_unitarias` | Canônico | Grafo canônico de $K$ gotas de grau 1: resolução direta com $Q_{\max} = K$. |
| 13 | `test_grafo_canonico_cascata_triangular` | Canônico | Cascata triangular de peeling ($d_i = \{i-1, i\}$): resolução sequencial com $Q_{\max} = 1$. |
| 14 | `test_falha_stopping_set_sem_grau1` | Negativo | Grafo cíclico de grau $\ge 2$ sem gotas unitárias: interrupção imediata (0 recuperados). |
| 15 | `test_descarte_crc_corrompido` | Negativo | Mutação pontual de 1 base no DNA recebido: 100% das gotas corrompidas descartadas pelo CRC-8. |
| 16 | `test_esgotamento_espaco_seeds` | Negativo | Demanda de gotas superior ao espaço de sementes disponíveis: disparo de `RuntimeError`. |
| 17 | `test_validacao_cruzada_peeling_gaussiana` | Integração | Validação cruzada: sucesso do *peeling* implica posto completo $K$ na Eliminação Gaussiana em $\mathrm{GF}(2)$. |
| 18 | `test_aderencia_qui_quadrado_espaco_seeds` | Estatístico | Teste de aderência qui-quadrado sobre todo o espaço de 65.536 sementes ($\chi^2 = 30{,}01$, $\text{df} = 31$, $p = 0{,}5169 > 0{,}05$) e desvio máximo $< 0{,}004$. |

---

## 7. Guia de Execução e Reprodução dos Resultados

Todos os scripts de reprodução encontram-se na raiz do repositório:

### 7.1 Execução da Suíte de Testes
```bash
python -m pytest tests/
```
*(Saída esperada: 18 passed em ~0.6 segundos)*.

### 7.2 Regeneração das Figuras Vetoriais da Monografia (Figuras 1 a 5)
Regenera todos os gráficos em formato PDF vetorial e PNG em alta resolução (300 DPI) em `Relatório/figuras/`, com notação numérica em ponto decimal:
```bash
python gerar_graficos.py
```

### 7.3 Auditoria de Falhas em $K = 16$ com Eliminação Gaussiana (Seção 6.1)
Reproduz a análise das 153 falhas do baseline de $K = 16$, classificando-as entre falta de cobertura inicial e *stopping sets* cíclicos em `falhas_k16.json`:
```bash
python analisar_falhas_k16.py
```

### 7.4 Teste $\chi^2$ do Screening e Estresse de Baixa Entropia (Seção 6.2)
Executa o teste de homogeneidade de 100.000 oligonucleotídeos candidatos aceitos versus rejeitados e o teste de estresse do *scrambler* com carga nula (`0x00000000`), salvando em `screening_qui_quadrado.json`:
```bash
python salvar_qui_quadrado.py
```

### 7.5 Estatísticas Isoladas de Custo Computacional (Tabela 6.7)
Gera as estatísticas de tempo de decodificação sequencial isolado (média, desvio padrão, mediana e extremos) e operações XOR:
```bash
# Apenas exibe a tabela a partir das 1.000 amostras individuais arquivadas:
python salvar_estatisticas_custo.py

# Para executar uma nova bateria ativa de 1.000 ensaios por K:
python salvar_estatisticas_custo.py --novo-benchmark
```

### 7.6 Execução da Grade Fatorial Completa (216.000 ensaios)
```bash
# Execução paralelizada (recomendado, utiliza todos os núcleos da CPU):
python reexecutar_grade.py

# Ou execução sequencial tradicional:
python run_experiments.py
```

### 7.7 Compilação da Monografia em LaTeX (UnB / ABNT)
Para compilar o relatório completo em PDF a partir do código-fonte em `Relatório/`:

**No Windows (PowerShell):**
```powershell
powershell -ExecutionPolicy Bypass -File "Relatório\compilar.ps1"
```

**No Linux/macOS (Make/Bash):**
```bash
cd Relatório
pdflatex main
bibtex main
pdflatex main
pdflatex main
```


## 9. Como Citar

Caso utilize este simulador, a suíte de testes ou os dados experimentais em pesquisas acadêmicas, cite:

```bibtex
@monography{monteiro2026dna,
  author    = {Gabriel Monteiro},
  title     = {Implementation and Experimental Evaluation of a DNA Fountain Pipeline for Synthetic DNA Data Storage in Finite Block Regimes},
  school    = {Faculdade de Ci{\^e}ncias e Tecnologias em Engenharia, Universidade de Bras{\'\i}lia (UnB)},
  address   = {Bras{\'\i}lia, DF, Brazil},
  year      = {2026},
  type      = {Undergraduate Thesis (Bachelor in Software Engineering)},
  note      = {Release v1.0-defesa. Available from: \url{https://github.com/GabrielSMonteiro/DNA_Fountain}}
}
```

---

## 10. Licença

Este projeto é disponibilizado sob a **Licença MIT**. Consulte o arquivo [LICENSE](LICENSE) para obter detalhes sobre permissões e condições de uso.
