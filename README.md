# Simulador do Pipeline DNA Fountain

Este repositório contém o código-fonte, a suíte de testes determinísticos e os dados experimentais do Trabalho de Conclusão de Curso (TCC) em Engenharia de Software na Universidade de Brasília (UnB / FCTE), intitulado:

> **Implementação e Avaliação Experimental de um Pipeline DNA Fountain para Armazenamento de Dados em DNA Sintético em Regime de Blocos Finitos**  
> **Autor:** Gabriel Monteiro  
> **Versão de Referência do Relatório:** Commit `bb2f0ab` (consolidado em `adeb76c`)

---

## 1. Visão Geral da Arquitetura

O simulador modela computacionalmente em nível lógico o pipeline de armazenamento do DNA Fountain (Erlich & Zielinski, *Science*, 2017), dividido em módulos especializados:

- `dna_fountain/config.py`: Centralização de hiperparâmetros e configurações do sistema.
- `dna_fountain/preprocessing.py`: Segmentação em blocos de $L$ bytes, zero-padding determinístico e reconstrução.
- `dna_fountain/soliton.py`: Distribuição Robusta de Soliton (DRS), cálculo analítico de $\rho(d)$, $\tau(d)$, $\mu(d)$ e amostragem determinística via CDF.
- `dna_fountain/lt_codec.py`: Codificador (`DNAFountainEncoder`) com screening com reposição e decodificador peeling (`DNAFountainDecoder`) com indexação reversa em grafo bipartido e monitoramento do *ripple*.
- `dna_fountain/screening.py`: Validação físico-química de conteúdo GC ($40\% \le \text{GC} \le 60\%$) e limite de homopolímeros ($\le 3$ nt).
- `dna_fountain/channel.py`: Canal de apagamento estocástico de pacotes ($p_{\text{loss}} \in \{0{,}0;\; 0{,}10;\; 0{,}20\}$) com PRNG independente.
- `dna_fountain/dna_conversion.py`: Mapeamento biunívoco de 2 bits por base (`00` $\to$ A, `01` $\to$ C, `10` $\to$ G, `11` $\to$ T).
- `dna_fountain/utils.py`: Operações XOR entre blocos, CRC-8 (polinômio `0x07`) e embaralhamento (*scrambling*).
- `dna_fountain/experiments.py`: Automação da execução de ensaios fatoriais em lote.

---

## 2. Separação de Ambientes e Dependências

Para assegurar longevidade e portabilidade acadêmica, o projeto estabelece uma separação rigorosa entre núcleo de execução e ferramentas de desenvolvimento:

1. **Núcleo do Simulador (`dna_fountain/` e `run_experiments.py`)**:
   - Depende **exclusivamente da biblioteca padrão do Python** (CPython 3.12.7).
   - Não requer a instalação de nenhum pacote externo (`pip`). Módulos utilizados: `math`, `random`, `struct`, `collections`, `typing`, `csv`, `json`, `time`, `dataclasses`, `os`, `sys`.
2. **Ambiente de Desenvolvimento, Testes e Visualização**:
   - As dependências de instrumentação estão descritas em `requirements-dev.txt`:
     - `pytest>=8.0.0`: Execução da suíte de testes unitários e de integração.
     - `matplotlib>=3.8.0`: Regeneração dos gráficos vetoriais (.pdf) e rasterizados (.png).

Para instalar as dependências de desenvolvimento:
```bash
pip install -r requirements-dev.txt
```

---

## 3. Especificação do Algoritmo de Sementes e Reprodutibilidade

A reprodutibilidade estrita no ambiente CPython 3.12.7 é governada por três políticas determinísticas:

### 3.1 Derivação da Semente Mestre (`seed_base`)
Em ensaios fatoriais em lote (`dna_fountain/experiments.py`), cada ponto paramétrico deriva uma semente mestre única através da fórmula:
$$\text{seed\_base} = \lfloor K \times 10^5 \rfloor + \lfloor c \times 10^4 \rfloor + \lfloor \delta \times 10^3 \rfloor + N_{\text{tx}}$$
Essa transformação mapeia os valores fracionários de $c \in \{0.1, 0.3, 0.5, 0.8\}$ e $\delta \in \{0.01, 0.10\}$ para inteiros únicos, prevenindo colisões entre combinações da grade fatorial. A semente `seed_base` inicializa o PRNG responsável pela carga útil da mensagem e pelos eventos estocásticos de apagamento no canal.

### 3.2 Sementes de Droplets e Espaço Amostral
- O codificador inicia a amostragem com `seed = 0` e incrementa `seed += 1` sequencialmente.
- Cada candidato avalia `RobustSoliton.sortear_grau_e_indices(seed)` com uma instância isolada de `random.Random(seed)`, garantindo que o grau e os índices dependam exclusivamente da semente do droplet.
- Candidatos reprovados no screening biológico são descartados com reposição determinística (avanço de `seed`).
- O espaço de endereçamento é limitado a 16 bits ($2^{16} = 65.536$ valores, intervalo $[0, 65.535]$). Caso a busca atinja $65.536$, o codificador encerra a execução com `RuntimeError("Espaço de seeds esgotado.")`.

---

## 4. Guia de Execução e Reprodução de Resultados

Todos os comandos devem ser executados a partir da raiz do repositório:

### 4.1 Suíte de Testes Automatizados (10 testes)
Executa a validação determinística de corretude algorítmica:
```bash
python -m pytest tests/
```
*(Saída esperada: 10 passed em ~0.05 segundos)*.

### 4.2 Regeneração de Tabelas e Estatísticas de Custo (Tabela 15)
- **Leitura do conjunto arquivado** (regenera estatísticas e gráficos idênticos aos publicados no TCC sem variações de escalonamento do SO):
  ```bash
  python salvar_estatisticas_custo.py
  ```
- **Execução de novo benchmark ativo** (executa 100 novas medições por configuração e atualiza o arquivo JSON):
  ```bash
  python salvar_estatisticas_custo.py --novo-benchmark
  ```

### 4.3 Auditoria de Falhas em $K=16$ (Seção 6.1)
Gera e analisa as 10 falhas observadas no cenário baseline de $K=16$:
```bash
python analisar_falhas_k16.py
```
*(Classifica as falhas entre isolamento inicial e stopping sets em `falhas_k16.json`)*.

### 4.4 Histograma e Teste $\chi^2$ do Screening Biológico (Seção 6.2)
Gera as contagens observadas e esperadas das 18 classes sob 2.000 candidatos:
```bash
python salvar_qui_quadrado.py
```
*(Gera `screening_qui_quadrado.json` com $\chi^2 = 25{,}66$, df = 17, $p = 0{,}081$)*.

### 4.5 Regeneração das Figuras e Gráficos da Monografia
```bash
python gerar_graficos.py
python gerar_diagrama_arquitetura.py
```

### 4.6 Execução Completa da Bateria Fatorial (21.600 ensaios)
Para reexecutar a grade experimental completa (3 matrizes de 216 linhas):
```bash
python run_experiments.py
```

---

## 5. Estrutura dos Dados e Disponibilidade

- **Matrizes Agregadas da Grade Fatorial**:
  - `resultados_10.csv`, `resultados_30.csv`, `resultados_100.csv`: Contêm 216 linhas cada, com médias, desvios e taxas observadas sobre 10, 30 e 100 repetições.
- **Registros Individuais de Auditoria**:
  - `tempos_brutos_decodificacao.json`: 100 medições individuais de tempo de decodificação e contagem de operações XOR para $K \in \{16, 32, 64\}$.
  - `falhas_k16.json`: Registros detalhados das 10 repetições com falha em $K=16$.
  - `screening_qui_quadrado.json`: Contagens observadas $O_j$ e esperadas $E_j$ das 18 classes do teste $\chi^2$ para 2.000 oligonucleotídeos candidatos.

