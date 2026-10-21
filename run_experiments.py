import csv
import json
import sys
import os

# Configurar stdout para UTF-8 no Windows
if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

from dna_fountain.experiments import rodar_experimentos
from dna_fountain.main import pipeline
from dna_fountain.config import DNAFountainConfig
from dna_fountain.preprocessing import MessageProcessor
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder

config = DNAFountainConfig()
processor = MessageProcessor(config)
encoder = DNAFountainEncoder(config)
decoder = DNAFountainDecoder(config)

print('--- 1. Demonstração ---')
config.num_droplets = 200  # Força overhead alto para compensar a ausência de Reed-Solomon em K pequeno
pipeline("Demonstracao inicial com uma mensagem longa o suficiente para evitar a altissima variancia de K pequenos.", config)

print('\n--- 2. Experimentos (10 repetições) ---')
rodar_experimentos('resultados_10.csv', 10)

print('\n--- 3. Experimentos (30 repetições) ---')
rodar_experimentos('resultados_30.csv', 30)

print('\n--- 4. Experimentos (100 repetições) ---')
rodar_experimentos('resultados_100.csv', 100)

print('\nExperimentos concluídos! Gerando gráficos em HTML...')

# Lendo CSV com módulo built-in
linhas = []
with open('resultados_100.csv', 'r', encoding='utf-8-sig') as f:
    leitor = csv.DictReader(f, delimiter=';')
    for linha in leitor:
        linhas.append(linha)

# Prepara dados para gráfico 1: Sucesso x Overhead, separados por K
dados_por_k = {'16': {}, '32': {}, '64': {}}
for K in ['16', '32', '64']:
    for c in ['0.1', '0.3', '0.5', '0.8']:
        label = f'c={c}'
        pontos = []
        for l in linhas:
            # Filtrar por delta=0.1 e p_perda=0.0 como recomendado
            if l['K'] == K and l['c'] == c and float(l['p_perda']) == 0.0 and float(l['delta']) == 0.1:
                pontos.append({
                    'x': float(l['overhead_enviado']),
                    'y': float(l['taxa_sucesso'])
                })
        pontos.sort(key=lambda p: p['x'])
        if pontos:
            dados_por_k[K][label] = pontos

# Preparando gráfico 2: Ripple por passo
MENSAGEM = "Demonstracao inicial com uma mensagem longa o suficiente para evitar a altissima variancia de K pequenos."
dados_originais, blocos, K = processor.preprocessar(MENSAGEM)

# Parâmetros para extrair um bom gráfico do ripple
config.soliton_c = 0.5
config.soliton_delta = 0.1
config.num_droplets = 70
soliton = RobustSoliton(config, K)

# Gera alguns oligos para análise do ripple
oligos, _, _ = encoder.codificar(blocos, soliton)

# Decodificação extrai automaticamente as métricas
dnas = [o["dna"] for o in oligos]
descobertos, metricas = decoder.decodificar(dnas, soliton)

# Pega o histórico real do tamanho do ripple que já vem calculado da nova arquitetura
dados_grafico2 = metricas["ripple_historico"]

# Gerando HTML
html_content = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Gráficos TCC</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        body {{ font-family: sans-serif; padding: 20px; }}
        .chart-container {{ width: 100%; max-width: 800px; margin: 0 auto 50px auto; }}
        h2 {{ text-align: center; }}
    </style>
</head>
<body>
    <h2>Taxa de Sucesso × Overhead Enviado (K=16, delta=0.1, p_perda=0.0)</h2>
    <div class="chart-container">
        <canvas id="chart16"></canvas>
    </div>
    
    <h2>Taxa de Sucesso × Overhead Enviado (K=32, delta=0.1, p_perda=0.0)</h2>
    <div class="chart-container">
        <canvas id="chart32"></canvas>
    </div>
    
    <h2>Taxa de Sucesso × Overhead Enviado (K=64, delta=0.1, p_perda=0.0)</h2>
    <div class="chart-container">
        <canvas id="chart64"></canvas>
    </div>

    <h2>Evolução do ripple (Demonstração K={K})</h2>
    <div class="chart-container">
        <canvas id="chartRipple"></canvas>
    </div>
    <script>
        const colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown'];
        
        function criarGraficoK(idCanvas, dadosK) {{
            const ctx = document.getElementById(idCanvas).getContext('2d');
            const datasets = [];
            let cIdx = 0;
            for (const [label, data] of Object.entries(dadosK)) {{
                datasets.push({{
                    label: label,
                    data: data,
                    borderColor: colors[cIdx % colors.length],
                    fill: false,
                    showLine: true,
                    tension: 0.1
                }});
                cIdx++;
            }}
            new Chart(ctx, {{
                type: 'scatter',
                data: {{ datasets: datasets }},
                options: {{
                    scales: {{
                        x: {{ title: {{ display: true, text: 'Overhead Enviado' }} }},
                        y: {{ title: {{ display: true, text: 'Taxa de Sucesso' }}, min: 0, max: 1 }}
                    }}
                }}
            }});
        }}

        const dadosK = {json.dumps(dados_por_k)};
        criarGraficoK('chart16', dadosK['16']);
        criarGraficoK('chart32', dadosK['32']);
        criarGraficoK('chart64', dadosK['64']);

        const ctxRipple = document.getElementById('chartRipple').getContext('2d');
        const dadosRipple = {json.dumps(dados_grafico2)};
        new Chart(ctxRipple, {{
            type: 'line',
            data: {{
                labels: dadosRipple.map((_, i) => i + 1),
                datasets: [{{
                    label: 'Tamanho do Ripple',
                    data: dadosRipple,
                    borderColor: 'rgb(75, 192, 192)',
                    tension: 0.1
                }}]
            }},
            options: {{
                scales: {{
                    x: {{ title: {{ display: true, text: 'Passo de processamento do ripple' }} }},
                    y: {{ title: {{ display: true, text: 'Tamanho do ripple' }} }}
                }}
            }}
        }});
    </script>
</body>
</html>
"""

with open('graficos.html', 'w', encoding='utf-8') as f:
    f.write(html_content)

print('Gráficos gerados com sucesso no arquivo graficos.html')
