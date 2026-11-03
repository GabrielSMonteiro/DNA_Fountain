import csv
import os
import matplotlib.pyplot as plt
import numpy as np

# Configurar estilo visual academico
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.titlesize': 12,
    'lines.linewidth': 1.8,
    'lines.markersize': 6,
    'figure.autolayout': True
})

output_dir = os.path.join('Relatório', 'figuras')
os.makedirs(output_dir, exist_ok=True)

# 1. Carregar dados do CSV
rows = []
with open('resultados_100.csv', 'r', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f, delimiter=','):
        rows.append(r)

print(f'Total de registros carregados: {len(rows)}')

# -------------------------------------------------------------
# FIGURA 1: Dinâmica do Ripple (Passo a Passo)
# -------------------------------------------------------------
from dna_fountain.config import DNAFountainConfig
from dna_fountain.preprocessing import MessageProcessor
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder

msg = "Demonstracao experimental de decodificacao peeling para DNA Fountain simulando dinamica de ripple."
cfg = DNAFountainConfig()
proc = MessageProcessor(cfg)
_, blocos, K = proc.preprocessar(msg)

# Caso 1: c=0.1, delta=0.1 (estavel)
cfg.soliton_c = 0.1
cfg.soliton_delta = 0.1
cfg.num_droplets = int(K * 2.0)
sol1 = RobustSoliton(cfg, K)
enc1 = DNAFountainEncoder(cfg)
dec1 = DNAFountainDecoder(cfg)
oligos1, _, _ = enc1.codificar(blocos, sol1)
desc1, met1 = dec1.decodificar([o["dna"] for o in oligos1], sol1)
hist_ripple1 = met1["ripple_historico"]

# Caso 2: c=0.5, delta=0.1 (degenerado / inflado)
cfg.soliton_c = 0.5
cfg.soliton_delta = 0.1
cfg.num_droplets = int(K * 2.0)
sol2 = RobustSoliton(cfg, K)
enc2 = DNAFountainEncoder(cfg)
dec2 = DNAFountainDecoder(cfg)
oligos2, _, _ = enc2.codificar(blocos, sol2)
desc2, met2 = dec2.decodificar([o["dna"] for o in oligos2], sol2)
hist_ripple2 = met2["ripple_historico"]

fig, ax = plt.subplots(figsize=(6.5, 3.8))
ax.plot(range(1, len(hist_ripple1) + 1), hist_ripple1, marker='o', color='#1f77b4', label=f'Parametrização Ótima ($c=0.1$, $\\delta=0.1$) - Sucesso: {len(desc1)}/{K}')
if hist_ripple2:
    ax.plot(range(1, len(hist_ripple2) + 1), hist_ripple2, marker='s', color='#d62728', linestyle='--', label=f'Parametrização Degenerada ($c=0.5$, $\\delta=0.1$) - Recuperados: {len(desc2)}/{K}')

ax.set_xlabel('Passo de Processamento do Peeling')
ax.set_ylabel('Tamanho da Fila do Ripple')
ax.set_title('Evolução do Tamanho do Ripple Durante a Decodificação ($K=32$, $\\epsilon=1.0$)')
ax.legend(loc='upper right', frameon=True)
ax.grid(True, linestyle=':', alpha=0.6)
fig.savefig(os.path.join(output_dir, 'grafico_ripple.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_ripple.png'), dpi=300)
plt.close(fig)
print('Figura 1 (grafico_ripple) gerada.')

# -------------------------------------------------------------
# FIGURA 2: Overhead vs Taxa de Sucesso (Painéis para K=16, 32, 64)
# -------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8), sharey=True)
ks = ['16', '32', '64']
c_vals = ['0.1', '0.3', '0.5', '0.8']
cores = {'0.1': '#1f77b4', '0.3': '#2ca02c', '0.5': '#ff7f0e', '0.8': '#d62728'}
marcadores = {'0.1': 'o', '0.3': 's', '0.5': '^', '0.8': 'D'}

for idx, k in enumerate(ks):
    ax = axes[idx]
    for c in c_vals:
        pontos = []
        for r in rows:
            if r['K'] == k and r['c'] == c and r['delta'] == '0.1' and float(r['p_perda']) == 0.0:
                pontos.append((float(r['overhead_enviado']) * 100, float(r['taxa_sucesso']) * 100))
        pontos.sort(key=lambda x: x[0])
        if pontos:
            xs, ys = zip(*pontos)
            ax.plot(xs, ys, marker=marcadores[c], color=cores[c], label=f'$c={c}$')
    
    ax.set_title(f'Tamanho de Bloco $K = {k}$')
    ax.set_xlabel('Overhead Enviado $\\epsilon$ (%)')
    ax.set_xticks([25, 50, 100])
    ax.grid(True, linestyle=':', alpha=0.6)
    if idx == 0:
        ax.set_ylabel('Taxa de Sucesso (%)')
    if idx == 1:
        ax.legend(loc='center left', frameon=True)

fig.suptitle('Taxa de Sucesso na Decodificação em Função do Overhead ($\\delta = 0.10$, $p_{\\text{loss}} = 0.0$)', fontsize=12)
fig.savefig(os.path.join(output_dir, 'grafico_overhead_sucesso.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_overhead_sucesso.png'), dpi=300)
plt.close(fig)
print('Figura 2 (grafico_overhead_sucesso) gerada.')

# -------------------------------------------------------------
# FIGURA 3: Impacto de c e delta no Grau Médio e Sucesso (K=32)
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

for d, cor, est in [('0.1', '#1f77b4', '-'), ('0.01', '#ff7f0e', '--')]:
    pts_grau = []
    pts_succ = []
    for c in ['0.1', '0.3', '0.5', '0.8']:
        for r in rows:
            if r['K'] == '32' and r['c'] == c and r['delta'] == d and r['overhead_enviado'] == '0.5' and float(r['p_perda']) == 0.0:
                pts_grau.append((float(c), float(r['grau_medio'])))
                pts_succ.append((float(c), float(r['taxa_sucesso']) * 100))
    pts_grau.sort()
    pts_succ.sort()
    
    c_x, g_y = zip(*pts_grau)
    _, s_y = zip(*pts_succ)
    
    ax1.plot(c_x, g_y, marker='o', color=cor, linestyle=est, label=f'$\\delta = {d}$')
    ax2.plot(c_x, s_y, marker='s', color=cor, linestyle=est, label=f'$\\delta = {d}$')

ax1.set_title('Grau Médio Amostrado ($\\bar{d}$)')
ax1.set_xlabel('Parâmetro $c$')
ax1.set_ylabel('Grau Médio')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(frameon=True)

ax2.set_title('Taxa de Sucesso (%)')
ax2.set_xlabel('Parâmetro $c$')
ax2.set_ylabel('Sucesso (%)')
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(frameon=True)

fig.suptitle('Impacto dos Parâmetros $c$ e $\\delta$ na DRS e no Sucesso ($K=32$, $\\epsilon=0.50$, $p_{\\text{loss}}=0.0$)', fontsize=12)
fig.savefig(os.path.join(output_dir, 'grafico_impacto_c_delta.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_impacto_c_delta.png'), dpi=300)
plt.close(fig)
print('Figura 3 (grafico_impacto_c_delta) gerada.')

# -------------------------------------------------------------
# FIGURA 4: Resiliência sob Perdas no Canal (K=32 e K=64)
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8), sharey=True)

oh_cores = {'0.25': '#1f77b4', '0.5': '#ff7f0e', '1.0': '#2ca02c'}
oh_labels = {'0.25': '$\\epsilon = 0.25$ (25%)', '0.5': '$\\epsilon = 0.50$ (50%)', '1.0': '$\\epsilon = 1.00$ (100%)'}
oh_marks = {'0.25': '^', '0.5': 's', '1.0': 'o'}

for k, ax, titulo in [('32', ax1, '$K=32$'), ('64', ax2, '$K=64$')]:
    for oh in ['0.25', '0.5', '1.0']:
        pts = []
        for r in rows:
            if r['K'] == k and r['c'] == '0.1' and r['delta'] == '0.1' and r['overhead_enviado'] == oh:
                pts.append((float(r['p_perda']) * 100, float(r['taxa_sucesso']) * 100))
        pts.sort()
        if pts:
            px, py = zip(*pts)
            ax.plot(px, py, marker=oh_marks[oh], color=oh_cores[oh], label=oh_labels[oh])
    
    ax.set_title(f'Blocos Fonte {titulo}')
    ax.set_xlabel('Probabilidade de Perda $p_{\\text{loss}}$ (%)')
    ax.set_xticks([0, 10, 20])
    ax.grid(True, linestyle=':', alpha=0.6)

ax1.set_ylabel('Taxa de Sucesso (%)')
ax1.legend(loc='lower left', frameon=True)

fig.suptitle('Resiliência do Sistema sob Perdas Estocásticas no Canal ($c=0.1$, $\\delta=0.10$)', fontsize=12)
fig.savefig(os.path.join(output_dir, 'grafico_resiliencia.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_resiliencia.png'), dpi=300)
plt.close(fig)
print('Figura 4 (grafico_resiliencia) gerada.')

# -------------------------------------------------------------
# FIGURA 5: Custo Computacional (Tempo e XORs)
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

ks_num = [16, 32, 64]
tempos_med = [1.37, 2.34, 3.84]
tempos_min = [0.86, 0.96, 1.83]
tempos_max = [2.27, 4.63, 8.29]
erros_inf = [m - mi for m, mi in zip(tempos_med, tempos_min)]
erros_sup = [ma - m for m, ma in zip(tempos_med, tempos_max)]

xors_totais = [107.6, 270.8, 688.4]
xors_por_bloco = [107.6 / 16, 270.8 / 32, 688.4 / 64]

# Painel 1: Tempo
ax1.errorbar(ks_num, tempos_med, yerr=[erros_inf, erros_sup], fmt='-o', color='#1f77b4', ecolor='#aec7e8', elinewidth=2, capsize=5, label='Tempo Médio (com mín./máx.)')
ax1.set_title('Tempo de Decodificação (ms)')
ax1.set_xlabel('Número de Blocos Fonte ($K$)')
ax1.set_ylabel('Tempo (ms)')
ax1.set_xticks(ks_num)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(frameon=True)

# Painel 2: XORs
ax2.plot(ks_num, xors_totais, '-s', color='#ff7f0e', label='Total de XORs')
ax2_twin = ax2.twinx()
ax2_twin.plot(ks_num, xors_por_bloco, '--^', color='#2ca02c', label='XORs por Bloco ($/K$)')
ax2.set_title('Operações XOR na Decodificação')
ax2.set_xlabel('Número de Blocos Fonte ($K$)')
ax2.set_ylabel('Total de XORs', color='#ff7f0e')
ax2_twin.set_ylabel('XORs por Bloco Fonte', color='#2ca02c')
ax2.set_xticks(ks_num)
ax2.grid(True, linestyle=':', alpha=0.6)

# Combinar legendas
lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True)

fig.suptitle('Custo Computacional e Escalabilidade do Peeling Decoder ($c=0.1$, $\\delta=0.10$, $\\epsilon=1.00$)', fontsize=12)
fig.savefig(os.path.join(output_dir, 'grafico_custo_computacional.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_custo_computacional.png'), dpi=300)
plt.close(fig)
print('Figura 5 (grafico_custo_computacional) gerada.')

print('Todos os 5 graficos foram gerados com sucesso!')

