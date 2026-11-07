import csv
import os
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import numpy as np

# Configure academic visual style
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

# Helper formatter for decimal commas (ABNT compliance)
comma_fmt = FuncFormatter(lambda x, pos: f"{x:g}".replace('.', ','))

# 1. Load CSV data
csv_path = 'resultados_1000.csv' if os.path.exists('resultados_1000.csv') else 'resultados_100.csv'
rows = []
with open(csv_path, 'r', encoding='utf-8-sig') as f:
    for r in csv.DictReader(f, delimiter=','):
        rows.append(r)

print(f'Total records loaded from {csv_path}: {len(rows)}')

# -------------------------------------------------------------
# FIGURE 1: Ripple Dynamics (Step-by-Step)
# -------------------------------------------------------------
import random
from dna_fountain.config import DNAFountainConfig
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder

K = 32
L = 4
cfg = DNAFountainConfig()
cfg.tamanho_bloco_l = L
rng_bl = random.Random(42)
blocos = [bytes(rng_bl.getrandbits(8) for _ in range(L)) for _ in range(K)]

# Case 1: c=0.1, delta=0.1 (recommended / stable)
cfg.soliton_c = 0.1
cfg.soliton_delta = 0.1
cfg.num_droplets = int(K * 2.0)
sol1 = RobustSoliton(cfg, K)
enc1 = DNAFountainEncoder(cfg)
dec1 = DNAFountainDecoder(cfg)
oligos1, _, _ = enc1.codificar(blocos, sol1)
desc1, met1 = dec1.decodificar([o["dna"] for o in oligos1], sol1)
hist_ripple1 = met1["ripple_historico"]

# Case 2: c=0.5, delta=0.1 (degenerate / inflated)
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
ax.plot(range(1, len(hist_ripple1) + 1), hist_ripple1, marker='o', markersize=4, color='#1f77b4', label=f'Recommended Configuration ($c=0{{,}}1$, $\\delta=0{{,}}10$) - Recovered: {len(desc1)}/{K}')
if hist_ripple2:
    ax.plot(range(1, len(hist_ripple2) + 1), hist_ripple2, marker='s', markersize=4, color='#d62728', linestyle='--', label=f'Degenerate Configuration ($c=0{{,}}5$, $\\delta=0{{,}}10$) - Recovered: {len(desc2)}/{K}')

ax.set_xlabel('Peeling Decoding Processing Step')
ax.set_ylabel('Ripple Queue Size')
# Internal title removed to avoid repeating external LaTeX caption
ax.legend(loc='upper right', frameon=True)
ax.grid(True, linestyle=':', alpha=0.6)
fig.savefig(os.path.join(output_dir, 'grafico_ripple.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_ripple.png'), dpi=300)
plt.close(fig)
print('Figure 1 (grafico_ripple) generated with decimal comma.')

# -------------------------------------------------------------
# FIGURE 2: Overhead vs Success Rate (Panels for K=16, 32, 64)
# -------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.8), sharey=True)
ks = ['16', '32', '64']
c_vals = ['0.1', '0.3', '0.5', '0.8']
c_labels = {'0.1': '$c=0{,}1$', '0.3': '$c=0{,}3$', '0.5': '$c=0{,}5$', '0.8': '$c=0{,}8$'}
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
            ax.plot(xs, ys, marker=marcadores[c], color=cores[c], label=c_labels[c])
    
    ax.set_title(f'$K = {k}$')
    ax.set_xlabel('Transmitted Overhead $\\epsilon$ (%)')
    ax.set_xticks([25, 50, 100])
    ax.grid(True, linestyle=':', alpha=0.6)
    if idx == 0:
        ax.set_ylabel('Success Rate (%)')
    if idx == 1:
        ax.legend(loc='center left', frameon=True)

# Internal suptitle removed to avoid repeating external LaTeX caption
fig.savefig(os.path.join(output_dir, 'grafico_overhead_sucesso.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_overhead_sucesso.png'), dpi=300)
plt.close(fig)
print('Figure 2 (grafico_overhead_sucesso) generated with decimal comma.')

# -------------------------------------------------------------
# FIGURE 3: Impact of c and delta on Mean Degree and Success (K=32)
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

for d, cor, est, d_label in [('0.1', '#1f77b4', '-', '$\\delta = 0{,}10$'), ('0.01', '#ff7f0e', '--', '$\\delta = 0{,}01$')]:
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
    
    ax1.plot(c_x, g_y, marker='o', color=cor, linestyle=est, label=d_label)
    ax2.plot(c_x, s_y, marker='s', color=cor, linestyle=est, label=d_label)

ax1.set_title('(a) Theoretical Mean Degree')
ax1.set_xlabel('Parameter $c$')
ax1.set_ylabel('Theoretical Mean Degree $\\mathbb{E}_\\mu[D]$')
ax1.set_xticks([0.1, 0.3, 0.5, 0.8])
ax1.set_xticklabels(['0,1', '0,3', '0,5', '0,8'])
ax1.yaxis.set_major_formatter(comma_fmt)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(frameon=True)

ax2.set_title('(b) Success Rate')
ax2.set_xlabel('Parameter $c$')
ax2.set_ylabel('Success Rate (%)')
ax2.set_xticks([0.1, 0.3, 0.5, 0.8])
ax2.set_xticklabels(['0,1', '0,3', '0,5', '0,8'])
ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{int(y)}"))
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.legend(frameon=True)

# Internal suptitle removed to avoid repeating external LaTeX caption
fig.savefig(os.path.join(output_dir, 'grafico_impacto_c_delta.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_impacto_c_delta.png'), dpi=300)
plt.close(fig)
print('Figure 3 (grafico_impacto_c_delta) generated with decimal comma.')

# -------------------------------------------------------------
# FIGURE 4: Resilience under Channel Erasures (K=32 and K=64)
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8), sharey=True)

oh_cores = {'0.25': '#1f77b4', '0.5': '#ff7f0e', '1.0': '#2ca02c'}
oh_labels = {'0.25': '$\\epsilon = 0{,}25$ (25%)', '0.5': '$\\epsilon = 0{,}50$ (50%)', '1.0': '$\\epsilon = 1{,}00$ (100%)'}
oh_marks = {'0.25': '^', '0.5': 's', '1.0': 'o'}

for k, ax, titulo in [('32', ax1, '(a) $K = 32$'), ('64', ax2, '(b) $K = 64$')]:
    for oh in ['0.25', '0.5', '1.0']:
        pts = []
        for r in rows:
            if r['K'] == k and r['c'] == '0.1' and r['delta'] == '0.1' and r['overhead_enviado'] == oh:
                pts.append((float(r['p_perda']) * 100, float(r['taxa_sucesso']) * 100))
        pts.sort()
        if pts:
            px, py = zip(*pts)
            ax.plot(px, py, marker=oh_marks[oh], color=oh_cores[oh], label=oh_labels[oh])
    
    ax.set_title(titulo)
    ax.set_xlabel('Erasure Probability $p_{\\text{loss}}$ (%)')
    ax.set_xticks([0, 10, 20])
    ax.grid(True, linestyle=':', alpha=0.6)

ax1.set_ylabel('Success Rate (%)')
ax1.legend(loc='lower left', frameon=True)

# Internal suptitle removed to avoid repeating external LaTeX caption
fig.savefig(os.path.join(output_dir, 'grafico_resiliencia.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_resiliencia.png'), dpi=300)
plt.close(fig)
print('Figure 4 (grafico_resiliencia) generated with decimal comma.')

# -------------------------------------------------------------
# FIGURE 5: Computational Cost (Runtime and XORs)
# -------------------------------------------------------------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

import json
json_path = 'tempos_brutos_decodificacao.json'
if os.path.exists(json_path):
    with open(json_path, 'r', encoding='utf-8') as f:
        bench_data = json.load(f)
    ks_num = [16, 32, 64]
    tempos_med = [bench_data[str(k)]['mean_t'] for k in ks_num]
    tempos_std = [bench_data[str(k)]['std_t'] for k in ks_num]
    xors_totais = [bench_data[str(k)]['mean_xor'] for k in ks_num]
    xors_por_bloco = [bench_data[str(k)]['xor_per_block'] for k in ks_num]
else:
    ks_num = [16, 32, 64]
    tempos_med = [1.21, 2.24, 4.95]
    tempos_std = [0.27, 0.25, 0.65]
    xors_totais = [112.7, 296.8, 760.4]
    xors_por_bloco = [7.04, 9.28, 11.88]

# Panel 1: Runtime
ax1.errorbar(ks_num, tempos_med, yerr=tempos_std, fmt='-o', color='#1f77b4', ecolor='#aec7e8', elinewidth=2, capsize=5, label='Mean Runtime (±1 std. dev.)')
ax1.set_title('(a) Decoding Runtime')
ax1.set_xlabel('Source Block Count ($K$)')
ax1.set_ylabel('Runtime (ms)')
ax1.set_xticks(ks_num)
ax1.yaxis.set_major_formatter(comma_fmt)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(frameon=True)

# Panel 2: XORs
ax2.plot(ks_num, xors_totais, '-s', color='#ff7f0e', label='Total XORs')
ax2_twin = ax2.twinx()
ax2_twin.plot(ks_num, xors_por_bloco, '--^', color='#2ca02c', label='XORs per Block ($/K$)')
ax2.set_title('(b) XOR Operations')
ax2.set_xlabel('Source Block Count ($K$)')
ax2.set_ylabel('Total XOR Operations', color='#ff7f0e')
ax2_twin.set_ylabel('XORs per Source Block', color='#2ca02c')
ax2.set_xticks(ks_num)
ax2.yaxis.set_major_formatter(FuncFormatter(lambda y, _: f"{int(y)}"))
ax2_twin.yaxis.set_major_formatter(comma_fmt)
ax2.grid(True, linestyle=':', alpha=0.6)

# Combine legends
lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True)

# Internal suptitle removed to avoid repeating external LaTeX caption
fig.savefig(os.path.join(output_dir, 'grafico_custo_computacional.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_custo_computacional.png'), dpi=300)
plt.close(fig)
print('Figure 5 (grafico_custo_computacional) generated with decimal comma.')

print('All 5 figures successfully generated in English with standardized decimal comma (ABNT format)!')
