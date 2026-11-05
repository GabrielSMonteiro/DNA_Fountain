import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

fig, ax = plt.subplots(figsize=(10, 6.2), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Cores
c_box_bg = '#F8FAFC'
c_box_border = '#334155'
c_header_bg = '#1E293B'
c_header_txt = '#FFFFFF'
c_proc_bg = '#EFF6FF'
c_proc_border = '#2563EB'
c_dna_bg = '#F0FDF4'
c_dna_border = '#16A34A'
c_chan_bg = '#FEF2F2'
c_chan_border = '#DC2626'
c_dec_bg = '#FAF5FF'
c_dec_border = '#9333EA'
c_arrow = '#475569'

def draw_box(ax, x, y, w, h, title, subtitle, items, bg_color, border_color, title_color='#0F172A'):
    # Caixa principal com cantos arredondados
    rect = patches.FancyBboxPatch((x, y), w, h,
                                  boxstyle="round,pad=0.5,rounding_size=1.5",
                                  facecolor=bg_color, edgecolor=border_color, linewidth=1.5)
    ax.add_patch(rect)
    
    # Titulo
    ax.text(x + w/2, y + h - 2.8, title, ha='center', va='center', fontsize=9.5, fontweight='bold', color=title_color)
    if subtitle:
        ax.text(x + w/2, y + h - 5.5, subtitle, ha='center', va='center', fontsize=8, fontstyle='italic', color='#64748B')
        y_start = y + h - 8.2
    else:
        y_start = y + h - 5.8
        
    # Itens
    for i, it in enumerate(items):
        ax.text(x + 1.8, y_start - i * 3.0, it, ha='left', va='center', fontsize=7.8, color='#334155')

def draw_arrow(ax, x1, y1, x2, y2, label=''):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="->", color=c_arrow, lw=1.6, mutation_scale=12))
    if label:
        ax.text((x1 + x2)/2, (y1 + y2)/2 + 1.2, label, ha='center', va='bottom', fontsize=7.5, fontweight='semibold', color='#1E293B')

# Titulo Geral
ax.text(50, 97.5, 'Arquitetura e Fluxo de Execução do Simulador DNA Fountain', 
        ha='center', va='center', fontsize=12, fontweight='bold', color='#0F172A')

# 1. Entrada e Preprocessamento
draw_box(ax, 3, 68, 26, 24, '1. Pré-Processamento', 'preprocessing.py', [
    '• Entrada: bytes (UTF-8 ou binário)',
    '• Particionamento em K blocos',
    '• Tamanho fixo de L bytes',
    '• Zero-padding no bloco final',
    '• Registro do tamanho original'
], c_proc_bg, c_proc_border, '#1D4ED8')

# 2. Distribuicao e Soliton
draw_box(ax, 3, 36, 26, 26, '2. DRS e Amostragem', 'soliton.py & config.py', [
    '• Parâmetros c e delta (DRS)',
    '• Cálculo de R e m = floor(K/R)',
    '• Salvaguardas para R > K e R < delta',
    '• Normalização mu(d) = (rho+tau)/beta',
    '• Amostragem via CDF determinística',
    '• Seed PRNG (16 bits) por droplet'
], c_proc_bg, c_proc_border, '#1D4ED8')

# Seta entre 1 e 3
draw_arrow(ax, 29, 80, 36, 80, 'Blocos Fonte (K)')

# 3. Codificador LT
draw_box(ax, 36, 68, 28, 24, '3. Codificação LT & Framing', 'lt_codec.py & utils.py', [
    '• Sorteio do grau d e vizinhos via seed',
    '• Combinação XOR dos d blocos',
    '• Scrambler aditivo (XOR determinístico)',
    '• Cabeçalho: Seed (16b) ofuscada',
    '• Checksum CRC-8 (polinômio 0x07)',
    '• Pacote binário: L + 3 bytes'
], c_proc_bg, c_proc_border, '#1D4ED8')

# Seta entre 2 e 3
draw_arrow(ax, 16, 62, 36, 73, 'DRS / Graus')

# Seta entre 3 e 4
draw_arrow(ax, 50, 68, 50, 58, 'Bytes do Droplet')

# 4. Conversao DNA & Screening
draw_box(ax, 36, 32, 28, 26, '4. Screening Biológico', 'screening.py & dna_conversion.py', [
    '• Mapeamento 2 bits/nt (A, C, G, T)',
    '• Validação GC (40% a 60%)',
    '• Limite de homopolímeros (<= 3 nt)',
    '• Rejeição de candidatos inválidos',
    '• Reposição: nova seed até N_tx válidos'
], c_dna_bg, c_dna_border, '#15803D')

# Seta de retorno do screening (reposicao)
ax.annotate('', xy=(36, 77), xytext=(36, 48),
            arrowprops=dict(arrowstyle="->", color='#EF4444', lw=1.3, linestyle='--', connectionstyle="arc3,rad=0.3"))
ax.text(28.5, 60, 'Rejeitado (nova seed)', ha='center', va='center', fontsize=7, color='#B91C1C', rotation=90)

# Seta entre 4 e 5
draw_arrow(ax, 64, 45, 71, 45, 'N_tx oligos')

# 5. Canal de Apagamento
draw_box(ax, 71, 33, 26, 24, '5. Canal Estocástico', 'channel.py', [
    '• Modelo BEC sem memória (i.i.d.)',
    '• Probabilidade de perda p_loss',
    '• PRNG dedicado (seed independente)',
    '• N_rx ~ Binomial(N_tx, 1 - p_loss)',
    '• Oligos íntegros transmitidos'
], c_chan_bg, c_chan_border, '#B91C1C')

# Seta entre 5 e 6
draw_arrow(ax, 84, 57, 84, 68, 'N_rx oligos')

# 6. Decodificador Peeling
draw_box(ax, 71, 68, 26, 24, '6. Peeling Decoder', 'lt_codec.py', [
    '• Conversão DNA -> bytes & CRC-8',
    '• Extração de seed e regeneração',
    '• Inversão do scrambling aditivo',
    '• Fila de processamento FIFO (deque)',
    '• Propagação XOR e resolução de nós'
], c_dec_bg, c_dec_border, '#7E22CE')

# 7. Metricas e Saida (inferior)
draw_box(ax, 18, 4, 64, 20, '7. Instrumentação, Avaliação Experimental e Reprodutibilidade', 'experiments.py & run_experiments.py', [
    '• Verificação byte a byte contra os blocos fonte originais (taxa de sucesso)',
    '• Rastreamento da dinâmica do Ripple (pico operacional R_max e tamanho da fila)',
    '• Contagem de operações XOR na decodificação e tempo de execução (time.perf_counter)',
    '• Bateria fatorial: 216 combinações (K, c, delta, overhead, p_loss) x 100 repetições',
    '• Exportação estruturada em CSV e geração de figuras vetoriais para documentação'
], '#F8FAFC', '#475569', '#0F172A')

# Setas para o bloco 7
draw_arrow(ax, 84, 68, 84, 24, 'Métricas')
draw_arrow(ax, 50, 32, 50, 24, '')

output_dir = os.path.join('Relatório', 'figuras')
os.makedirs(output_dir, exist_ok=True)
pdf_path = os.path.join(output_dir, 'diagrama_arquitetura.pdf')
png_path = os.path.join(output_dir, 'diagrama_arquitetura.png')
plt.savefig(pdf_path, bbox_inches='tight')
plt.savefig(png_path, bbox_inches='tight', dpi=300)
plt.close()
print(f"Diagrama gerado com sucesso em {pdf_path}")

