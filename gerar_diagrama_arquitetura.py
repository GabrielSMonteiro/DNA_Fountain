import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Configuração da figura (12.5 x 8.5 polegadas, 300 DPI) para máxima nitidez e legibilidade
fig, ax = plt.subplots(figsize=(12.5, 8.5), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Cores profissionais de alto contraste e harmonia visual
colors = {
    'bg_prep': '#F0F9FF', 'border_prep': '#0284C7', 'title_prep': '#0369A1',
    'bg_drs':  '#F0FDF4', 'border_drs':  '#16A34A', 'title_drs':  '#15803D',
    'bg_enc':  '#EEF2FF', 'border_enc':  '#4F46E5', 'title_enc':  '#3730A3',
    'bg_scr':  '#ECFDF5', 'border_scr':  '#059669', 'title_scr':  '#047857',
    'bg_chan': '#FFF1F2', 'border_chan': '#E11D48', 'title_chan': '#BE123C',
    'bg_dec':  '#FAF5FF', 'border_dec':  '#9333EA', 'title_dec':  '#7E22CE',
    'bg_rec':  '#F8FAFC', 'border_rec':  '#475569', 'title_rec':  '#0F172A',
    'arrow':   '#1E293B', 'arrow_rej':   '#DC2626'
}

def draw_box(ax, x, y, w, h, title, subtitle, items, bg_col, border_col, title_col):
    # Retângulo principal com cantos arredondados
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=1.5",
        facecolor=bg_col, edgecolor=border_col, linewidth=1.5, zorder=2
    )
    ax.add_patch(rect)
    
    # Barra de cabeçalho destacada
    header_h = 3.7
    header_rect = patches.FancyBboxPatch(
        (x + 0.12, y + h - header_h - 0.12), w - 0.24, header_h,
        boxstyle="round,pad=0.0,rounding_size=1.3",
        facecolor=border_col, edgecolor="none", alpha=0.14, zorder=3
    )
    ax.add_patch(header_rect)
    
    # Título do Componente
    ax.text(x + w / 2, y + h - 2.0, title,
            ha='center', va='center', fontsize=8.8, fontweight='bold', color=title_col, zorder=4)
    
    # Nome do módulo / subtítulo
    ax.text(x + w / 2, y + h - 4.7, subtitle,
            ha='center', va='center', fontsize=7.2, fontstyle='italic', color='#64748B', zorder=4)
    
    # Linha divisória sutil
    ax.plot([x + 1.2, x + w - 1.2], [y + h - 5.8, y + h - 5.8], color=border_col, alpha=0.22, lw=0.7, zorder=3)
    
    # Lista de tópicos (bullets) com entrelinha regular e margens amplas
    n_items = len(items)
    avail_h = h - 7.0
    line_h = avail_h / max(n_items, 1)
    
    for i, it in enumerate(items):
        item_y = (y + h - 7.0) - i * line_h - (line_h / 2)
        ax.text(x + 1.2, item_y, it,
                ha='left', va='center', fontsize=7.0, color='#1E293B', zorder=4)

def draw_arrow(ax, x1, y1, x2, y2, label='', label_pos=(0.5, 0.5), color=colors['arrow'], style='->', lw=1.5, ls='-', rad=0.0):
    connectionstyle = f"arc3,rad={rad}" if rad != 0.0 else "arc3,rad=0"
    ax.annotate(
        '', xy=(x2, y2), xytext=(x1, y1),
        arrowprops=dict(
            arrowstyle=style, color=color, lw=lw, linestyle=ls,
            mutation_scale=11, shrinkA=2, shrinkB=2,
            connectionstyle=connectionstyle
        ),
        zorder=5
    )
    if label:
        xm = x1 + (x2 - x1) * label_pos[0]
        ym = y1 + (y2 - y1) * label_pos[1]
        text_color = color if color != colors['arrow'] else '#0F172A'
        border_color = color if color != colors['arrow'] else '#94A3B8'
        ax.text(
            xm, ym, label, ha='center', va='center',
            fontsize=6.5, fontweight='bold', color=text_color,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec=border_color, lw=0.7, alpha=0.98),
            zorder=6
        )

# Título Geral do Diagrama
ax.text(50, 97.5, 'Arquitetura Modular e Fluxo de Execução do Simulador DNA Fountain',
        ha='center', va='center', fontsize=12.0, fontweight='bold', color='#0F172A')

# ==============================================================================
# DIMENSÕES: GRID COM GAPS AMPLOS (X_GAP = 13 UNIDADES, Y_GAP = 12 UNIDADES)
# Coluna 1: x = 2 a 25  (w=23) | Gap: 25 a 38 (13u)
# Coluna 2: x = 38 a 63 (w=25) | Gap: 63 a 76 (13u)
# Coluna 3: x = 76 a 98 (w=22)
# ==============================================================================

# 1. Pré-Processamento (Coluna 1, Linha Superior)
draw_box(
    ax, x=2, y=71, w=23, h=22,
    title='1. Pré-Processamento', subtitle='preprocessing.py',
    items=[
        '• Carga digital (UTF-8/bytes)',
        '• Segmentação em K blocos',
        '• Tamanho fixo de L bytes',
        '• Zero-padding no bloco final',
        '• Guarda do comprimento original'
    ],
    bg_col=colors['bg_prep'], border_col=colors['border_prep'], title_col=colors['title_prep']
)

# 2. DRS e Amostragem (Coluna 1, Linha Média)
draw_box(
    ax, x=2, y=37, w=23, h=22,
    title='2. DRS & Amostragem', subtitle='soliton.py & config.py',
    items=[
        '• Parâmetros de sintonia: c e δ',
        '• Distribuições: ρ(d) e τ(d)',
        '• Normalização μ(d) e CDF',
        '• Salvaguardas: m < 1 e R < δ',
        '• Seed PRNG isolada (16 bits)',
        '• Sorteio determinístico de d'
    ],
    bg_col=colors['bg_drs'], border_col=colors['border_drs'], title_col=colors['title_drs']
)

# 3. Codificação LT e Framing (Coluna 2, Linha Superior)
draw_box(
    ax, x=38, y=71, w=25, h=22,
    title='3. Codificação LT & Framing', subtitle='lt_codec.py & utils.py',
    items=[
        '• Combinação linear XOR: ⨁ b_i',
        '• Scrambling aditivo da carga',
        '• Máscara XOR na seed (0x1B1B)',
        '• Checksum CRC-8 (polinômio 0x07)',
        '• Framing: [Seed|Payload|CRC]'
    ],
    bg_col=colors['bg_enc'], border_col=colors['border_enc'], title_col=colors['title_enc']
)

# 4. Screening Biológico e DNA (Coluna 3, Linha Superior)
draw_box(
    ax, x=76, y=71, w=22, h=22,
    title='4. Screening & DNA', subtitle='screening.py & dna_conversion.py',
    items=[
        '• Mapeamento 2 bits/nt (A,C,G,T)',
        '• Janela GC: 40% ≤ GC ≤ 60%',
        '• Limite homopolímeros: ≤ 3 nt',
        '• Rejeição prévia de inválidos',
        '• Reposição: avanço de seed'
    ],
    bg_col=colors['bg_scr'], border_col=colors['border_scr'], title_col=colors['title_scr']
)

# 5. Canal de Apagamento (Coluna 3, Linha Média)
draw_box(
    ax, x=76, y=37, w=22, h=22,
    title='5. Canal de Apagamento (BEC)', subtitle='channel.py',
    items=[
        '• Abstração BEC por pacote',
        '• Descarte Bernoulli i.i.d. (p_loss)',
        '• PRNG dedicado (independente)',
        '• N_rx ~ Binomial(N_tx, 1 - p_loss)',
        '• Sobreviventes chegam íntegros'
    ],
    bg_col=colors['bg_chan'], border_col=colors['border_chan'], title_col=colors['title_chan']
)

# 6. Peeling Decoder (Coluna 3, Linha Inferior)
draw_box(
    ax, x=76, y=3, w=22, h=22,
    title='6. Peeling Decoder', subtitle='lt_codec.py',
    items=[
        '• DNA → bytes e checagem CRC-8',
        '• Extração de seed e vizinhos',
        '• Inversão do scrambling',
        '• Fila operacional FIFO (deque)',
        '• Propagação XOR iterativa',
        '• Monitoramento (Q_max, Q_min)'
    ],
    bg_col=colors['bg_dec'], border_col=colors['border_dec'], title_col=colors['title_dec']
)

# 7. Reconstrução, Instrumentação & Métricas (Base Esquerda e Central)
draw_box(
    ax, x=2, y=3, w=60, h=22,
    title='7. Reconstrução, Instrumentação Experimental & Reprodutibilidade',
    subtitle='preprocessing.py, experiments.py & run_experiments.py',
    items=[
        '• Reversão do preenchimento e truncamento exato para o comprimento original',
        '• Confronto estrito byte a byte contra os blocos fonte originais (taxa de sucesso)',
        '• Coleta de métricas: taxa de sucesso, pico da fila (Q_max), tempo total e XORs de propagação',
        '• Grade fatorial: 216 combinações (K, c, δ, ε, p_loss) × 100 repetições (21.600 ensaios)',
        '• Exportação estruturada em matrizes agregadas (CSV) e registros de auditoria (JSON)'
    ],
    bg_col=colors['bg_rec'], border_col=colors['border_rec'], title_col=colors['title_rec']
)

# ==============================================================================
# SETAS COM ROTULAGEM EM DUAS LINHAS (PERFEITAMENTE CONTIDAS NOS GAPS LIVRES)
# ==============================================================================

# Seta 1 -> 3: Blocos fonte K (Gap 25 a 38 = 13u)
draw_arrow(ax, 25, 82.0, 38, 82.0, label='K Blocos\nFonte')

# Seta 2 -> 3: Grau d e vizinhança (Diagonal no espaço aberto)
draw_arrow(ax, 25, 48.0, 38, 73.0, label='Grau d e\nVizinhos', label_pos=(0.48, 0.48))

# Seta 3 -> 4 (Ida superior, y=84): Candidato gerado enviado ao screening (Gap 63 a 76 = 13u)
draw_arrow(ax, 63, 84.0, 76, 84.0, label='Candidato\n(L + 3)')

# Seta 4 -> 3 (Retorno inferior, y=74): Feedback de rejeição com incremento de seed (Gap 63 a 76 = 13u)
draw_arrow(
    ax, 76, 75.0, 63, 75.0,
    label='Rejeição\n(seed++)',
    color=colors['arrow_rej'], style='->', lw=1.5, ls='--'
)

# Seta 4 -> 5: Oligos aprovados enviados ao canal (Gap vertical de 12u: 71 a 59)
draw_arrow(ax, 87, 71.0, 87, 59.0, label='$N_{\\text{tx}}$ Oligos\nAprovados')

# Seta 5 -> 6: Oligos que chegam íntegros ao decodificador (Gap vertical de 12u: 37 a 25)
draw_arrow(ax, 87, 37.0, 87, 25.0, label='$N_{\\text{rx}}$ Oligos\nÍntegros')

# Seta 6 -> 7: Blocos reconstruídos enviados à validação (Gap horizontal de 14u: 76 a 62)
draw_arrow(ax, 76, 14.0, 62, 14.0, label='Blocos\nRecuperados')

# Salva nos formatos vetorial (PDF) e rasterizado (PNG)
output_dir = os.path.join('Relatório', 'figuras')
os.makedirs(output_dir, exist_ok=True)
pdf_path = os.path.join(output_dir, 'diagrama_arquitetura.pdf')
png_path = os.path.join(output_dir, 'diagrama_arquitetura.png')

plt.savefig(pdf_path, bbox_inches='tight')
plt.savefig(png_path, bbox_inches='tight', dpi=300)
plt.close()
print(f"Diagrama gerado com sucesso em {pdf_path} e {png_path}")
