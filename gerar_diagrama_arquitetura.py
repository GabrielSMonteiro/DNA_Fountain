import os
import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Figure configuration (12.5 x 8.5 inches, 300 DPI) for clarity and readability
fig, ax = plt.subplots(figsize=(12.5, 8.5), dpi=300)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis('off')

# Professional colors
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
    rect = patches.FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.0,rounding_size=1.5",
        facecolor=bg_col, edgecolor=border_col, linewidth=1.5, zorder=2
    )
    ax.add_patch(rect)
    
    header_h = 3.7
    header_rect = patches.FancyBboxPatch(
        (x + 0.12, y + h - header_h - 0.12), w - 0.24, header_h,
        boxstyle="round,pad=0.0,rounding_size=1.3",
        facecolor=border_col, edgecolor="none", alpha=0.14, zorder=3
    )
    ax.add_patch(header_rect)
    
    ax.text(x + w / 2, y + h - 2.0, title,
            ha='center', va='center', fontsize=8.8, fontweight='bold', color=title_col, zorder=4)
    
    ax.text(x + w / 2, y + h - 4.7, subtitle,
            ha='center', va='center', fontsize=7.2, fontstyle='italic', color='#64748B', zorder=4)
    
    ax.plot([x + 1.2, x + w - 1.2], [y + h - 5.8, y + h - 5.8], color=border_col, alpha=0.22, lw=0.7, zorder=3)
    
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

# Main Title
ax.text(50, 97.5, 'Modular Architecture and Execution Flow of the DNA Fountain Simulator',
        ha='center', va='center', fontsize=12.0, fontweight='bold', color='#0F172A')

# 1. Preprocessing (Col 1, Top)
draw_box(
    ax, x=2, y=71, w=23, h=22,
    title='1. Preprocessing', subtitle='preprocessing.py',
    items=[
        '• Digital payload (UTF-8/bytes)',
        '• Segmentation into K blocks',
        '• Fixed length of L bytes',
        '• Zero-padding on final block',
        '• Store original length metadata'
    ],
    bg_col=colors['bg_prep'], border_col=colors['border_prep'], title_col=colors['title_prep']
)

# 2. RSD & Sampling (Col 1, Middle)
draw_box(
    ax, x=2, y=37, w=23, h=22,
    title='2. RSD & Sampling', subtitle='soliton.py & config.py',
    items=[
        '• Tuning parameters: c and δ',
        '• Distributions: ρ(d) and τ(d)',
        '• Normalization μ(d) and CDF',
        '• Safeguards: m < 1 and R < δ',
        '• Isolated PRNG seed (16 bits)',
        '• Deterministic sampling of d'
    ],
    bg_col=colors['bg_drs'], border_col=colors['border_drs'], title_col=colors['title_drs']
)

# 3. LT Encoding & Framing (Col 2, Top)
draw_box(
    ax, x=38, y=71, w=25, h=22,
    title='3. LT Encoding & Framing', subtitle='lt_codec.py & utils.py',
    items=[
        '• Linear XOR combination: ⨁ b_i',
        '• Additive payload scrambling',
        '• Seed XOR mask (0x1B1B)',
        '• CRC-8 checksum (0x07 poly)',
        '• Framing: [Seed|Payload|CRC]'
    ],
    bg_col=colors['bg_enc'], border_col=colors['border_enc'], title_col=colors['title_enc']
)

# 4. Biological Screening & DNA (Col 3, Top)
draw_box(
    ax, x=76, y=71, w=22, h=22,
    title='4. Screening & DNA', subtitle='screening.py & dna_conversion.py',
    items=[
        '• 2 bits/nt mapping (A,C,G,T)',
        '• GC window: 40% ≤ GC ≤ 60%',
        '• Homopolymer limit: ≤ 3 nt',
        '• Pre-synthesis rejection',
        '• Replacement: seed increment'
    ],
    bg_col=colors['bg_scr'], border_col=colors['border_scr'], title_col=colors['title_scr']
)

# 5. Erasure Channel (Col 3, Middle)
draw_box(
    ax, x=76, y=37, w=22, h=22,
    title='5. Packet Erasure Channel (BEC)', subtitle='channel.py',
    items=[
        '• Packet BEC abstraction',
        '• Bernoulli i.i.d. loss (p_loss)',
        '• Dedicated independent PRNG',
        '• N_rx ~ Binomial(N_tx, 1 - p_loss)',
        '• Surviving oligos arrive intact'
    ],
    bg_col=colors['bg_chan'], border_col=colors['border_chan'], title_col=colors['title_chan']
)

# 6. Peeling Decoder (Col 3, Bottom)
draw_box(
    ax, x=76, y=3, w=22, h=22,
    title='6. Peeling Decoder', subtitle='lt_codec.py',
    items=[
        '• DNA → bytes and CRC-8 check',
        '• Extract seed and neighbors',
        '• Self-inversion descrambling',
        '• FIFO operational queue (deque)',
        '• Iterative XOR propagation',
        '• Monitoring (Q_max, Q_min)'
    ],
    bg_col=colors['bg_dec'], border_col=colors['border_dec'], title_col=colors['title_dec']
)

# 7. Reconstruction & Metrics (Bottom Left & Center)
draw_box(
    ax, x=2, y=3, w=60, h=22,
    title='7. Reconstruction, Experimental Instrumentation & Reproducibility',
    subtitle='preprocessing.py, experiments.py & run_experiments.py',
    items=[
        '• Padding reversal and exact truncation to original file length',
        '• Strict byte-by-byte comparison against source blocks (success rate)',
        '• Metrics logging: success rate, queue peak (Q_max), runtime, XORs',
        '• Factorial grid: 216 combinations (K, c, δ, ε, p_loss) × 1,000 runs (216,000 trials)',
        '• Structured export in aggregated matrices (CSV) and audit logs (JSON)'
    ],
    bg_col=colors['bg_rec'], border_col=colors['border_rec'], title_col=colors['title_rec']
)

# Arrows
draw_arrow(ax, 25, 82.0, 38, 82.0, label='K Source\nBlocks')
draw_arrow(ax, 25, 48.0, 38, 73.0, label='Degree d and\nNeighbors', label_pos=(0.48, 0.48))
draw_arrow(ax, 63, 84.0, 76, 84.0, label='Candidate\n(L + 3)')
draw_arrow(
    ax, 76, 75.0, 63, 75.0,
    label='Rejection\n(seed++)',
    color=colors['arrow_rej'], style='->', lw=1.5, ls='--'
)
draw_arrow(ax, 87, 71.0, 87, 59.0, label='$N_{\\text{tx}}$ Approved\nOligos')
draw_arrow(ax, 87, 37.0, 87, 25.0, label='$N_{\\text{rx}}$ Intact\nOligos')
draw_arrow(ax, 76, 14.0, 62, 14.0, label='Recovered\nBlocks')

output_dir = os.path.join('Relatório', 'figuras')
os.makedirs(output_dir, exist_ok=True)
pdf_path = os.path.join(output_dir, 'diagrama_arquitetura.pdf')
png_path = os.path.join(output_dir, 'diagrama_arquitetura.png')

plt.savefig(pdf_path, bbox_inches='tight')
plt.savefig(png_path, bbox_inches='tight', dpi=300)
plt.close()
print(f"Diagram successfully generated in English at {pdf_path} and {png_path}")
