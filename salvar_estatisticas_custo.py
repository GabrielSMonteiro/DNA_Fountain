import json
import statistics
import time
import random
import os
import matplotlib.pyplot as plt
from dna_fountain.config import DNAFountainConfig
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder
from dna_fountain.channel import ErasureChannel

def rodar_benchmark_k(K, c=0.1, delta=0.1, folga=1.0, p_perda=0.0, repeticoes=100):
    n_droplets = int(K * (1 + folga))
    seed_base = int(K * 100000) + int(c * 10000) + int(delta * 1000) + n_droplets
    
    config = DNAFountainConfig()
    config.soliton_c = c
    config.soliton_delta = delta
    config.num_droplets = n_droplets
    config.prob_perda = p_perda
    
    rng = random.Random(seed_base)
    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    decoder = DNAFountainDecoder(config)
    canal = ErasureChannel(config, rng)
    
    sucessos = 0
    xors = []
    t_dec_ms = []
    
    for _ in range(repeticoes):
        blocos = [bytes(rng.getrandbits(8) for _ in range(config.tamanho_bloco_l)) for _ in range(K)]
        oligos, rejeitados, _ = encoder.codificar(blocos, soliton)
        recebidos, _ = canal.transmitir(oligos)
        
        t0 = time.perf_counter()
        desc, met = decoder.decodificar([o["dna"] for o in recebidos], soliton)
        t_ms = (time.perf_counter() - t0) * 1000.0
        t_dec_ms.append(t_ms)
        
        sucesso = len(desc) == K and all(desc[i] == blocos[i] for i in range(K))
        sucessos += int(sucesso)
        xors.append(met["xors"])
        
    mean_t = statistics.mean(t_dec_ms)
    std_t = statistics.stdev(t_dec_ms)
    med_t = statistics.median(t_dec_ms)
    min_t = min(t_dec_ms)
    max_t = max(t_dec_ms)
    
    mean_xor = statistics.mean(xors)
    xor_per_block = mean_xor / K
    
    return {
        "K": K,
        "N": n_droplets,
        "seed_base": seed_base,
        "repeticoes": repeticoes,
        "sucessos": sucessos,
        "tempos_brutos_ms": t_dec_ms,
        "mean_t": mean_t,
        "std_t": std_t,
        "med_t": med_t,
        "min_t": min_t,
        "max_t": max_t,
        "mean_xor": mean_xor,
        "xor_per_block": xor_per_block
    }

print("Executando baterias para K in (16, 32, 64)...")
res = {}
for K in (16, 32, 64):
    res[K] = rodar_benchmark_k(K)
    print(f"K={K}: Média={res[K]['mean_t']:.2f} ms, Desvio={res[K]['std_t']:.2f} ms, Mediana={res[K]['med_t']:.2f} ms, Min={res[K]['min_t']:.2f} ms, Max={res[K]['max_t']:.2f} ms, XORs={res[K]['mean_xor']:.1f}")

# Salva dados brutos em JSON para reprodutibilidade
with open("tempos_brutos_decodificacao.json", "w", encoding="utf-8") as f:
    json.dump({k: {m: v for m, v in val.items() if m != 'tempos_brutos_ms'} for k, val in res.items()}, f, indent=2)

# Atualiza grafico_custo_computacional.pdf e png
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8.5, 3.8))

ks_num = [16, 32, 64]
tempos_med = [res[k]['mean_t'] for k in ks_num]
tempos_std = [res[k]['std_t'] for k in ks_num]
xors_totais = [res[k]['mean_xor'] for k in ks_num]
xors_por_bloco = [res[k]['xor_per_block'] for k in ks_num]

ax1.errorbar(ks_num, tempos_med, yerr=tempos_std, fmt='-o', color='#1f77b4', ecolor='#aec7e8', elinewidth=2, capsize=5, label='Tempo Médio (±1 desvio-padrão)')
ax1.set_title('Tempo de Decodificação (ms)')
ax1.set_xlabel('Quantidade de Blocos Fonte ($K$)')
ax1.set_ylabel('Tempo (ms)')
ax1.set_xticks(ks_num)
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.legend(frameon=True)

ax2.plot(ks_num, xors_totais, '-s', color='#ff7f0e', label='Total de XORs')
ax2_twin = ax2.twinx()
ax2_twin.plot(ks_num, xors_por_bloco, '--^', color='#2ca02c', label='XORs por Bloco ($/K$)')
ax2.set_title('Operações XOR na Decodificação')
ax2.set_xlabel('Número de Blocos Fonte ($K$)')
ax2.set_ylabel('Total de XORs', color='#ff7f0e')
ax2_twin.set_ylabel('XORs por Bloco Fonte', color='#2ca02c')
ax2.set_xticks(ks_num)
ax2.grid(True, linestyle=':', alpha=0.6)

lines1, labels1 = ax2.get_legend_handles_labels()
lines2, labels2 = ax2_twin.get_legend_handles_labels()
ax2.legend(lines1 + lines2, labels1 + labels2, loc='upper left', frameon=True)

fig.suptitle('Custo Computacional e Escalabilidade do Peeling Decoder ($c=0.1$, $\\delta=0.10$, $\\epsilon=1.00$)', fontsize=12)
output_dir = os.path.join('Relatório', 'figuras')
fig.savefig(os.path.join(output_dir, 'grafico_custo_computacional.pdf'))
fig.savefig(os.path.join(output_dir, 'grafico_custo_computacional.png'), dpi=300)
plt.close(fig)
print("Gráfico atualizado em Relatório/figuras/grafico_custo_computacional.pdf")

