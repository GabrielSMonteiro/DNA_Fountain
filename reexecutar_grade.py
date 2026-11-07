"""Reexecução paralela da grade fatorial completa com seeds independentes por repetição.

Gera:
- resultados_1000.csv (grade principal com 1000 repetições por ponto = 216.000 ensaios)
- resultados_100.csv (cópia mantida para compatibilidade direta de scripts existentes)
- tempos_brutos_decodificacao.json (benchmark atualizado de custo computacional)
"""
import csv
import json
import os
import sys
import time
import random
import statistics
import concurrent.futures

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from dna_fountain.experiments import experimento
from dna_fountain.config import DNAFountainConfig
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder

REPETICOES = 1000
MAX_WORKERS = 8

def worker_tarefa(params):
    K, c, delta, folga, p_perda, repeticoes = params
    n_droplets = int(K * (1 + folga))
    t0 = time.time()
    res = experimento(K, c, delta, n_droplets, p_perda, repeticoes)
    dt = time.time() - t0
    return (params, res, dt)

def main():
    print(f"{'='*70}")
    print(f"  REEXECUÇÃO DA GRADE FATORIAL — SEEDS INDEPENDENTES")
    print(f"  {REPETICOES} repetições por ponto × 216 combinações")
    print(f"  = {REPETICOES * 216:,} ensaios no total (Workers: {MAX_WORKERS})")
    print(f"{'='*70}\n", flush=True)

    t_inicio = time.time()

    # Prepara lista de tarefas
    tarefas = []
    for K in (16, 32, 64):
        for c in (0.1, 0.3, 0.5, 0.8):
            for delta in (0.1, 0.01):
                for folga in (0.25, 0.5, 1.0):
                    for p_perda in (0.0, 0.1, 0.2):
                        tarefas.append((K, c, delta, folga, p_perda, REPETICOES))

    total = len(tarefas)
    concluidos = 0
    resultados_dict = {}

    with concurrent.futures.ProcessPoolExecutor(max_workers=MAX_WORKERS) as executor:
        futs = {executor.submit(worker_tarefa, t): t for t in tarefas}
        for fut in concurrent.futures.as_completed(futs):
            params, res, dt = fut.result()
            concluidos += 1
            K, c, delta, folga, p_perda, _ = params
            chave = (K, c, delta, folga, p_perda)
            resultados_dict[chave] = res
            
            if concluidos % 10 == 0 or concluidos == total:
                tempo_decorrido = time.time() - t_inicio
                taxa = concluidos / tempo_decorrido if tempo_decorrido > 0 else 1
                tempo_restante = (total - concluidos) / taxa
                print(f"  [{concluidos:3d}/{total}] Concluído {concluidos/total:6.1%} | "
                      f"Último: K={K:2d} c={c:.1f} δ={delta:.2f} ε={folga:.2f} p={p_perda:.1f} → {res['taxa_sucesso']*100:5.1f}% | "
                      f"Tempo decorrido: {tempo_decorrido/60:.1f}m (est. restante: {tempo_restante/60:.1f}m)",
                      flush=True)

    # Reordenar resultados na ordem canônica da grade
    linhas = []
    for t in tarefas:
        K, c, delta, folga, p_perda, _ = t
        linhas.append(resultados_dict[(K, c, delta, folga, p_perda)])

    # Salvar em resultados_1000.csv e resultados_100.csv
    for arq in ("resultados_1000.csv", "resultados_100.csv"):
        with open(arq, "w", newline="", encoding="utf-8-sig") as f:
            w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()), delimiter=',')
            w.writeheader()
            w.writerows(linhas)
        print(f"  Resultados salvos com sucesso em {arq}")

    t_total = time.time() - t_inicio
    print(f"\n{'='*70}")
    print(f"  Grade completa em {t_total:.1f}s ({t_total/60:.2f} min)")
    print(f"{'='*70}\n", flush=True)

    # ─────────────────────────────────────────────────────────────
    # Benchmark de custo computacional com seeds independentes
    # ─────────────────────────────────────────────────────────────
    print(f"  Recalculando benchmark de custo computacional com {REPETICOES} repetições...")
    json_path = "tempos_brutos_decodificacao.json"
    bench = {}

    for K in (16, 32, 64):
        n_droplets = int(K * 2.0)  # overhead = 1.00
        config = DNAFountainConfig()
        config.soliton_c = 0.1
        config.soliton_delta = 0.1
        config.num_droplets = n_droplets
        config.prob_perda = 0.0

        seed_base = int(K * 100000) + 1000 + n_droplets
        rng_master = random.Random(seed_base)
        rng_payload = random.Random(rng_master.randint(0, 2**31))
        rng_seeds = random.Random(rng_master.randint(0, 2**31))

        soliton = RobustSoliton(config, K)
        encoder = DNAFountainEncoder(config)
        decoder = DNAFountainDecoder(config)
        max_seed = 2 ** config.bits_seed
        seed_headroom = max_seed - n_droplets * 2

        t_dec_ms = []
        xors_list = []
        sucessos = 0

        for _ in range(REPETICOES):
            blocos = [bytes(rng_payload.getrandbits(8) for _ in range(config.tamanho_bloco_l)) for _ in range(K)]
            seed_inicial = rng_seeds.randint(0, max(0, seed_headroom))
            oligos, _, _ = encoder.codificar(blocos, soliton, seed_inicial=seed_inicial)

            t0 = time.perf_counter()
            desc, met = decoder.decodificar([o["dna"] for o in oligos], soliton)
            t_ms = (time.perf_counter() - t0) * 1000.0
            t_dec_ms.append(t_ms)
            xors_list.append(met["xors"])

            sucesso = (len(desc) == K and all(desc[i] == blocos[i] for i in range(K)))
            sucessos += int(sucesso)

        mean_t = statistics.mean(t_dec_ms)
        std_t = statistics.stdev(t_dec_ms)
        med_t = statistics.median(t_dec_ms)
        min_t = min(t_dec_ms)
        max_t = max(t_dec_ms)
        mean_xor = statistics.mean(xors_list)
        xor_per_block = mean_xor / K

        bench[K] = {
            "K": K,
            "N": n_droplets,
            "repeticoes": REPETICOES,
            "sucessos": sucessos,
            "taxa_sucesso": round(sucessos / REPETICOES, 4),
            "mean_t": round(mean_t, 4),
            "std_t": round(std_t, 4),
            "med_t": round(med_t, 4),
            "min_t": round(min_t, 4),
            "max_t": round(max_t, 4),
            "mean_xor": round(mean_xor, 1),
            "xor_per_block": round(xor_per_block, 2)
        }
        print(f"    K={K:2d}: tempo_médio={mean_t:.2f}ms (desvio={std_t:.2f}ms) | XORs_médio={mean_xor:.1f} | sucesso={sucessos}/{REPETICOES} ({sucessos/REPETICOES*100:.1f}%)")

    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(bench, f, indent=2)
    print(f"  Benchmark salvo com sucesso em {json_path}")

if __name__ == "__main__":
    main()
