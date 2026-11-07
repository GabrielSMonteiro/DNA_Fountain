import csv
import random
import time
from .config import DNAFountainConfig
from .preprocessing import MessageProcessor
from .soliton import RobustSoliton
from .lt_codec import DNAFountainEncoder, DNAFountainDecoder
from .channel import ErasureChannel

def experimento(K, c, delta, n_droplets, p_perda, repeticoes, seed_base=None):
    """Run a batch of independent repetitions for one factorial grid point.

    Each repetition draws its own random initial seed for the encoder,
    ensuring that the bipartite graph structure varies across runs and
    the measured success rate reflects the LT code ensemble rather than
    a single deterministic seed block.

    Three independent PRNG streams are maintained:
      - rng_payload : generates random source blocks (payload data)
      - rng_seeds   : draws the encoder's starting seed per repetition
      - rng_channel : drives the stochastic erasure channel

    Including p_perda in seed_base guarantees that different loss-rate
    columns are not implicitly paired (they use distinct PRNG states).
    """
    if seed_base is None:
        seed_base = (int(K * 100000)
                     + int(c * 10000)
                     + int(delta * 1000)
                     + n_droplets
                     + int(p_perda * 100))

    # Master RNG only used to derive independent child seeds
    rng_master = random.Random(seed_base)
    child_seed_payload = rng_master.randint(0, 2**31)
    child_seed_seeds   = rng_master.randint(0, 2**31)
    child_seed_channel = rng_master.randint(0, 2**31)

    rng_payload = random.Random(child_seed_payload)
    rng_seeds   = random.Random(child_seed_seeds)
    rng_channel = random.Random(child_seed_channel)

    config = DNAFountainConfig()
    config.soliton_c = c
    config.soliton_delta = delta
    config.num_droplets = n_droplets
    config.prob_perda = p_perda

    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    decoder = DNAFountainDecoder(config)
    canal = ErasureChannel(config, rng_channel)

    max_seed = 2 ** config.bits_seed
    # Reserve enough room so that seed_inicial + num_droplets*2 (worst case
    # with ~40% screening rejection) stays within the 16-bit seed space.
    seed_headroom = max_seed - n_droplets * 2

    sucessos, xors, rmax, rmin, t_dec = 0, [], [], [], []

    for _ in range(repeticoes):
        # Independent random payload per repetition
        blocos = [bytes(rng_payload.getrandbits(8)
                        for _ in range(config.tamanho_bloco_l))
                  for _ in range(K)]

        # Independent random starting seed per repetition
        seed_inicial = rng_seeds.randint(0, max(0, seed_headroom))

        oligos, rejeitados, _ = encoder.codificar(blocos, soliton,
                                                  seed_inicial=seed_inicial)
        recebidos, _ = canal.transmitir(oligos)

        t0 = time.perf_counter()
        desc, met = decoder.decodificar([o["dna"] for o in recebidos], soliton)
        t_dec.append(time.perf_counter() - t0)

        sucesso = (len(desc) == K
                   and all(desc[i] == blocos[i] for i in range(K)))
        sucessos += int(sucesso)

        xors.append(met["xors"])
        rmax.append(met["ripple_max"])
        rmin.append(met["ripple_min"])

    n = repeticoes
    return {
        "K": K,
        "c": c,
        "delta": delta,
        "N": n_droplets,
        "p_perda": p_perda,
        "overhead_enviado": round((n_droplets - K) / K, 3),
        "R": round(soliton.R, 3),
        "grau_medio": round(soliton.grau_medio, 3),
        "taxa_sucesso": sucessos / n,
        "xors_medio": sum(xors) / n if n else 0,
        "ripple_max_medio": sum(rmax) / n if n else 0,
        "ripple_min_medio": sum(rmin) / n if n else 0,
        "tempo_decod_ms": 1000 * sum(t_dec) / n if n else 0
    }

def rodar_experimentos(arquivo="resultados_experimentos.csv", repeticoes=10):
    linhas = []
    print(f"Iniciando simulações para {arquivo} (repetições: {repeticoes})...")

    for K in (16, 32, 64):
        for c in (0.1, 0.3, 0.5, 0.8):
            for delta in (0.1, 0.01):
                for folga in (0.25, 0.5, 1.0):
                    for p_perda in (0.0, 0.1, 0.2):
                        linhas.append(experimento(K, c, delta,
                                                  int(K * (1 + folga)),
                                                  p_perda, repeticoes))

    with open(arquivo, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()), delimiter=',')
        w.writeheader()
        w.writerows(linhas)

    print(f"Resultados salvos em {arquivo}.")
    return linhas
