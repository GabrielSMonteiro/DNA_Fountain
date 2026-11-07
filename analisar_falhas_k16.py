"""Auditoria detalhada de falhas de decodificação para K = 16.

Avalia 1.000 repetições independentes com K = 16, c = 0.1, delta = 0.10,
overhead epsilon = 1.00 (N = 32 droplets) e canal ideal (p_loss = 0.0).

Para cada falha do peeling decoder:
1. Verifica se existem blocos-fonte isolados (grau 0 de cobertura inicial).
2. Executa a Eliminação Gaussiana em GF(2) para determinar se a falha decorre
   de posto deficiente (posto < K) ou de um conjunto de parada (stopping set com posto = K).
"""
import json
import random
from dna_fountain.config import DNAFountainConfig
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder
from dna_fountain.gaussian_decoder import GaussianGF2Decoder

def analisar(repeticoes: int = 1000):
    K, c, delta, n_droplets, p_perda = 16, 0.1, 0.1, 32, 0.0
    seed_base = int(K * 100000) + int(c * 10000) + int(delta * 1000) + n_droplets
    config = DNAFountainConfig()
    config.soliton_c = c
    config.soliton_delta = delta
    config.num_droplets = n_droplets
    config.prob_perda = p_perda

    rng_master = random.Random(seed_base)
    rng_payload = random.Random(rng_master.randint(0, 2**31))
    rng_seeds = random.Random(rng_master.randint(0, 2**31))

    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    decoder_peeling = DNAFountainDecoder(config)
    decoder_gauss = GaussianGF2Decoder(config)

    max_seed = 2 ** config.bits_seed
    seed_headroom = max_seed - n_droplets * 2

    falhas = []
    sucessos_peeling = 0
    sucessos_gauss = 0

    for rep in range(repeticoes):
        blocos = [bytes(rng_payload.getrandbits(8) for _ in range(config.tamanho_bloco_l)) for _ in range(K)]
        seed_inicial = rng_seeds.randint(0, max(0, seed_headroom))
        oligos, _, _ = encoder.codificar(blocos, soliton, seed_inicial=seed_inicial)
        dnas = [o["dna"] for o in oligos]

        desc_p, met_p = decoder_peeling.decodificar(dnas, soliton)
        desc_g, met_g = decoder_gauss.decodificar(dnas, soliton)

        if len(desc_p) == K:
            sucessos_peeling += 1
        if met_g["posto_completo"]:
            sucessos_gauss += 1

        if len(desc_p) < K:
            cobertura = set()
            for o in oligos:
                _, indices = soliton.sortear_grau_e_indices(o["seed"])
                cobertura.update(indices)
            isolados = [i for i in range(K) if i not in cobertura]

            # Categorização da falha
            if isolados:
                categoria = "isolamento_inicial_grau0"
            elif met_g["posto_completo"]:
                categoria = "stopping_set_posto_completo"
            else:
                categoria = "dependencia_linear_posto_deficiente"

            falhas.append({
                "repeticao": rep,
                "seed_inicial": seed_inicial,
                "recuperados_peeling": len(desc_p),
                "posto_gaussiano": met_g["posto"],
                "posto_completo": met_g["posto_completo"],
                "blocos_isolados_grau0": isolados,
                "categoria": categoria
            })

    resumo = {
        "configuracao": {"K": K, "c": c, "delta": delta, "N": n_droplets, "overhead": 1.0, "repeticoes": repeticoes},
        "taxa_sucesso_peeling": round(sucessos_peeling / repeticoes, 4),
        "taxa_sucesso_gaussiana": round(sucessos_gauss / repeticoes, 4),
        "total_falhas_peeling": len(falhas),
        "categorias": {
            "stopping_set_posto_completo": sum(1 for f in falhas if f["categoria"] == "stopping_set_posto_completo"),
            "dependencia_linear_posto_deficiente": sum(1 for f in falhas if f["categoria"] == "dependencia_linear_posto_deficiente"),
            "isolamento_inicial_grau0": sum(1 for f in falhas if f["categoria"] == "isolamento_inicial_grau0")
        },
        "detalhe_primeiras_20_falhas": falhas[:20]
    }

    with open("falhas_k16.json", "w", encoding="utf-8") as f:
        json.dump(resumo, f, indent=2, ensure_ascii=False)

    print(f"Auditoria de K=16 concluída ({repeticoes} repetições):")
    print(f"  Sucesso Peeling:    {sucessos_peeling}/{repeticoes} ({sucessos_peeling/repeticoes*100:.1f}%)")
    print(f"  Sucesso Gaussiana:  {sucessos_gauss}/{repeticoes} ({sucessos_gauss/repeticoes*100:.1f}%)")
    print(f"  Total de falhas peeling: {len(falhas)}")
    print(f"    - Stopping sets (posto completo K=16): {resumo['categorias']['stopping_set_posto_completo']}")
    print(f"    - Dependência linear (posto < 16):    {resumo['categorias']['dependencia_linear_posto_deficiente']}")
    print(f"    - Blocos sem cobertura (grau 0):      {resumo['categorias']['isolamento_inicial_grau0']}")
    print(f"  Resultados salvos em falhas_k16.json")

if __name__ == "__main__":
    analisar()
