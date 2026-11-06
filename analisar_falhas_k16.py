import json
import random
from dna_fountain.config import DNAFountainConfig
from dna_fountain.soliton import RobustSoliton
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder
from dna_fountain.channel import ErasureChannel

def analisar():
    K, c, delta, n_droplets, p_perda, repeticoes = 16, 0.1, 0.1, 32, 0.0, 100
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

    falhas = []
    for rep in range(repeticoes):
        blocos = [bytes(rng.getrandbits(8) for _ in range(config.tamanho_bloco_l)) for _ in range(K)]
        oligos, _, _ = encoder.codificar(blocos, soliton)
        recebidos, _ = canal.transmitir(oligos)
        desc, met = decoder.decodificar([o['dna'] for o in recebidos], soliton)
        sucesso = len(desc) == K
        if not sucesso:
            cobertura = set()
            for o in recebidos:
                _, indices = soliton.sortear_grau_e_indices(o['seed'])
                cobertura.update(indices)
            isolados = [i for i in range(K) if i not in cobertura]
            falhas.append({
                'repeticao': rep,
                'seed_base_cenario': seed_base,
                'recuperados': len(desc),
                'blocos_isolados_grau0': isolados,
                'categoria': 'isolamento_inicial' if isolados else 'stopping_set_esgotamento_ripple'
            })

    with open('falhas_k16.json', 'w', encoding='utf-8') as f:
        json.dump(falhas, f, indent=2, ensure_ascii=False)
    print(f"Total de falhas registradas: {len(falhas)} salvas em falhas_k16.json")

if __name__ == '__main__':
    analisar()
