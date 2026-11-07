import pytest
import math
import random
import csv
import io
from dna_fountain.config import DNAFountainConfig
from dna_fountain.dna_conversion import bytes_para_dna, dna_para_bytes
from dna_fountain.soliton import RobustSoliton
from dna_fountain.screening import BiologicalScreening
from dna_fountain.preprocessing import MessageProcessor
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder
from dna_fountain.channel import ErasureChannel
from dna_fountain.experiments import experimento
from dna_fountain.utils import crc8

def test_conversao_bytes_dna_reversibilidade():
    """Valida reversibilidade exata da conversao bytes -> DNA -> bytes para todos os 256 bytes e payload arbitrario."""
    # Todos os 256 bytes
    todos_bytes = bytes(range(256))
    dna = bytes_para_dna(todos_bytes)
    recuperados = dna_para_bytes(dna)
    assert recuperados == todos_bytes

    # Vetor aleatorio
    rng = random.Random(42)
    vetor = bytes(rng.getrandbits(8) for _ in range(1024))
    dna_vetor = bytes_para_dna(vetor)
    rec_vetor = dna_para_bytes(dna_vetor)
    assert rec_vetor == vetor

def test_somatoria_drs():
    """Valida que sum(mu(d)) == 1.0 para multiplas configuracoes de K, c e delta."""
    config = DNAFountainConfig()
    for K in [16, 32, 64, 100]:
        for c in [0.05, 0.1, 0.2]:
            for delta in [0.01, 0.05, 0.1]:
                config.soliton_c = c
                config.soliton_delta = delta
                soliton = RobustSoliton(config, K)
                soma = sum(soliton.mu[1:])
                assert math.isclose(soma, 1.0, rel_tol=1e-5)

def test_screening_gc_limits():
    """Valida que o screening rejeita conteudos GC fora de [40%, 60%]."""
    config = DNAFountainConfig()
    config.gc_min = 0.40
    config.gc_max = 0.60
    config.max_homopolimero = 3
    screening = BiologicalScreening(config)

    # GC = 0% (invalido)
    assert not screening.verificar("ATATATATATATATAT")
    # GC = 100% (invalido)
    assert not screening.verificar("CGCGCGCGCGCGCGCG")
    # GC = 50% sem homopolimeros longos (valido)
    assert screening.verificar("ACGTACGTACGTACGT")

def test_screening_homopolimeros():
    """Valida rejeicao de sequencias com 4 ou mais nucleotideos consecutivos identicos."""
    config = DNAFountainConfig()
    config.gc_min = 0.30
    config.gc_max = 0.70
    config.max_homopolimero = 3
    screening = BiologicalScreening(config)

    # Homopolimeros aceitos (tamanho <= 3)
    assert screening.verificar("AAACCCGGGTTTAAAC")
    # Homopolimeros rejeitados (tamanho >= 4)
    assert not screening.verificar("AAAACCCGGGTTTAAA")
    assert not screening.verificar("ACGTCCCCACGTACGT")
    assert not screening.verificar("ACGTACGTGGGGACGT")
    assert not screening.verificar("ACGTACGTTTTTACGT")

def test_peeling_canal_ideal():
    """Valida recuperacao 100% exata dos blocos fonte pelo peeling decoder em canal ideal com overhead adequado."""
    config = DNAFountainConfig()
    config.tamanho_bloco_l = 8
    config.soliton_c = 0.1
    config.soliton_delta = 0.1
    K = 32
    config.num_droplets = int(K * 2.0)  # overhead = 1.0 (N = 64)

    rng = random.Random(12345)
    blocos = [bytes(rng.getrandbits(8) for _ in range(config.tamanho_bloco_l)) for _ in range(K)]

    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    decoder = DNAFountainDecoder(config)

    oligos, _, _ = encoder.codificar(blocos, soliton)
    dnas = [o["dna"] for o in oligos]

    descobertos, metricas = decoder.decodificar(dnas, soliton)
    assert len(descobertos) == K
    for i in range(K):
        assert descobertos[i] == blocos[i]

def test_segmentacao_padding_truncamento():
    """RF01: Valida que mensagens com comprimento nao multiplo de L recebem zero-padding e sao truncadas exatamente na reconstrucao."""
    config = DNAFountainConfig()
    config.tamanho_bloco_l = 8
    config.padding_aleatorio = False
    proc = MessageProcessor(config)
    
    # 25 bytes (nao multiplo de 8 -> ceil(25/8) = 4 blocos de 8 bytes, ultimo com 7 zeros de padding)
    msg_bytes = b"Mensagem de teste de 25 b"
    assert len(msg_bytes) == 25
    
    blocos = proc.segmentar(msg_bytes)
    assert len(blocos) == 4
    for b in blocos:
        assert len(b) == 8
    assert blocos[3] == msg_bytes[24:25] + b"\x00" * 7
    
    # Simulando decodificacao completa
    descobertos = {i: blocos[i] for i in range(4)}
    recuperada = proc.reconstruir(descobertos, K=4, tamanho_original=len(msg_bytes))
    assert recuperada == msg_bytes

def test_amostragem_deterministica_seed():
    """RF02: Valida que a mesma seed gera deterministicamente o mesmo grau e os mesmos indices de blocos fonte na DRS."""
    config = DNAFountainConfig()
    config.soliton_c = 0.1
    config.soliton_delta = 0.10
    K = 32
    soliton = RobustSoliton(config, K)
    
    for seed in [0, 42, 100, 2024, 65535]:
        grau1, ind1 = soliton.sortear_grau_e_indices(seed)
        grau2, ind2 = soliton.sortear_grau_e_indices(seed)
        assert grau1 == grau2
        assert ind1 == ind2
        assert 1 <= grau1 <= K
        assert len(ind1) == grau1
        assert len(set(ind1)) == grau1  # sem indices duplicados
        assert all(0 <= idx < K for idx in ind1)

def test_consistencia_vizinhancas_encoder_decoder():
    """RF03: Valida que os indices amostrados pelo codificador coincidem estritamente com os regenerados pelo decodificador a partir da seed."""
    config = DNAFountainConfig()
    config.tamanho_bloco_l = 8
    K = 16
    config.num_droplets = 20
    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    
    blocos = [bytes([i] * config.tamanho_bloco_l) for i in range(K)]
    oligos, _, _ = encoder.codificar(blocos, soliton)
    
    for oligo in oligos:
        seed = oligo["seed"]
        grau_dec, ind_dec = soliton.sortear_grau_e_indices(seed)
        assert grau_dec == oligo["grau"]
        assert ind_dec == oligo["indices"]

def test_canal_apagamento_extremos():
    """RF06: Valida o comportamento do canal nos limites teoricos p_loss = 0.0 (preservacao total) e p_loss = 1.0 (apagamento total)."""
    config = DNAFountainConfig()
    oligos = [{"seed": i, "dna": "ACGT"} for i in range(50)]
    
    # p_loss = 0.0
    config.prob_perda = 0.0
    canal_sem_perda = ErasureChannel(config, rng=random.Random(123))
    recebidos, perdidos = canal_sem_perda.transmitir(oligos)
    assert len(recebidos) == 50
    assert len(perdidos) == 0
    assert recebidos == oligos
    
    # p_loss = 1.0
    config.prob_perda = 1.0
    canal_perda_total = ErasureChannel(config, rng=random.Random(123))
    recebidos, perdidos = canal_perda_total.transmitir(oligos)
    assert len(recebidos) == 0
    assert len(perdidos) == 50
    assert perdidos == oligos

def test_exportacao_metricas_csv():
    """RF08: Valida que a execucao experimental consolida e exporta corretamente as metricas com as colunas esperadas."""
    res = experimento(K=16, c=0.1, delta=0.1, n_droplets=32, p_perda=0.0, repeticoes=2, seed_base=42)
    campos_esperados = [
        "K", "c", "delta", "N", "p_perda", "overhead_enviado",
        "R", "grau_medio", "taxa_sucesso", "xors_medio",
        "ripple_max_medio", "ripple_min_medio", "tempo_decod_ms"
    ]
    for c in campos_esperados:
        assert c in res
        assert isinstance(res[c], (int, float))
    assert res["taxa_sucesso"] in (0.0, 0.5, 1.0)
    assert res["N"] == 32
    assert res["K"] == 16
    
    # Valida gravacao e leitura de CSV
    buffer = io.StringIO()
    writer = csv.DictWriter(buffer, fieldnames=campos_esperados, delimiter=";")
    writer.writeheader()
    writer.writerow(res)
    buffer.seek(0)
    reader = csv.DictReader(buffer, delimiter=";")
    linhas = list(reader)
    assert len(linhas) == 1
    assert int(linhas[0]["K"]) == 16
    assert float(linhas[0]["overhead_enviado"]) == 1.0

def test_somatoria_drs_salvaguarda_m_menor_1():
    """RF02: Valida que a DRS com valores elevados de c (c >= 0.5) ativa a salvaguarda m_raw < 1 -> m = 1 e mantem sum(mu) == 1.0."""
    config = DNAFountainConfig()
    casos = [
        (16, 0.8, 0.10),  # R ~ 16.24 > 16 -> m_raw = 0 -> m = 1
        (16, 0.8, 0.01),  # R ~ 22.75 > 16 -> m_raw = 0 -> m = 1
        (32, 0.8, 0.01),  # R ~ 36.52 > 32 -> m_raw = 0 -> m = 1
        (32, 0.5, 0.10),  # R ~ 16.32 > 16 -> m = 1
    ]
    for K, c, delta in casos:
        config.soliton_c = c
        config.soliton_delta = delta
        soliton = RobustSoliton(config, K)
        assert math.isclose(sum(soliton.mu[1:]), 1.0, rel_tol=1e-5)
        # No regime c >= 0.5 com m=1, a massa corretiva concentra no grau 1
        assert soliton.mu[1] > 0.50

def test_grafo_canonico_gotas_unitarias():
    """RF07: Valida o peeling decoder sobre grafo canonico composto estritamente por K gotas de grau 1."""
    config = DNAFountainConfig()
    config.crc_ativo = True
    config.scrambler_ativo = False
    K = 8
    L = 4
    config.tamanho_bloco_l = L
    
    blocos_esperados = [bytes([i + 1] * L) for i in range(K)]
    dnas = []
    mapping = {}
    for i in range(K):
        seed = i
        cabecalho = (seed ^ config.seed_xor_mask).to_bytes(2, "big")
        dados = cabecalho + blocos_esperados[i]
        dados += bytes([crc8(dados)])
        dnas.append(bytes_para_dna(dados))
        mapping[seed] = (1, [i])

    class MockSoliton:
        def __init__(self, K, mapping):
            self.K = K
            self.mapping = mapping
        def sortear_grau_e_indices(self, seed):
            return self.mapping[seed]

    sol = MockSoliton(K, mapping)
    dec = DNAFountainDecoder(config)
    desc, met = dec.decodificar(dnas, sol)

    assert len(desc) == K
    assert met["decodificacao_completa"] is True
    assert met["ripple_max"] == K
    assert all(desc[i] == blocos_esperados[i] for i in range(K))

def test_grafo_canonico_cascata_triangular():
    """RF07: Valida o peeling decoder sobre grafo canonico em cascata triangular (gotas encadeadas d_i = {i-1, i})."""
    config = DNAFountainConfig()
    config.crc_ativo = True
    config.scrambler_ativo = False
    K = 6
    L = 4
    config.tamanho_bloco_l = L

    blocos = [bytes([i + 1] * L) for i in range(K)]
    dnas = []
    mapping = {}
    for i in range(K):
        seed = i
        if i == 0:
            indices = [0]
            payload = blocos[0]
        else:
            indices = [i - 1, i]
            payload = bytes(a ^ b for a, b in zip(blocos[i - 1], blocos[i]))

        cabecalho = (seed ^ config.seed_xor_mask).to_bytes(2, "big")
        dados = cabecalho + payload
        dados += bytes([crc8(dados)])
        dnas.append(bytes_para_dna(dados))
        mapping[seed] = (len(indices), indices)

    class MockSoliton:
        def __init__(self, K, mapping):
            self.K = K
            self.mapping = mapping
        def sortear_grau_e_indices(self, seed):
            return self.mapping[seed]

    sol = MockSoliton(K, mapping)
    dec = DNAFountainDecoder(config)
    desc, met = dec.decodificar(dnas, sol)

    assert len(desc) == K
    assert met["decodificacao_completa"] is True
    assert met["ripple_max"] == 1
    assert met["ripple_min"] == 1
    assert all(desc[i] == blocos[i] for i in range(K))

def test_falha_stopping_set_sem_grau1():
    """RNF: Valida parada imediata do peeling decoder na presenca de um stopping set sem gotas de grau 1."""
    config = DNAFountainConfig()
    config.crc_ativo = True
    config.scrambler_ativo = False
    K = 3
    L = 4
    config.tamanho_bloco_l = L
    blocos = [bytes([i + 1] * L) for i in range(K)]

    arestas = [[0, 1], [1, 2], [2, 0]]
    dnas = []
    mapping = {}
    for seed, ind in enumerate(arestas):
        payload = bytes(a ^ b for a, b in zip(blocos[ind[0]], blocos[ind[1]]))
        cabecalho = (seed ^ config.seed_xor_mask).to_bytes(2, "big")
        dados = cabecalho + payload
        dados += bytes([crc8(dados)])
        dnas.append(bytes_para_dna(dados))
        mapping[seed] = (2, ind)

    class MockSoliton:
        def __init__(self, K, mapping):
            self.K = K
            self.mapping = mapping
        def sortear_grau_e_indices(self, seed):
            return self.mapping[seed]

    sol = MockSoliton(K, mapping)
    dec = DNAFountainDecoder(config)
    desc, met = dec.decodificar(dnas, sol)

    assert len(desc) == 0
    assert met["decodificacao_completa"] is False
    assert met["ripple_max"] == 0
    assert met["passos"] == 0

def test_descarte_crc_corrompido():
    """RF03/RF07: Valida descarte obrigatorio de oligonucleotideos com CRC invalido por erro de sequenciamento simulado."""
    config = DNAFountainConfig()
    config.tamanho_bloco_l = 4
    K = 4
    config.num_droplets = 8
    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    decoder = DNAFountainDecoder(config)

    blocos = [bytes([i] * 4) for i in range(K)]
    oligos, _, _ = encoder.codificar(blocos, soliton)

    dnas_corrompidos = []
    for o in oligos:
        mutado = o["dna"][:-1] + ("A" if o["dna"][-1] != "A" else "C")
        dnas_corrompidos.append(mutado)

    desc, met = decoder.decodificar(dnas_corrompidos, soliton)
    assert met["crc_falhas"] == len(oligos)
    assert len(desc) == 0
    assert met["decodificacao_completa"] is False

def test_esgotamento_espaco_seeds():
    """RNF: Valida lancamento de excecao quando o espaco de enderecamento de seeds de 16 bits e esgotado."""
    config = DNAFountainConfig()
    config.num_droplets = 10
    soliton = RobustSoliton(config, K=4)
    encoder = DNAFountainEncoder(config)
    blocos = [bytes([i] * 4) for i in range(4)]

    # Inicia a 3 passos do limite do espaco de 16 bits (65536)
    with pytest.raises(RuntimeError, match="Espaço de seeds esgotado"):
        encoder.codificar(blocos, soliton, seed_inicial=65533)

def test_validacao_cruzada_peeling_gaussiana():
    """Validacao cruzada: o algoritmo de peeling e verificado contra a Eliminacao Gaussiana em GF(2)."""
    from dna_fountain.gaussian_decoder import GaussianGF2Decoder

    config = DNAFountainConfig()
    config.tamanho_bloco_l = 4
    config.soliton_c = 0.1
    config.soliton_delta = 0.1
    K = 16
    config.num_droplets = 32

    soliton = RobustSoliton(config, K)
    encoder = DNAFountainEncoder(config)
    decoder_peeling = DNAFountainDecoder(config)
    decoder_gauss = GaussianGF2Decoder(config)

    rng = random.Random(999)
    for seed_exp in [42, 100, 2024]:
        blocos = [bytes(rng.getrandbits(8) for _ in range(4)) for _ in range(K)]
        oligos, _, _ = encoder.codificar(blocos, soliton, seed_inicial=seed_exp * 10)
        dnas = [o["dna"] for o in oligos]

        desc_p, met_p = decoder_peeling.decodificar(dnas, soliton)
        desc_g, met_g = decoder_gauss.decodificar(dnas, soliton)

        # Propriedade 1: Posto completo e condicao necessaria para decodificacao completa
        if met_p["decodificacao_completa"]:
            assert met_g["posto_completo"] is True
            assert met_g["decodificacao_completa"] is True
            assert all(desc_p[i] == desc_g[i] == blocos[i] for i in range(K))

        # Propriedade 2: Posto da Gaussiana limita superiormente os blocos recuperados
        assert met_p["blocos_recuperados"] <= met_g["posto"]


