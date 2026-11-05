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

