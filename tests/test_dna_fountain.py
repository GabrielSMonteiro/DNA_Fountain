import pytest
import math
import random
from dna_fountain.config import DNAFountainConfig
from dna_fountain.dna_conversion import bytes_para_dna, dna_para_bytes
from dna_fountain.soliton import RobustSoliton
from dna_fountain.screening import BiologicalScreening
from dna_fountain.preprocessing import MessageProcessor
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder

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
