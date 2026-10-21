import random
from collections import deque

from .config import DNAFountainConfig
from .soliton import RobustSoliton
from .screening import BiologicalScreening
from .dna_conversion import bytes_para_dna, dna_para_bytes
from .utils import xor_blocos, crc8, scramble

class DNAFountainEncoder:

    def __init__(self, config: DNAFountainConfig):
        self.config = config
        self.screening = BiologicalScreening(config)

    def codificar(self, blocos: list[bytes], soliton: RobustSoliton) -> tuple[list[dict], int, list[dict]]:
        oligos = []
        log_rejeicoes = []
        rejeitados = 0
        seed = 0

        while len(oligos) < self.config.num_droplets:
            if seed >= 2 ** self.config.bits_seed:
                raise RuntimeError("Espaço de seeds esgotado.")

            grau, indices = soliton.sortear_grau_e_indices(seed)
            payload = xor_blocos([blocos[i] for i in indices])

            if self.config.scrambler_ativo:
                payload = scramble(payload, seed)

            cabecalho = (seed ^ self.config.seed_xor_mask).to_bytes(self.config.num_bytes_seed, "big")
            dados_oligo = cabecalho + payload

            if self.config.crc_ativo:
                dados_oligo += bytes([crc8(dados_oligo)])

            dna = bytes_para_dna(dados_oligo)

            if self.screening.verificar(dna):
                oligos.append({
                    "seed": seed, "grau": grau, "indices": indices,
                    "dna": dna, "payload": payload,
                })
            else:
                rejeitados += 1
                log_rejeicoes.append({
                    "seed": seed, "grau": grau, "dna": dna, "motivo": "screening"
                })
            seed += 1

        return oligos, rejeitados, log_rejeicoes


class DNAFountainDecoder:

    def __init__(self, config: DNAFountainConfig):
        self.config = config

    def decodificar(self, dnas: list[str], soliton: RobustSoliton) -> tuple[dict, dict]:
        crc_falhas = 0
        pacotes = []

        # Extrair e verificar oligos
        for dna in dnas:
            bruto = dna_para_bytes(dna)
            
            if self.config.crc_ativo:
                crc_recebido = bruto[-1]
                dados = bruto[:-1]
                if crc8(dados) != crc_recebido:
                    crc_falhas += 1
                    continue
            else:
                dados = bruto

            seed = int.from_bytes(dados[:self.config.num_bytes_seed], "big") ^ self.config.seed_xor_mask
            payload = dados[self.config.num_bytes_seed:]

            if self.config.scrambler_ativo:
                payload = scramble(payload, seed)

            _, indices = soliton.sortear_grau_e_indices(seed)
            pacotes.append({"seed": seed, "indices": set(indices), "valor": payload})

        # Índice reverso
        por_bloco = {i: [] for i in range(soliton.K)}
        for p in pacotes:
            for idx in p["indices"]:
                por_bloco[idx].append(p)

        # Peeling
        ripple = deque(p for p in pacotes if len(p["indices"]) == 1)
        descobertos = {}
        total_xors = 0
        ripple_historico = []
        passos = 0

        while ripple:
            ripple_historico.append(len(ripple))
            p = ripple.popleft()

            if len(p["indices"]) != 1: continue
            
            idx_bloco = next(iter(p["indices"]))
            if idx_bloco in descobertos: continue

            descobertos[idx_bloco] = p["valor"]
            passos += 1

            for outro in por_bloco[idx_bloco]:
                if idx_bloco in outro["indices"]:
                    outro["valor"] = xor_blocos([outro["valor"], descobertos[idx_bloco]])
                    total_xors += 1
                    outro["indices"].remove(idx_bloco)
                    if len(outro["indices"]) == 1:
                        ripple.append(outro)

        metricas = {
            "xors": total_xors,
            "passos": passos,
            "ripple_historico": ripple_historico,
            "ripple_max": max(ripple_historico, default=0),
            "ripple_min": min(ripple_historico, default=0),
            "blocos_recuperados": len(descobertos),
            "decodificacao_completa": len(descobertos) == soliton.K,
            "crc_falhas": crc_falhas,
        }

        return descobertos, metricas
