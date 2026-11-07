"""Decodificador exato por Eliminação Gaussiana sobre o corpo finito GF(2).

Usado para validação cruzada do simulador DNA Fountain e determinação do limite
superior analítico de decodificabilidade (posto matricial = K).

Propriedades fundamentais:
- O algoritmo de peeling é uma heurística de tempo linear O(N) que falha na presença
  de conjuntos de parada (stopping sets).
- A eliminação gaussiana em GF(2) é um solucionador exato O(K^2 * N) que encontra
  a solução sempre que a matriz de incidência possui posto completo (posto = K).
- Se a eliminação gaussiana decodifica mas o peeling falha, a falha decorre
  estritamente de um stopping set estrutural no grafo de conexões.
- Se a eliminação gaussiana falha (posto < K), os pacotes recebidos não cobrem
  o espaço linear dos blocos-fonte (falha de cobertura ou redundância insuficiente).
"""

from .config import DNAFountainConfig
from .soliton import RobustSoliton
from .dna_conversion import dna_para_bytes
from .utils import crc8, scramble


class GaussianGF2Decoder:
    """Decodificador ótimo por eliminação de Gauss-Jordan em GF(2)."""

    def __init__(self, config: DNAFountainConfig):
        self.config = config

    def decodificar(self, dnas: list[str], soliton: RobustSoliton) -> tuple[dict, dict]:
        """Executa a decodificação por eliminação gaussiana em GF(2).

        Args:
            dnas: lista de sequências de DNA recebidas.
            soliton: gerador da DRS configurado para o K em questão.

        Returns:
            descobertos: dicionário {idx_bloco: bytes} com os blocos recuperados.
            metricas: dicionário com estatísticas (posto, posto_completo, etc.).
        """
        K = soliton.K
        L = self.config.tamanho_bloco_l
        crc_falhas = 0
        linhas_matriz = []
        vetor_payloads = []

        # 1. Parse dos oligonucleotídeos válidos
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

            # Representa a linha da matriz como um bitmask inteiro de K bits
            # Bit i ligado se o bloco i participa do droplet
            mask = 0
            for idx in indices:
                mask |= (1 << idx)

            linhas_matriz.append(mask)
            vetor_payloads.append(bytearray(payload))

        M = len(linhas_matriz)

        # 2. Eliminação Gaussiana em GF(2) (Gauss-Jordan)
        # Montamos a matriz aumentada [A | y]
        matriz_A = list(linhas_matriz)
        matriz_Y = list(vetor_payloads)

        pivos = {}  # col_bloco -> row_idx
        linha_atual = 0

        for col in range(K):
            col_mask = (1 << col)
            # Procura pivô na coluna atual a partir de linha_atual
            pivo_row = -1
            for r in range(linha_atual, M):
                if matriz_A[r] & col_mask:
                    pivo_row = r
                    break

            if pivo_row == -1:
                # Coluna linearmente dependente ou descoberta
                continue

            # Troca linha atual com a linha do pivô
            matriz_A[linha_atual], matriz_A[pivo_row] = matriz_A[pivo_row], matriz_A[linha_atual]
            matriz_Y[linha_atual], matriz_Y[pivo_row] = matriz_Y[pivo_row], matriz_Y[linha_atual]

            # Elimina todas as outras linhas (acima e abaixo para Gauss-Jordan)
            pivot_A = matriz_A[linha_atual]
            pivot_Y = matriz_Y[linha_atual]

            for r in range(M):
                if r != linha_atual and (matriz_A[r] & col_mask):
                    matriz_A[r] ^= pivot_A
                    # XOR nos bytes do payload
                    for byte_i in range(L):
                        matriz_Y[r][byte_i] ^= pivot_Y[byte_i]

            pivos[col] = linha_atual
            linha_atual += 1
            if linha_atual >= M:
                break

        posto = len(pivos)
        posto_completo = (posto == K)

        # 3. Extrai blocos recuperados (aqueles cujo vetor unitário isolou um único bloco)
        descobertos = {}
        for col, row in pivos.items():
            # Se a linha resultante contém exatamente 1 bit (o da coluna do pivô)
            if matriz_A[row] == (1 << col):
                descobertos[col] = bytes(matriz_Y[row])

        metricas = {
            "posto": posto,
            "posto_completo": posto_completo,
            "blocos_recuperados": len(descobertos),
            "decodificacao_completa": len(descobertos) == K,
            "total_gotas_validas": M,
            "crc_falhas": crc_falhas
        }

        return descobertos, metricas

