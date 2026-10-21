from dataclasses import dataclass

@dataclass
class DNAFountainConfig:
    # 1. Pré-processamento
    tamanho_bloco_l: int = 4
    padding_aleatorio: bool = True
    padding_seed: int = 7777

    # 2. Distribuição Robusta de Sóliton
    soliton_c: float = 0.5
    soliton_delta: float = 0.1

    # 3. Codificação LT
    num_droplets: int = 30
    bits_seed: int = 16
    seed_xor_mask: int = 0x1B1B

    # 4. Scrambler e CRC
    scrambler_ativo: bool = True
    scrambler_offset: int = 0x5A5A
    crc_ativo: bool = True
    crc_polinomio: int = 0x07

    # 5. Screening Biológico
    gc_min: float = 0.40
    gc_max: float = 0.60
    max_homopolimero: int = 3

    # 6. Canal e Experimentos
    prob_perda: float = 0.25
    seed_experimento: int = 42

    @property
    def num_bytes_seed(self) -> int:
        return self.bits_seed // 8
