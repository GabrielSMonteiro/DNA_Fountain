import random

class C:
    RESET = '\033[0m'
    RED = '\033[91m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    BLUE = '\033[94m'
    MAGENTA = '\033[95m'
    CYAN = '\033[96m'
    BOLD = '\033[1m'


def xor_blocos(blocos: list[bytes]) -> bytes:
    if not blocos:
        raise ValueError("É necessário pelo menos um bloco para XOR.")

    resultado = bytearray(len(blocos[0]))
    for bloco in blocos:
        for i, byte_val in enumerate(bloco):
            resultado[i] ^= byte_val

    return bytes(resultado)


def crc8(dados: bytes, polinomio: int = 0x07) -> int:
    crc = 0x00
    for byte in dados:
        crc ^= byte
        for _ in range(8):
            if crc & 0x80:
                crc = ((crc << 1) ^ polinomio) & 0xFF
            else:
                crc = (crc << 1) & 0xFF
    return crc


def scramble(payload: bytes, seed: int, offset: int = 0x5A5A) -> bytes:
    rng = random.Random(seed ^ offset)
    seq = bytes(rng.getrandbits(8) for _ in range(len(payload)))
    return xor_blocos([payload, seq])
