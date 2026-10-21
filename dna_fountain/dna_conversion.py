MAPA_BITS_PARA_DNA = {"00": "A", "01": "C", "10": "G", "11": "T"}
MAPA_DNA_PARA_BITS = {v: k for k, v in MAPA_BITS_PARA_DNA.items()}


def bytes_para_dna(dados: bytes) -> str:
    bits = "".join(f"{byte:08b}" for byte in dados)
    return "".join(MAPA_BITS_PARA_DNA[bits[i:i+2]] for i in range(0, len(bits), 2))


def dna_para_bytes(dna: str) -> bytes:
    bits = "".join(MAPA_DNA_PARA_BITS[base] for base in dna)
    return bytes(int(bits[i:i+8], 2) for i in range(0, len(bits), 8))


def exibir_conversao(dados: bytes):
    """Mostra passo a passo a conversão de cada byte para DNA."""
    print(f"\n{'='*50}")
    print(f"  Conversão Bytes → DNA")
    print(f"{'='*50}")
    print(f"  {'Byte':>6}  {'Hex':>5}  {'Bits':>10}  {'DNA':>6}")
    print(f"  {'----':>6}  {'---':>5}  {'--------':>10}  {'----':>6}")
    for byte_val in dados:
        bits = f"{byte_val:08b}"
        bases = "".join(MAPA_BITS_PARA_DNA[bits[i:i+2]] for i in range(0, 8, 2))
        print(f"  {byte_val:>6d}  0x{byte_val:02X}   {bits:>10}  {bases:>6}")
    print()
