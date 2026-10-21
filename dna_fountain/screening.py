from .config import DNAFountainConfig
from .utils import C

class BiologicalScreening:
    def __init__(self, config: DNAFountainConfig):
        self.gc_min = config.gc_min
        self.gc_max = config.gc_max
        self.max_homopolimero = config.max_homopolimero

    @staticmethod
    def calcular_gc(dna: str) -> float:
        if not dna: return 0.0
        return (dna.count("G") + dna.count("C")) / len(dna)

    @staticmethod
    def maior_homopolimero(dna: str) -> int:
        if len(dna) <= 1: return len(dna)
        maior, atual = 1, 1
        for i in range(1, len(dna)):
            if dna[i] == dna[i - 1]:
                atual += 1
                if atual > maior: maior = atual
            else:
                atual = 1
        return maior

    def verificar(self, dna: str) -> bool:
        gc = self.calcular_gc(dna)
        homo = self.maior_homopolimero(dna)
        return self.gc_min <= gc <= self.gc_max and homo <= self.max_homopolimero

    def exibir_resultado(self, dna: str):
        gc = self.calcular_gc(dna)
        homo = self.maior_homopolimero(dna)
        gc_ok = self.gc_min <= gc <= self.gc_max
        homo_ok = homo <= self.max_homopolimero
        
        status = f"{C.GREEN}✓ APROVADO{C.RESET}" if (gc_ok and homo_ok) else f"{C.RED}✗ REJEITADO{C.RESET}"
        cor_gc = C.GREEN if gc_ok else C.RED
        cor_homo = C.GREEN if homo_ok else C.RED

        print(f"  DNA: {dna[:40]}{'...' if len(dna) > 40 else ''}")
        print(f"    Comprimento: {len(dna)} bases")
        print(f"    GC: {cor_gc}{gc:.2%}{C.RESET} (faixa [{self.gc_min}, {self.gc_max}]) → {'OK' if gc_ok else 'FALHA'}")
        print(f"    Homopolímero: {cor_homo}{homo}{C.RESET} (máx {self.max_homopolimero}) → {'OK' if homo_ok else 'FALHA'}")
        print(f"    Resultado: {status}")
