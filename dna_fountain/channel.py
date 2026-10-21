import random
from .config import DNAFountainConfig
from .utils import C

class ErasureChannel:
    def __init__(self, config: DNAFountainConfig, rng: random.Random = None):
        self.prob_perda = config.prob_perda
        self.rng = rng or random.Random()

    def transmitir(self, oligos: list[dict]) -> tuple[list[dict], list[dict]]:
        recebidos = []
        perdidos = []
        for oligo in oligos:
            if self.rng.random() >= self.prob_perda:
                recebidos.append(oligo)
            else:
                perdidos.append(oligo)
        return recebidos, perdidos

    def exibir_resultado(self, total: int, recebidos: int, perdidos: int):
        taxa_real = perdidos / total if total > 0 else 0.0
        print(f"\n{C.CYAN}{'='*50}")
        print(f"  Canal de Apagamento")
        print(f"{'='*50}{C.RESET}")
        print(f"  Prob. de perda configurada: {self.prob_perda:.2%}")
        print(f"  Oligos enviados:           {total}")
        print(f"  Oligos recebidos:          {C.GREEN}{recebidos} ({recebidos/total:.2%}){C.RESET}")
        print(f"  Oligos perdidos:           {C.RED}{perdidos} ({taxa_real:.2%}){C.RESET}\n")
