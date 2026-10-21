import random
from .config import DNAFountainConfig

class MessageProcessor:
    def __init__(self, config: DNAFountainConfig):
        self.config = config

    def string_para_bytes(self, mensagem: str) -> bytes:
        return mensagem.encode("utf-8")

    def segmentar(self, dados: bytes) -> list[bytes]:
        blocos = []
        tamanho = self.config.tamanho_bloco_l
        
        for i in range(0, len(dados), tamanho):
            bloco = dados[i : i + tamanho]
            faltam = tamanho - len(bloco)
            
            if faltam > 0:
                if self.config.padding_aleatorio:
                    rng = random.Random(self.config.padding_seed)
                    padding = bytes(rng.getrandbits(8) for _ in range(faltam))
                    bloco += padding
                else:
                    bloco = bloco.ljust(tamanho, b"\x00")
            blocos.append(bloco)
            
        return blocos

    def preprocessar(self, mensagem: str) -> tuple[bytes, list[bytes], int]:
        dados = self.string_para_bytes(mensagem)
        blocos = self.segmentar(dados)
        return dados, blocos, len(blocos)

    def reconstruir(self, descobertos: dict, K: int, tamanho_original: int) -> bytes | None:
        if len(descobertos) != K:
            return None
        dados = b"".join(descobertos[i] for i in range(K))
        return dados[:tamanho_original]
