import math
import random
from .config import DNAFountainConfig
from .utils import C

class RobustSoliton:
    def __init__(self, config: DNAFountainConfig, K: int):
        self.K = K
        self.c = config.soliton_c
        self.delta = config.soliton_delta
        
        self.rho = self._soliton_ideal()
        self.tau, self.R, self.m = self._funcao_tau()
        
        self.Z = sum(self.rho[d] + self.tau[d] for d in range(1, self.K + 1))
        self.mu = [(self.rho[d] + self.tau[d]) / self.Z for d in range(self.K + 1)]
        self.grau_medio = sum(d * self.mu[d] for d in range(1, self.K + 1))
        self.cdf = self._calcular_cdf()

    def _soliton_ideal(self) -> list[float]:
        rho = [0.0] * (self.K + 1)
        rho[1] = 1.0 / self.K
        for d in range(2, self.K + 1):
            rho[d] = 1.0 / (d * (d - 1))
        return rho

    def _funcao_tau(self) -> tuple[list[float], float, int]:
        R = self.c * math.log(self.K / self.delta) * math.sqrt(self.K)
        m = max(1, min(self.K, int(self.K / R)))

        tau = [0.0] * (self.K + 1)
        for d in range(1, m):
            tau[d] = R / (d * self.K)
        tau[m] = max(0.0, R * math.log(R / self.delta) / self.K)
        return tau, R, m

    def _calcular_cdf(self) -> list[float]:
        acumulada = []
        soma = 0.0
        for d in range(1, len(self.mu)):
            soma += self.mu[d]
            acumulada.append(soma)
        return acumulada

    def sortear_grau_e_indices(self, seed: int) -> tuple[int, list[int]]:
        rng = random.Random(seed)
        u = rng.random()
        grau = next((d + 1 for d, p in enumerate(self.cdf) if u <= p), self.K)
        indices = sorted(rng.sample(range(self.K), grau))
        return grau, indices

    def exibir(self):
        print(f"\n{C.CYAN}{'='*60}")
        print(f"  Distribuição Robusta de Sóliton  —  K={self.K}, c={self.c}, δ={self.delta}")
        print(f"  R = {self.R:.4f}  |  m = {self.m}")
        print(f"  Z = {self.Z:.6f}  |  Grau médio = {self.grau_medio:.4f}")
        print(f"{'='*60}{C.RESET}")
        print(f"  {'d':>3}  {'ρ(d)':>10}  {'τ(d)':>10}  {'μ(d)':>10}  {'μ(d) %':>8}")
        print(f"  {'---':>3}  {'----------':>10}  {'----------':>10}  {'----------':>10}  {'------':>8}")
        for d in range(1, self.K + 1):
            print(f"  {d:>3}  {self.rho[d]:>10.6f}  {self.tau[d]:>10.6f}  {C.YELLOW}{self.mu[d]:>10.6f}{C.RESET}  {self.mu[d]*100:>7.2f}%")
        print(f"  {'':>3}  {'':>10}  {'':>10}  {'Σ =':>10}  {sum(self.mu):>7.4f}\n")
