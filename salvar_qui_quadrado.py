import json
import math
from dna_fountain.config import DNAFountainConfig
from dna_fountain.soliton import RobustSoliton

def salvar():
    config = DNAFountainConfig()
    config.soliton_c = 0.1
    config.soliton_delta = 0.10
    K = 32
    soliton = RobustSoliton(config, K)
    
    classes_labels = [
        "d=1", "d=2", "d=3", "d=4", "d=5", "d=6", "d=7", "d=8", "d=9", "d=10", "d=11", "d=12",
        "d in {13,14}", "d in {15,16}", "d in {17,18}", "d in {19-21}", "d in {22-25}", "d in {26-32}"
    ]
    classes = [
        [1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12],
        [13, 14], [15, 16], [17, 18], [19, 20, 21], [22, 23, 24, 25], list(range(26, 33))
    ]
    
    obs = [115, 412, 162, 74, 50, 32, 20, 23, 244, 7, 6, 8, 7, 9, 3, 6, 3, 14]
    n_aceitos = sum(obs)
    esp = [round(n_aceitos * sum(soliton.mu[d] for d in c), 2) for c in classes]
    
    chi2 = sum((o - e)**2 / e for o, e in zip(obs, esp))
    
    dados = {
        "ensaio": "Screening Biológico - Aderência Qui-Quadrado",
        "parametros": {"K": 32, "c": 0.1, "delta": 0.10, "candidatos": 2000, "seed_ensaio": 41},
        "aceitos": n_aceitos,
        "rejeitados": 2000 - n_aceitos,
        "graus_de_liberdade": len(classes) - 1,
        "chi2_calculado": round(chi2, 4),
        "tabela_classes": [
            {"classe": lbl, "observado": o, "esperado": e}
            for lbl, o, e in zip(classes_labels, obs, esp)
        ]
    }
    
    with open("screening_qui_quadrado.json", "w", encoding="utf-8") as f:
        json.dump(dados, f, indent=2, ensure_ascii=False)
    print("Salvo com sucesso em screening_qui_quadrado.json")

if __name__ == "__main__":
    salvar()
