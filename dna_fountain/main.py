import sys
import os
import random

if sys.stdout.encoding != "utf-8":
    sys.stdout.reconfigure(encoding="utf-8")

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from dna_fountain.config import DNAFountainConfig
from dna_fountain.preprocessing import MessageProcessor
from dna_fountain.soliton import RobustSoliton
from dna_fountain.screening import BiologicalScreening
from dna_fountain.lt_codec import DNAFountainEncoder, DNAFountainDecoder
from dna_fountain.channel import ErasureChannel
from dna_fountain.dna_conversion import exibir_conversao
from dna_fountain.utils import C
import json

def gerar_grafico_ripple(historico: list[int], arquivo: str = "grafico_ripple.html"):
    if not historico: return
    html = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Gr&aacute;fico de Ripple</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>body {{ font-family: sans-serif; padding: 20px; text-align: center; }}</style>
</head>
<body>
    <h2>Evolu&ccedil;&atilde;o do Tamanho do Ripple (Peeling Decoder)</h2>
    <div style="width: 100%; max-width: 800px; margin: 0 auto;">
        <canvas id="chartRipple"></canvas>
    </div>
    <script>
        new Chart(document.getElementById('chartRipple').getContext('2d'), {{
            type: 'line',
            data: {{
                labels: {json.dumps(list(range(1, len(historico)+1)))},
                datasets: [{{
                    label: 'Tamanho do Ripple',
                    data: {json.dumps(historico)},
                    borderColor: 'rgb(75, 192, 192)',
                    backgroundColor: 'rgba(75, 192, 192, 0.2)',
                    fill: true, tension: 0.1
                }}]
            }},
            options: {{ scales: {{ y: {{ min: 0 }} }} }}
        }});
    </script>
</body>
</html>"""
    with open(arquivo, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"\n  {C.MAGENTA}📊 Gráfico da evolução do ripple gerado em: {arquivo}{C.RESET}")

def separador(titulo: str):
    print(f"\n{C.MAGENTA}{C.BOLD}{'#'*70}")
    print(f"#  {titulo}")
    print(f"{'#'*70}{C.RESET}")

def exibir_configuracao(config: DNAFountainConfig):
    separador("CONFIGURAÇÃO DO SISTEMA (config.py)")
    print(f"""
  Pré-processamento:
    Tamanho do bloco (L):       {config.tamanho_bloco_l} bytes
    Padding aleatório:          {'SIM' if config.padding_aleatorio else 'NÃO (zeros)'}

  Distribuição Robusta de Sóliton:
    Constante c:                {config.soliton_c}
    Probabilidade de falha δ:   {config.soliton_delta}

  Codificação LT:
    Número de droplets (N):     {config.num_droplets}
    Bits da seed:               {config.bits_seed}
    Máscara XOR da seed:        0x{config.seed_xor_mask:04X}
    Scrambler:                  {'ATIVO' if config.scrambler_ativo else 'DESATIVADO'}
    CRC-8:                      {'ATIVO' if config.crc_ativo else 'DESATIVADO'}

  Screening:
    GC mínimo:                  {config.gc_min:.0%}
    GC máximo:                  {config.gc_max:.0%}
    Max homopolímero:           {config.max_homopolimero}

  Canal de apagamento:
    Probabilidade de perda:     {config.prob_perda:.0%}

  Reprodutibilidade:
    Seed do experimento:        {config.seed_experimento}
""")

def pipeline(mensagem: str, config: DNAFountainConfig):
    processor = MessageProcessor(config)
    screening = BiologicalScreening(config)
    canal = ErasureChannel(config, rng=random.Random(config.seed_experimento))
    encoder = DNAFountainEncoder(config)
    decoder = DNAFountainDecoder(config)

    # 1. Pré-processamento
    separador("ETAPA 1 — PRÉ-PROCESSAMENTO")
    dados_originais, blocos, K = processor.preprocessar(mensagem)
    
    margem_peeling = 0.80
    fator_canal = 1.0 - config.prob_perda
    num_droplets_necessarios = int((K * (1 + margem_peeling)) / fator_canal)
    if config.num_droplets < num_droplets_necessarios:
        print(f"  {C.YELLOW}Ajustando droplets de {config.num_droplets} para {num_droplets_necessarios} (K={K}, perda={config.prob_perda:.0%}){C.RESET}")
        config.num_droplets = num_droplets_necessarios
    
    print(f"\n  Mensagem original:  \"{mensagem}\"")
    print(f"  Tamanho em bytes:   {len(dados_originais)}")
    print(f"  Tamanho do bloco L: {config.tamanho_bloco_l} bytes")
    print(f"  Número de blocos K: {K}")
    print(f"\n  Blocos segmentados:")
    for i, bloco in enumerate(blocos):
        bits = " ".join(f"{b:08b}" for b in bloco)
        texto_seguro = "".join(chr(b) if 32 <= b <= 126 else '.' for b in bloco)
        print(f"    Bloco {i:>2}: {bloco.hex()}  |  {bits}  |  \"{texto_seguro}\"")

    # 2. Distribuição de Sóliton
    separador("ETAPA 2 — DISTRIBUIÇÃO ROBUSTA DE SÓLITON")
    soliton = RobustSoliton(config, K)
    soliton.exibir()

    # 3. Codificação LT
    separador("ETAPA 3 — CODIFICAÇÃO LT + SCREENING")
    oligos, rejeitados, log_rejeicoes = encoder.codificar(blocos, soliton)
    
    total_tentativas = len(oligos) + rejeitados
    taxa_rejeicao = rejeitados / total_tentativas if total_tentativas > 0 else 0

    print(f"  Droplets válidos gerados:   {len(oligos)}")
    print(f"  Droplets rejeitados:        {rejeitados} ({taxa_rejeicao:.1%})")
    print(f"  Seeds consumidas:           {total_tentativas}")

    print(f"\n  Primeiros 5 droplets aceitos:")
    for o in oligos[:5]:
        print(f"    Seed: {o['seed']:>4} | Índices: {o['indices']} | DNA: {o['dna'][:30]}...")

    if oligos:
        print(f"\n  Exemplo de screening detalhado (seed={oligos[0]['seed']}):")
        screening.exibir_resultado(oligos[0]["dna"])

    # 4. Canal
    separador("ETAPA 4 — CANAL DE APAGAMENTO")
    recebidos, perdidos = canal.transmitir(oligos)
    canal.exibir_resultado(len(oligos), len(recebidos), len(perdidos))

    # 5. Decodificação
    separador("ETAPA 5 — DECODIFICAÇÃO LT (PEELING)")
    dnas_recebidos = [o["dna"] for o in recebidos]
    descobertos, metricas = decoder.decodificar(dnas_recebidos, soliton)

    print(f"\n  Oligos recebidos:          {C.GREEN}{len(dnas_recebidos)}{C.RESET}")
    if config.crc_ativo: print(f"  CRC falhas (descartados):  {C.RED}{metricas['crc_falhas']}{C.RESET}")
    print(f"  Blocos recuperados:        {C.BLUE}{metricas['blocos_recuperados']}/{K}{C.RESET}")
    
    if metricas['decodificacao_completa']:
        print(f"  Decodificação completa:    {C.GREEN}SIM ✓{C.RESET}")
    else:
        print(f"  Decodificação completa:    {C.RED}NÃO ✗{C.RESET}")
        
    print(f"  Tamanho máx do ripple:     {metricas['ripple_max']}")

    # 6. Reconstrução
    separador("ETAPA 6 — RECONSTRUÇÃO DA MENSAGEM")
    mensagem_recuperada = processor.reconstruir(descobertos, K, len(dados_originais))

    if mensagem_recuperada:
        texto = mensagem_recuperada.decode("utf-8", errors="replace")
        print(f"\n  Mensagem original:    \"{mensagem}\"")
        print(f"  Mensagem recuperada:  \"{texto}\"")
        if mensagem_recuperada == dados_originais:
            print(f"  Resultado:            {C.GREEN}PERFEITA ✓{C.RESET}")
        else:
            print(f"  Resultado:            {C.RED}CORROMPIDA ✗{C.RESET}")
    else:
        print(f"\n  {C.RED}✗ FALHA na reconstrução! Recuperados {metricas['blocos_recuperados']}/{K}{C.RESET}")

    # Gerar gráfico ao final
    if metricas['ripple_historico']:
        gerar_grafico_ripple(metricas['ripple_historico'])

def main():
    print(f"{C.CYAN}{C.BOLD}  ╔══════════════════════════════════════════════════╗")
    print("  ║         DNA FOUNTAIN — Simulador Completo       ║")
    print(f"  ╚══════════════════════════════════════════════════╝{C.RESET}")

    config = DNAFountainConfig()
    exibir_configuracao(config)

    print("─" * 70)
    entrada = input("  Digite a mensagem (ou caminho para um arquivo .txt): ").strip()
    print("─" * 70)

    if not entrada:
        mensagem = "DNA Fountain - TCC"
    elif os.path.isfile(entrada) and entrada.endswith(".txt"):
        with open(entrada, "r", encoding="utf-8") as f:
            mensagem = f.read()
        print(f"  {C.CYAN}📥 Arquivo '{entrada}' carregado com sucesso! ({len(mensagem)} caracteres){C.RESET}")
    else:
        mensagem = entrada

    pipeline(mensagem, config)

if __name__ == "__main__":
    main()
