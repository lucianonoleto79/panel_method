---
title: "Gabarito — Trabalho de CFD: cilindro em Re = 100"
subtitle: "Resultados de referência obtidos com os três códigos (uso do professor)"
lang: pt-BR
---

Resultados das rodadas feitas em 27/11/2026 com as versões dos repositórios em que este
gabarito foi versionado. Tabela completa em `results/06_trabalho_cilindro/saida.txt`, gerada por
`python trabalho_cilindro/comparar.py`. Os valores da literatura foram citados de memória:
confira nas fontes antes de usá-los na correção.

# Rodadas

| Solução | Malha | Δt | Passos | Tempo medido |
|---|---|---|---|---|
| Cordel (FEM, CBS), `cylinder --tmax 200` | 16 043 nós, 31 852 triângulos | 5,3 × 10⁻⁴ | 377 012 | 7,8 h (dividindo a máquina com as rodadas do fvm-cfd; sozinho, ~45 ms/passo ≈ 4,7 h) |
| fvm-cfd, malha O, `cylinder.py --NT 128` (t = 150) | 128 × 128 = 16 384 volumes | 5,08 × 10⁻³ | 29 525 | 6 min sozinho; 17 min dividindo a máquina |
| fvm-cfd, malha uniforme, `cylinder2d.py --h 0.04` | 550 × 350 = 192 500 células | 0,01 | 20 001 | 1,8 h |
| painéis, N = 64 | 64 painéis | — | — | < 1 s |

# Tabela-resumo (Re = 100, médias na metade final)

| | Cordel | fvm-cfd, malha O | fvm-cfd, uniforme + FI | painéis | literatura |
|---|---|---|---|---|---|
| C<sub>d</sub> médio | 1,387 | 1,338 | 1,451 | 0 (6 × 10⁻¹⁶) | 1,33–1,35 |
| C<sub>l</sub>′ RMS | 0,243 | 0,228 | 0,258 | — | |
| C<sub>l</sub>′ amplitude | 0,344 | 0,323 | 0,366 | — | 0,32–0,34 |
| C<sub>d</sub>′ RMS | 0,0071 | 0,0063 | 0,0078 | — | |
| St (FFT de C<sub>l</sub>) | 0,170 | 0,163 | 0,175 | — | 0,164–0,165 |
| C<sub>p</sub> de estagnação | 1,135 | 1,061 | 1,218 | 0,990 (exato: 1) | |
| C<sub>p</sub> mínimo (ângulo da estagnação) | −1,168 (86°) | −1,129 (83°) | −1,226 (84°) | −2,990 (87°; exato −3 em 90°) | |
| C<sub>p</sub> de base | −0,696 | −0,705 | −0,728 | 0,990 | |
| C<sub>d,p</sub> (integral de C<sub>p</sub>) | 1,038 (75 % do C<sub>d</sub>) | 0,987 (74 %) | 1,093 (75 %) | 0 | |

O St da malha O pelos cruzamentos de zero de C<sub>l</sub> (11 ciclos, `resultados.json`) é 0,1637;
pela FFT, 0,1633. A diferença ilustra a resolução em frequência da FFT numa janela de 75 D/U.

# Respostas esperadas

**A1.**

- C<sub>p</sub> = 1 − 4 sen²θ.
- O arrasto é nulo porque C<sub>p</sub>(θ) é simétrico em relação ao eixo vertical (frente e trás iguais).
- No polígono regular de N painéis, os pontos de controle ficam no círculo inscrito
  ligeiramente menor, mas a solução de fontes reproduz o escoamento em torno do polígono,
  e nos pontos médios dos lados ela coincide com a do círculo até o arredondamento. O teste
  `tests/test_validacao.py` verifica isso.
- A hipótese responsável é a ausência de viscosidade: sem camada limite, não há
  separação nem esteira.
- A solução potencial continua útil só perto do ponto de estagnação. Na figura
  `comparacao_cp.png`, as curvas se afastam visivelmente a partir de ~20–25°; o C<sub>p</sub>
  potencial cruza zero em 30°, e o viscoso, em ~40°.

**A3.**

- O transitório dura até t ≈ 50 no Cordel e na malha uniforme, que partem de uma
  perturbação pequena (a instabilidade cresce exponencialmente). Na malha O ele dura até
  t ≈ 20, porque a rotação inicial do cilindro dispara a emissão logo.
- A janela t ≥ 100 (Cordel e malha uniforme) ou t ≥ 75 (malha O) contém 10 a 17 ciclos.
  Uma boa justificativa mostra que a média de C<sub>d</sub> muda menos de 0,1 % ao mudar o início da
  janela.
- C<sub>d</sub> oscila com 2f porque o arrasto é o mesmo quando o vórtice sai de cima ou de baixo,
  enquanto C<sub>l</sub> troca de sinal. O espectro de C<sub>d</sub> tem pico em 2 St ≈ 0,33.
- A malha O cai dentro da faixa da literatura em C<sub>d</sub> e St. O Cordel fica 3 % acima em C<sub>d</sub> e
  St. A malha uniforme fica 8 % acima em C<sub>d</sub> e 6 % em St (causas em A7).

**A4.**

- Concordância com a solução potencial só até ~20°.
- C<sub>p</sub> de estagnação maior que 1 tem duas causas, e o aluno deve separá-las:
    - **Efeito viscoso real**: a malha O, referenciada a p∞, dá 1,06; em Re finito a
      desaceleração viscosa na estagnação soma tensão normal e pressão.
    - **Referência de pressão**: Cordel (p = 0 na saída, x = 25) e malha uniforme (p = 0 em
      x = 16) medem a pressão em relação à saída, que está na esteira com p < p∞. As curvas
      ficam deslocadas para cima por uma constante: +0,07 no Cordel e +0,16 na malha
      uniforme no ponto de estagnação.
- Uma constante não altera C<sub>d,p</sub> (∫ cos θ dθ = 0). Por isso o C<sub>d,p</sub> continua comparável.
- C<sub>p</sub> mínimo ≈ −1,13 a −1,23, entre 83° e 86° da estagnação, contra −3 a 90° no potencial.
- Base em C<sub>p</sub> ≈ −0,70 nas três soluções, quase constante de ~110° até 180°. Esse patamar é a
  região separada; a separação média em Re = 100 fica por volta de 115–120° da estagnação
  (confira na literatura).
- C<sub>d,p</sub> ≈ 74–75 % do C<sub>d</sub>: em Re = 100 a pressão domina; o atrito responde por ~25 %.

**A5.**

- **Malha O** (R<sub>ext</sub> = 40 D), C<sub>d</sub> = 1,3385 / 1,3353 / 1,3336 e St = 0,1637 / 0,1641 /
  0,1643 em NT = 128 / 192 / 256.
    - St converge com ordem observada ≈ 1,9 (extrapolado: 0,1645).
    - C<sub>d</sub> converge devagar (ordem aparente ≈ 0,7), provavelmente limitado pelo cálculo das
      forças na parede e pelo Δt, que muda com a malha.
- **Malha uniforme + FI**, C<sub>d</sub> = 1,456 / 1,398 / 1,451 e St = 0,1725 / 0,1762 / 0,1749 em
  D/12,5 / D/20 / D/25: **não monótona**.
    - A fronteira imersa de fração volumétrica é de 1ª ordem na parede.
    - Cada malha corta o círculo de um jeito: o "cilindro efetivo" muda de malha para
      malha, e a extrapolação de Richardson não se aplica.
- **Custo:** a malha O chega a C<sub>d</sub> ≈ 1,34 com 16 mil volumes em 6 min. A malha uniforme não
  chega lá nem com 190 mil células.
- **Resposta esperada:** a ordem de 2ª do esquema só aparece quando a parede também é de
  2ª ordem e a malha resolve a camada limite. A camada limite tem espessura ~ D/√Re =
  0,1 D: são 2–3 células na malha uniforme D/25, contra a 1ª célula a 0,017 D na malha O.

**A6.**

- **Cordel.** A menor aresta é ~0,005 D (1ª camada) e Δt = 5,3 × 10⁻⁴.
    - CFL ≈ U Δt / h ≈ 0,1.
    - ν Δt / h² ≈ 0,21: é o limite difusivo que manda.
    - O CBS explícito usa o Δt local de estabilidade dos menores elementos.
- **Malha O.** Δr na parede é 0,017 D e r Δθ = 0,025 D; Δt = 5,1 × 10⁻³.
    - CFL medido de 0,24.
    - A difusão é implícita (Crank–Nicolson), sem limite: ν Δt / Δr² ≈ 0,17 já
      complicaria um esquema explícito.
- **Malha uniforme.** h = 0,04 e Δt = 0,01.
    - CFL ≈ 0,4 (|u|<sub>máx</sub> ≈ 1,6).
    - ν Δt / h² = 0,0625 ≤ 1/8.
    - Em Re = 40 o limite difusivo passa a mandar (Δt = 0,0064). A rodada com 0,1 h²/ν
      → 0,2 h²/ν divergiu; esse erro foi corrigido e serve de exemplo.
- **Por que o Δt do Cordel é 10–20× menor:** os elementos da parede são 3 a 8× menores, e o
  limite difusivo explícito escala com h². Custo total: 377 mil passos × ~45 ms contra
  29 mil × ~12 ms.

**A7.** Fontes de diferença, em ordem aproximada de importância:

- **Domínio e bloqueio.**
    - O estudo de domínio do fvm-cfd (`results/estudo_dominio_cilindro_Re100/`) mostra
      que o retângulo [−10, 20] × [±10] D (bloqueio 5 %, laterais com deslizamento) dá
      C<sub>d</sub> = 1,354 e **St = 0,1706**, praticamente o St do Cordel (0,170), que usa um domínio
      parecido, com u = 1 imposto no topo e no fundo.
    - A malha uniforme tem bloqueio de 7 %.
    - Esse é o efeito dominante em St.
- **Representação da parede e resolução da camada limite.** Esse é o efeito dominante no
  C<sub>d</sub> da malha uniforme: em domínio igual, o retangular com fronteira imersa ainda difere
  2–3 % da malha O (README do fvm-cfd).
- **Referência de pressão.** Afeta o C<sub>p</sub> (deslocamento constante), não as forças.
- **Janela de média e perturbação inicial.** < 0,5 %.
- **Classificação:**
    - painéis × solução exata e estudos de malha = **verificação**;
    - comparação com St e C<sub>d</sub> da literatura experimental e numérica = **validação**,
      lembrando que os valores numéricos da literatura também têm erro de malha e domínio.
- **Mais confiável:** a malha O. Tem convergência de malha demonstrada, domínio grande
  (40 D) com condição de entrada potencial e parede ajustada ao corpo.

**B1 (Re = 40).** Malha uniforme D/25, t = 100:

- C<sub>d</sub> = 1,695, contra 1,5–1,6 na literatura: o bloqueio de 7 % e a fronteira imersa
  superestimam, como em Re = 100.
- C<sub>l</sub>′ ≈ 3 × 10⁻⁴, então St não é definido. Os três códigos agora devolvem nan (—) quando a
  amplitude de C<sub>l</sub> < 0,02. Antes, a FFT do resíduo do transitório dava um falso St = 0,1275,
  armadilha boa para discutir.
- L<sub>w</sub>/D = 2,24, contra 2,2–2,3 na literatura.
- C<sub>p</sub>: 1,351 na estagnação (com o deslocamento da referência), mínimo −1,02 a 88°, base
  −0,51.

**B3.** Malha O, R<sub>ext</sub> = 10 / 20 / 40 / 80 D:

- C<sub>d</sub> = 1,374 / 1,351 / 1,339 / 1,334;
- St = 0,1660 / 0,1646 / 0,1637 / 0,1634 (README do fvm-cfd).
