---
title: "Trabalho de CFD — Escoamento em torno de um cilindro e a esteira de von Kármán"
subtitle: "Três códigos, um mesmo problema: potencial, elementos finitos e volumes finitos"
lang: pt-BR
---

| | |
|---|---|
| **Disciplina** | [a definir pelo professor] |
| **Professor(a)** | [a definir] |
| **Grupos** | [a definir: número de alunos por grupo] |
| **Entrega** | [a definir: data e forma de entrega] |

# 1. Contexto e objetivos

O escoamento em torno de um cilindro circular é o problema clássico de corpo rombudo. Em
número de Reynolds baixo o escoamento é permanente e simétrico. A partir de Re ≈ 47 a
esteira fica instável e passa a soltar vórtices alternados, um de cada lado: é a
**esteira de von Kármán**. Esse desprendimento produz uma força lateral oscilante. É o que
faz cabos "cantarem" no vento e pode excitar a vibração de chaminés, pilares e
tubulações.

Neste trabalho você vai estudar o cilindro em **Re = 100** com três códigos que partem
de hipóteses diferentes:

- um **método dos painéis** (escoamento potencial: sem viscosidade);
- o **Cordel**, que resolve Navier–Stokes incompressível por **elementos finitos**;
- o **fvm-cfd**, que resolve as mesmas equações por **volumes finitos**, de duas formas: numa
  **malha O** ajustada ao corpo e numa **malha cartesiana uniforme** em que o cilindro entra
  por **fronteira imersa**.

Ao final do trabalho você deve ser capaz de:

1. escrever as equações que cada código resolve e explicar o que cada hipótese elimina da
   física (por exemplo: por que o escoamento potencial dá arrasto nulo);
2. extrair de uma simulação transiente as grandezas de engenharia: Cd médio, amplitude
   de Cl, número de Strouhal e distribuição de Cp. Você deve justificar a janela de tempo
   usada nas médias;
3. explicar, com números, de onde vêm as diferenças entre as três soluções viscosas
   (malha, representação da parede, domínio, condições de contorno, passo de tempo);
4. distinguir **verificação** (o código resolve bem as equações?) de **validação** (as
   equações e o modelo representam o experimento?).

# 2. O problema

| Item | Valor |
|---|---|
| Geometria | cilindro circular de diâmetro D, em escoamento uniforme U |
| Adimensionalização | D = 1, U = 1, ρ = 1; tempo em unidades de D/U |
| Número de Reynolds | Re = U D / ν = **100** (ν = 0,01); opcional: Re = 40 |
| Escoamento | incompressível, 2D, laminar, partindo do repouso relativo (partida impulsiva) |
| Perturbação inicial | para disparar a assimetria: v = 0,05 exp(−[(x−2)² + y²]) no Cordel e na malha uniforme do fvm-cfd; na malha O, o cilindro gira brevemente no início (u<sub>θ</sub> = 0,5 sen(πt/2) na parede, t < 2) |
| Duração | t = 200 D/U (Cordel, malha uniforme) e t = 150 (malha O); médias na metade final |

Grandezas de interesse (a força F é por unidade de comprimento do cilindro):

- C<sub>d</sub> = F<sub>x</sub> / (½ ρ U² D) e C<sub>l</sub> = F<sub>y</sub> / (½ ρ U² D);
- número de Strouhal St = f D / U, em que f é a frequência dominante de C<sub>l</sub>;
- C<sub>p</sub> = (p − p<sub>∞</sub>) / (½ ρ U²) na parede, em função do ângulo θ medido a partir do ponto
  de estagnação.

> **Para pensar antes de começar.** Que velocidade de água (ν ≈ 10⁻⁶ m²/s) em torno de um
> fio de 1 mm de diâmetro dá Re = 100? E de ar (ν ≈ 1,5 × 10⁻⁵ m²/s)? Qual seria a
> frequência de desprendimento em hertz em cada caso?

# 3. Os três códigos

| Código | Equações | Discretização | Repositório |
|---|---|---|---|
| **panel_method** | escoamento potencial: ∇²φ = 0, com V·n = 0 na parede | painéis retos de fontes de intensidade constante; sistema linear denso N × N | `github.com/lucianonoleto79/panel_method` |
| **Cordel** | Navier–Stokes incompressível 2D | elementos finitos lineares (P1) em triângulos; método CBS (*Characteristic-Based Split*), explícito, com Δt local de estabilidade | `github.com/lucianonoleto79/cordel` |
| **fvm-cfd** | Navier–Stokes incompressível 2D | volumes finitos em malha desencontrada, passo fracionado (projeção incremental), Adams–Bashforth 2 e QUICK na convecção. Duas abordagens: `cylinder.py`, **malha O** polar ajustada ao corpo, com Crank–Nicolson na difusão; `cylinder2d.py`, **malha cartesiana uniforme**, difusão explícita e cilindro por **fronteira imersa** (fração de volume) | `github.com/lucianonoleto79/fvm-cfd` |

Diferenças de configuração que importam para a análise:

| | Cordel | fvm-cfd, malha O (`cylinder.py`) | fvm-cfd, malha uniforme (`cylinder2d.py`) |
|---|---|---|---|
| Domínio (centro do cilindro na origem) | retângulo [−10, 25] × [−10, 10] | círculo de raio R<sub>ext</sub> = 40 D | retângulo [−6, 16] × [−7, 7] |
| Malha | não estruturada; anel de 12 camadas junto à parede (1ª camada 0,005 D) e zonas de refino; 16 043 nós | polar 128 × 128 (r × θ), r em progressão geométrica, células quase quadradas; 16 384 volumes | cartesiana uniforme h = D/25 (550 × 350 = 192 500 células); a parede corta as células |
| Parede | ajustada ao corpo (*body-fitted*), u = v = 0 nos nós | ajustada ao corpo: a parede é a linha r = D/2 da malha | fronteira imersa: a velocidade é multiplicada por (1 − α), com α a fração sólida da célula |
| Contorno externo | topo e fundo com u = 1, v = 0; saída com p = 0 | metade a montante: escoamento potencial em torno do cilindro; metade a jusante: saída convectiva | topo e fundo com deslizamento (v = 0, ∂u/∂y = 0); saída convectiva, p = 0 |
| Passo de tempo | Δt ≈ 5 × 10⁻⁴, limitado pelos elementos finos da parede | Δt ≈ 0,005, por CFL (difusão implícita) | Δt = 0,01, por CFL (difusão explícita: ν Δt/h² ≤ 1/8) |
| Forças | integral de −p n + τ·n nas arestas da parede | integral de −p n + τ·n na parede (p extrapolada; τ por parábola de 3 pontos) | soma da força da fronteira imersa sobre as células (F = −Σ f ΔV) |

O método dos painéis não tem tempo, viscosidade nem esteira. Entra no trabalho como
**referência invíscida**: a solução exata do cilindro em escoamento potencial é
C<sub>p</sub> = 1 − 4 sen²θ.

# 4. Material fornecido

- os três repositórios, com instruções de instalação nos `README.md` (Python ≥ 3.10; numpy,
  scipy, matplotlib);
- o caso do cilindro em cada código:
    - `cordel/casos/cylinder.py`, com a malha pronta em `casos/dados/cilindro_Re100.msh`;
    - `fvm-cfd/python/cylinder.py` (malha O; o solver está em `python/polar.py`) e
      `fvm-cfd/python/cylinder2d.py` (malha uniforme com fronteira imersa);
    - `panel_method/examples/02_cilindro_paineis.py`;
- o script de comparação `panel_method/trabalho_cilindro/comparar.py`, que lê as saídas das
  três soluções viscosas, resolve o escoamento potencial e gera as figuras de comparação e
  uma tabela;
- [a definir pelo professor: o repositório `fvm-cfd` é privado; dar acesso aos grupos ou
  torná-lo público];
- [a definir pelo professor: arquivos de estado `.npz` das rodadas longas, se forem
  disponibilizados para que os grupos não precisem rodá-las desde o início].

# 5. Como rodar

Clone os três repositórios lado a lado, na mesma pasta:

```bash
git clone https://github.com/lucianonoleto79/cordel.git
git clone https://github.com/lucianonoleto79/fvm-cfd.git
git clone https://github.com/lucianonoleto79/panel_method.git
pip install numpy scipy matplotlib pytest
```

**Método dos painéis** (segundos):

```bash
cd panel_method
python examples/02_cilindro_paineis.py
python -m pytest tests -q
```

**Cordel** (FEM). O passo é pequeno, por isso a rodada é longa:

```bash
cd cordel
python cordel2008v4.py cylinder --itmax 200 --saida teste    # meça o tempo por passo antes
python cordel2008v4.py cylinder --tmax 200 --saida cylinder_trabalho
# se precisar interromper: crie resultados/cylinder_trabalho/cylinder_trabalho.pausar
# e retome com o mesmo comando + --continuar
```

**fvm-cfd** (FVM), as duas abordagens:

```bash
cd fvm-cfd
python python/cylinder.py --NT 128                                # malha O, t = 150
python python/cylinder.py --NT 192 256                            # malha O: estudo de malha
python python/cylinder2d.py --Re 100 --h 0.04 --tmax 200          # malha uniforme
python python/cylinder2d.py --Re 100 --h 0.08 0.05 --tmax 200     # malha uniforme: estudo de malha
```

Cada caso fica em `fvm-cfd/results/<caso>/`, com figuras, `resultados.json` e o C<sub>p</sub> médio
na parede em `cp_parede.csv`.

**Comparação** (depois das rodadas):

```bash
cd panel_method
python trabalho_cilindro/comparar.py      # figuras e tabela em results/06_trabalho_cilindro/
```

Tempos medidos num computador pessoal, para t = 200 (só ordem de grandeza):

| Rodada | Passos | Tempo |
|---|---|---|
| Cordel, malha padrão | ≈ 380 000 | 5 a 8 h |
| fvm-cfd, malha O 128 × 128 (t = 150) | ≈ 29 000 | ≈ 6 min |
| fvm-cfd, malha O 256 × 256 (t = 150) | ≈ 59 000 | ≈ 1 h |
| fvm-cfd, malha uniforme h = D/25 | 20 000 | ≈ 2 h |
| fvm-cfd, malha uniforme h = D/12,5 | 10 000 | ≈ 5 min |
| painéis, N = 64 | — | < 1 s |

**Para explorar, use malhas grossas e rodadas curtas.** Exemplos: `--NT 64` ou `--h 0.08` no fvm-cfd e
`--tmax 100` com `--inicio` a partir de um estado já desenvolvido no Cordel. Deixe as
rodadas finas para os resultados finais.

# 6. Tarefas

## Parte A — obrigatória (todos os grupos)

**A1. Escoamento potencial e o paradoxo de d'Alembert.**

1. Deduza C<sub>p</sub> = 1 − 4 sen²θ para o cilindro sem circulação.
2. Mostre, integrando C<sub>p</sub> na superfície, que o arrasto é nulo.
3. Rode o método dos painéis com N = 8, 16, 32, 64 e compare o C<sub>p</sub> nos pontos de controle
    com a solução exata. Explique por que o erro fica no nível do arredondamento para o
    polígono regular.
4. Qual hipótese do escoamento potencial é responsável pelo arrasto nulo? Onde, ao longo
    da superfície, você espera que essa solução continue útil num escoamento real?

**A2. Formulações.** Escreva as equações resolvidas pelo Cordel e pelo fvm-cfd (na malha O,
em coordenadas polares) e descreva o algoritmo de cada passo de tempo (CBS e projeção).
Explique:

- o papel da equação de Poisson para a pressão nos dois métodos;
- por que a malha desencontrada evita o desacoplamento "tabuleiro de xadrez" da pressão;
- o que a fronteira imersa faz de diferente de uma malha ajustada ao corpo;
- por que a malha O usa difusão implícita (Crank–Nicolson) e a malha uniforme pode usar
  difusão explícita.

**A3. Forças e Strouhal.** Rode o caso Re = 100 no Cordel e nas duas abordagens do fvm-cfd.

1. Faça os gráficos de C<sub>d</sub>(t) e C<sub>l</sub>(t).
2. Identifique o transitório (partida impulsiva, crescimento da instabilidade) e o regime
    periódico. Justifique a janela usada nas médias.
3. Calcule:
    - C<sub>d</sub> médio;
    - a amplitude e o RMS de C<sub>l</sub>;
    - a amplitude de C<sub>d</sub>;
    - St, pela FFT de C<sub>l</sub>.
4. Por que C<sub>d</sub> oscila com o dobro da frequência de C<sub>l</sub>?
5. Compare St e C<sub>d</sub> com a literatura: St ≈ 0,164–0,165 e C<sub>d</sub> ≈ 1,33–1,35 em Re = 100 (ver
   referências; confira os valores nas fontes).

**A4. Distribuição de pressão.** Compare o C<sub>p</sub> médio na parede das três soluções viscosas e
do método dos painéis num mesmo gráfico. Discuta:

- até que ângulo a solução potencial representa bem a pressão;
- o valor de C<sub>p</sub> no ponto de estagnação (por que pode ser um pouco maior que 1 em Re = 100?);
- o C<sub>p</sub> mínimo e sua posição;
- a pressão de base e sua relação com a separação.

Calcule o **arrasto de pressão** integrando o C<sub>p</sub> médio,
C<sub>d,p</sub> = −½ ∫<sub>0</sub><sup>2π</sup> C<sub>p</sub> cos θ dθ (θ a partir de +x). A diferença para o C<sub>d</sub> total é o
arrasto viscoso. Qual parcela domina em Re = 100?

**A5. Estudo de malha (fvm-cfd).**

1. Malha O: rode NT = 128, 192 e 256. Malha uniforme: rode h = D/12,5, D/20 e D/25 (ou mais
   fina, se o tempo permitir).
2. Faça gráficos de C<sub>d</sub> e St contra o número de volumes para as duas abordagens e estime a
   ordem de convergência observada. Se for possível, use a extrapolação de Richardson.
3. Em qual abordagem a convergência é monótona? Que parte de cada método limita a ordem
   perto da parede?
4. Para chegar ao mesmo C<sub>d</sub>, quantos volumes e quanto tempo de computador cada abordagem
   precisa?

**A6. Passo de tempo e custo.**

1. Deduza os limites de estabilidade de um esquema explícito: CFL convectivo,
   U Δt / h ≤ C, e número de difusão, ν Δt / h² ≤ K.
2. Calcule esses limites para a menor célula de cada malha.
3. Explique por que o Δt do Cordel é 10 a 20 vezes menor que os do fvm-cfd, embora
   todos sejam explícitos na convecção.
4. Compare o custo total, em passos × custo por passo, com o tempo que você mediu.

**A7. Discussão crítica.**

1. Liste as fontes de diferença entre o Cordel e as duas soluções do fvm-cfd e estime a
    importância de cada uma:
    - malha;
    - representação da parede;
    - tamanho do domínio e bloqueio D/H;
    - condição no topo e no fundo;
    - referência de pressão;
    - janela de média.
2. Classifique cada comparação que você fez como verificação ou validação.
3. Qual resultado você considera mais confiável, e por quê?

## Parte B — tópicos à escolha (cada grupo escolhe [a definir] tópicos)

**B1. Re = 40, escoamento permanente.**

1. Rode `--Re 40` no fvm-cfd (`cylinder.py` e `cylinder2d.py`) e `--nu 0.025` no Cordel.
2. Verifique que o escoamento fica permanente e simétrico (C<sub>l</sub> → 0).
3. Meça o comprimento da bolha de recirculação L<sub>w</sub>/D e C<sub>d</sub>. Compare com a literatura
   (L<sub>w</sub>/D ≈ 2,2–2,3 e C<sub>d</sub> ≈ 1,5–1,6; confira nas referências).
4. O `cylinder_polar` do Cordel resolve esse caso num domínio muito maior: compare.

**B2. Strouhal contra Reynolds.** Rode Re = 60, 80, 100, 150 (malha grossa basta) e compare
St(Re) com o ajuste empírico de Williamson para o regime laminar.

**B3. Bloqueio e domínio.**

1. Na malha uniforme, varie a altura do domínio com `--domain`, por exemplo ±5, ±7, ±10.
   Na malha O, varie o raio externo com `--Rext 10 20 40 80`.
2. Quantifique o efeito do bloqueio D/H sobre C<sub>d</sub> e St. O repositório também tem um estudo de
   domínio retangular com malha esticada (`cylinder_rect.py`): compare com ele.
3. A diferença entre as soluções se explica por isso?

**B4. Esquema convectivo.** Compare `--scheme quick` e `--scheme central` na malha uniforme.
Calcule o número de Péclet de malha U h / ν e explique o que acontece com as diferenças
centrais.

**B5. Modelo de esteira no método dos painéis (desafio).**

1. Implemente um método de vórtices discretos: vórtices livres soltos de pontos de
   separação fixos, convectados pelo campo induzido, com conservação da circulação total.
2. Obtenha uma esteira qualitativa e um St.
3. Discuta o que esse modelo precisa supor que os códigos viscosos calculam.

**B6. Visualização.** Faça uma animação da vorticidade em um período de desprendimento
(na malha O: `python python/cylinder.py --NT 128 --resume --animate 3`).
Relacione a fase de C<sub>l</sub> à posição dos vórtices.

# 7. Entregáveis

1. **Relatório** em PDF (no máximo [a definir] páginas, sem anexos), com:
    - introdução;
    - formulações;
    - configuração das rodadas: malha, Δt, domínio, tempo de computador;
    - resultados;
    - discussão;
    - conclusões.

    Toda figura deve ter eixos com nome e unidades (ou "adimensional") e ser citada no texto.
2. **Tabela-resumo** com C<sub>d</sub>, C<sub>l</sub>′, St, C<sub>p</sub> de estagnação, C<sub>p</sub> mínimo e sua posição, C<sub>p</sub>
   de base e C<sub>d,p</sub> para as quatro soluções e a literatura.
3. **Scripts** de pós-processamento que o grupo escreveu ou modificou, e os comandos
   exatos que reproduzem cada rodada.
4. [a definir: apresentação oral de N minutos].

# 8. Critérios de avaliação

| Critério | Peso |
|---|---|
| Formulações e hipóteses corretas (A1, A2) | 20% |
| Resultados corretos e reprodutíveis (A3, A4, A5) | 30% |
| Análise das diferenças, com números (A4, A6, A7) | 30% |
| Parte B | 10% |
| Clareza do relatório e qualidade das figuras | 10% |

Respostas qualitativas ("os resultados ficaram próximos") sem números valem pouco. Uma
diferença de 5% entre dois códigos deve ser atribuída a causas concretas e, quando
possível, testada com uma rodada adicional.

# 9. Referências

As referências abaixo foram citadas de memória. Confira os dados bibliográficos e os
valores numéricos nas fontes.

- J. Katz & A. Plotkin, *Low-Speed Aerodynamics*, Cambridge University Press (método dos
  painéis).
- J. D. Anderson, *Fundamentals of Aerodynamics*, McGraw-Hill (escoamento potencial,
  paradoxo de d'Alembert).
- J. H. Ferziger & M. Perić, *Computational Methods for Fluid Dynamics*, Springer (volumes
  finitos, malha desencontrada, métodos de projeção).
- O. C. Zienkiewicz, R. L. Taylor & P. Nithiarasu, *The Finite Element Method for Fluid
  Dynamics*, Butterworth-Heinemann (método CBS).
- C. H. K. Williamson (1996), Vortex dynamics in the cylinder wake, *Annual Review of
  Fluid Mechanics* 28, 477–539.
- R. D. Henderson (1995), Details of the drag curve near the onset of vortex shedding,
  *Physics of Fluids* 7, 2102–2104.
- M. Braza, P. Chassaing & H. Ha Minh (1986), Numerical study and physical analysis of the
  pressure and velocity fields in the near wake of a circular cylinder, *Journal of Fluid
  Mechanics* 165, 79–130.
- J. Park, K. Kwon & H. Choi (1998), Numerical solutions of flow past a circular cylinder
  at Reynolds numbers up to 160, *KSME International Journal* 12, 1200–1205.
- T. Kajishima, S. Takiguchi, H. Hamasaki & Y. Miyake (2001), Turbulence structure of
  particle-laden flow in a vertical plane channel due to vortex shedding, *JSME
  International Journal B* 44 (fronteira imersa por fração de volume).
- J. Kim & P. Moin (1985), Application of a fractional-step method to incompressible
  Navier–Stokes equations, *Journal of Computational Physics* 59, 308–323.
