# Resultados

Cada subpasta corresponde a um script de `examples/` e contém as figuras e a
saída numérica completa (`saida.txt`). Para regenerar tudo:

```bash
for f in examples/0*.py; do python "$f"; done
```

| Pasta | Script | Conteúdo |
|---|---|---|
| `01_escoamentos_elementares/` | `examples/01_escoamentos_elementares.py` | `panorama.png` com os seis casos; cada caso em sua subpasta (abaixo) |
| `  ├ 01a_dipolo/` | 〃 | Linhas de corrente circulares, decaimento 1/r², circulação nula |
| `  ├ 01b_vortice/` | 〃 | V_θ = Γ/(2πr), circulação Γ em qualquer raio, nula fora do vórtice |
| `  ├ 01c_semicorpo_rankine/` | 〃 | Ponto de estagnação, linha do corpo ψ = m/2, assíntota m/(2U) |
| `  ├ 01d_oval_rankine/` | 〃 | Semicomprimento, semiespessura, corpo fechado |
| `  ├ 01e_cilindro/` | 〃 | ψ = 0 em r = R, Cp = 1 − 4 sen²θ, d'Alembert |
| `  └ 01f_cilindro_circulacao/` | 〃 | Γ = 4πUR: Cl por pressão = Kutta–Joukowski, pontos de estagnação fundidos |
| `02_cilindro_paineis/` | `examples/02_cilindro_paineis.py` | Cilindro por painéis × solução exata, convergência, cilindro girante |
| `03_corpo_rankine/` | `examples/03_corpo_rankine.py` | Oval de Rankine por painéis × superposição analítica, convergência |
| `04_aerofolios/` | `examples/04_aerofolios.py` | Kármán–Trefftz/Joukowski × solução exata, NACA Cl × α, efeito solo, linhas de corrente |
| `05_corpo_ahmed/` | `examples/05_corpo_ahmed.py` | Corpo de Ahmed 2D: corpo fechado, modelo de esteira, sensibilidade, campo de escoamento |
