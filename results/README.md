# Resultados

Cada subpasta corresponde a um script de `examples/` e contém as figuras e a
saída numérica completa (`saida.txt`). Para regenerar tudo:

```bash
for f in examples/0*.py; do python "$f"; done
```

| Pasta | Script | Conteúdo |
|---|---|---|
| `01_escoamentos_elementares/` | `examples/01_escoamentos_elementares.py` | Dipolo, vórtice, semicorpo e oval de Rankine, cilindro com/sem circulação (soluções analíticas) e verificações |
| `02_cilindro_paineis/` | `examples/02_cilindro_paineis.py` | Cilindro por painéis × solução exata, convergência, cilindro girante |
| `03_corpo_rankine/` | `examples/03_corpo_rankine.py` | Oval de Rankine por painéis × superposição analítica, convergência |
| `04_aerofolios/` | `examples/04_aerofolios.py` | Kármán–Trefftz/Joukowski × solução exata, NACA Cl × α, efeito solo, linhas de corrente |
| `05_corpo_ahmed/` | `examples/05_corpo_ahmed.py` | Corpo de Ahmed 2D: corpo fechado, modelo de esteira, sensibilidade, campo de escoamento |
