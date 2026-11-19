"""Corpo de Ahmed 2D (perfil lateral) por escoamento potencial: estimativa de
sustentacao sob hipoteses explicitas. Ver avisos em panelpy/ahmed.py e no README.

Coeficientes por unidade de envergadura, referidos a q * L (L = comprimento).
"""
import numpy as np
import matplotlib.pyplot as plt

from _common import out, rel, start
from panelpy import Panels, solve
from panelpy.ahmed import AHMED, ahmed_forces, ahmed_profile
from panelpy.plotting import BLUE, ORANGE, GRAY, streamlines

start("05_corpo_ahmed")

N = 600
L, H, hg = AHMED["L"], AHMED["H"], AHMED["clearance"]
GROUND = -hg


def run(slant, sep, ground=True, lifting=True, **kw):
    x, y, info = ahmed_profile(slant, N, separation=sep, **kw)
    p = Panels.from_points(x, y)
    s = solve(p, lifting=lifting, ground=GROUND if ground else None)
    return p, s, info


# ------------------------------------------------ 1) corpo fechado
print("== Corpo fechado (sem esteira) ==")
print(f"{'rampa':>6} {'Kutta':>8} {'Cl livre':>9} {'Cd livre':>9} {'Cl solo':>9} {'Cd solo':>9}")
for slant in (0, 25, 35):
    for kut in ("nenhum", "top", "bottom"):
        x, y, info = ahmed_profile(slant, N, separation=None,
                                   kutta_at="top" if kut == "nenhum" else kut)
        p = Panels.from_points(x, y)
        lift = kut != "nenhum"
        f1 = solve(p, lifting=lift).forces(L)
        f2 = solve(p, lifting=lift, ground=GROUND).forces(L)
        print(f"{slant:5d}° {kut:>8} {f1['cl']:9.3f} {f1['cd']:9.3f} {f2['cl']:9.3f} {f2['cd']:9.3f}")
print("-> sem circulacao: forca nula fora do solo (d'Alembert) e succao para o solo")
print("   (Cl<0) com solo. Com Kutta em UMA aresta o Cl muda de sinal conforme a")
print("   aresta escolhida: o problema esta mal posto para um corpo de base rombuda.")

# ------------------------------------------------ 2) corpo equivalente
print("\n== Modelo de corpo equivalente (esteira fechada, Kutta na ponta), com solo ==")
casos = [(0, "slant_end"), (12.5, "slant_end"), (20, "slant_end"), (25, "slant_end"),
         (30, "slant_end"), (35, "roof_edge"), (40, "roof_edge")]
print(f"{'rampa':>6} {'separacao':>10} {'Lw (m)':>7} {'Cl_eq':>7} {'Cl_corpo':>9} {'Cl livre':>9} {'Cp_base':>8}")
res = []
for slant, sep in casos:
    p, s, info = run(slant, sep, ground=True)
    f = ahmed_forces(s, info)
    _, s0, _ = run(slant, sep, ground=False)
    f0 = ahmed_forces(s0, info)
    res.append((slant, f["cl_body"], f0["cl_body"]))
    print(f"{slant:5.1f}° {sep:>10} {info['wake_length']:7.3f} {f['cl_eq']:7.3f} "
          f"{f['cl_body']:9.3f} {f0['cl_body']:9.3f} {f['cp_base']:8.3f}")

# sensibilidade as hipoteses da esteira (rampa 25 graus)
print("\n== Sensibilidade: rampa 25°, separacao no fim da rampa, com solo ==")
print(f"{'Lw (m)':>7} {'y_ponta (m)':>12} {'Cl_corpo':>9}")
sens = []
for Lw in (0.15, 0.29, 0.416, 0.6, 1.0):
    for yt in (0.0, 0.05, 0.10):
        p, s, info = run(25, "slant_end", wake_length=Lw, wake_tip_y=yt)
        cl = ahmed_forces(s, info)["cl_body"]
        sens.append((Lw, yt, cl))
        print(f"{Lw:7.3f} {yt:12.3f} {cl:9.3f}")

# correcao grosseira de razao de aspecto (asa de baixo alongamento, Helmbold)
W = 0.389  # largura do corpo de Ahmed (m) - a confirmar no artigo original
AR = W / L
k = AR / (np.sqrt(AR**2 + 4.0) + 2.0)  # (dCL/da)_3D / 2 pi  (Helmbold)
print(f"\nFator de Helmbold para AR = W/L = {AR:.3f}: {k:.3f}")
print("Cl 3D 'grosseiro' (Cl_corpo 2D com solo x fator) — apenas ordem de grandeza:")
for slant, cl, _ in res:
    print(f"  rampa {slant:5.1f}°: {cl * k:6.3f}")

# ------------------------------------------------ figuras
fig, axs = plt.subplots(3, 1, figsize=(12, 11))
for ax, (slant, sep) in zip(axs, ((0, None), (25, "slant_end"), (35, "roof_edge"))):
    if sep is None:
        p, s, info = run(slant, None, ground=True, lifting=False)
        titulo = "Rampa 0°, corpo fechado sem circulacao (com solo)"
    else:
        p, s, info = run(slant, sep, ground=True)
        titulo = f"Rampa {slant}°, corpo equivalente ({sep}), com solo"
    streamlines(ax, s.velocity, (-0.4, 2.0), (GROUND, 0.75), body=p.polygon(),
                ground=GROUND, n=320, density=2.2)
    w = ~info["wet"]
    if w.any():
        ax.plot(np.r_[p.xa[w], p.xb[w][-1]], np.r_[p.ya[w], p.yb[w][-1]], "--", color=ORANGE, lw=2.2, zorder=5)
    ax.set_title(titulo + " — tracejado laranja: fronteira modelada da esteira", fontsize=10)
fig.tight_layout()
fig.savefig(out("05_ahmed_campo.png"), dpi=130)

fig, axs = plt.subplots(1, 2, figsize=(14, 4.5))
axs[0].plot([r[0] for r in res], [r[1] for r in res], "o-", color=BLUE, label="com solo")
axs[0].plot([r[0] for r in res], [r[2] for r in res], "s--", color=ORANGE, label="sem solo")
axs[0].axvline(32.5, color=GRAY, lw=0.8, ls=":")
axs[0].text(33, 0.2, "troca de hipotese\nde separacao", fontsize=8, color=GRAY)
axs[0].set_xlabel("angulo da rampa (graus)")
axs[0].set_ylabel("$C_l$ 2D (ref. q L)")
axs[0].set_ylim(bottom=0)
axs[0].grid(alpha=0.3)
axs[0].legend(fontsize=8)
axs[0].set_title("Corpo de Ahmed 2D: Cl do corpo (modelo de esteira)")
for yt, cor in zip((0.0, 0.05, 0.10), (BLUE, ORANGE, "#2a9d55")):
    d = [(a, c) for a, b, c in sens if b == yt]
    axs[1].plot(*zip(*d), "o-", color=cor, label=f"y ponta = {yt:.2f} m")
axs[1].set_xlabel("comprimento da esteira modelada Lw (m)")
axs[1].set_ylabel("$C_l$ 2D")
axs[1].axhline(0, color=GRAY, lw=0.8)
axs[1].grid(alpha=0.3)
axs[1].legend(fontsize=8)
axs[1].set_title("Sensibilidade as hipoteses (rampa 25°, com solo)")
fig.tight_layout()
fig.savefig(out("05_ahmed_cl.png"), dpi=130)
print("Figuras:", rel(out("05_ahmed_campo.png")), rel(out("05_ahmed_cl.png")))
