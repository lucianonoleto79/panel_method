"""Casos classicos por superposicao analitica de escoamentos elementares:
dipolo, vortice, semi-corpo de Rankine, oval de Rankine e cilindro com e sem
circulacao (efeito Magnus)."""
import numpy as np
import matplotlib.pyplot as plt

from _common import out, rel, start
from panelpy import analytic as an
from panelpy.elementary import (Doublet, Source, Uniform, Vortex, cylinder_flow,
                                rankine_flow)

start("01_escoamentos_elementares")

U = 1.0
casos = [
    ("Dipolo ($\\kappa=2\\pi$)", Doublet(2 * np.pi), None),
    ("Vortice ($\\Gamma=2\\pi$, anti-horario)", Vortex(2 * np.pi), None),
    ("Semi-corpo de Rankine (uniforme + fonte)", Uniform(U) + Source(2 * np.pi), None),
    ("Oval de Rankine (m=2$\\pi$, a=1)", rankine_flow(U, 2 * np.pi, 1.0), "rankine"),
    ("Cilindro: uniforme + dipolo", cylinder_flow(U, 1.0), "cyl"),
    ("Cilindro com circulacao $\\Gamma=4\\pi U R$", cylinder_flow(U, 1.0, 4 * np.pi), "cyl"),
]

X, Y = np.meshgrid(np.linspace(-3.5, 3.5, 501), np.linspace(-2.5, 2.5, 361))
fig, axs = plt.subplots(2, 3, figsize=(15, 7.8))
for ax, (titulo, flow, corpo) in zip(axs.ravel(), casos):
    with np.errstate(all="ignore"):
        psi = flow.stream_function(X, Y)
    ax.contour(X, Y, psi, levels=np.linspace(-4, 4, 41), colors="k", linewidths=0.6)
    ax.contour(X, Y, psi, levels=[0.0], colors="#d9730d", linewidths=1.6)
    if corpo == "cyl":
        th = np.linspace(0, 2 * np.pi, 200)
        ax.fill(np.cos(th), np.sin(th), color="0.85", ec="k", zorder=3)
    ax.set_title(titulo, fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-2.5, 2.5)
fig.suptitle("Escoamentos potenciais elementares (linhas de corrente $\\psi$ = cte; laranja: $\\psi=0$)")
fig.tight_layout()
fig.savefig(out("01_elementares.png"), dpi=130)

# ---------------------------------------------------------------- validacoes
print("== Verificacoes analiticas ==")
th = np.linspace(0, 2 * np.pi, 721)
psi_cyl = cylinder_flow(U, 1.0).stream_function(np.cos(th), np.sin(th))
print(f"Cilindro: max |psi| na superficie r=R          = {np.abs(psi_cyl).max():.2e} (esperado 0)")

l = an.rankine_half_length(U, 2 * np.pi, 1.0)
h = an.rankine_half_thickness(U, 2 * np.pi, 1.0)
u_est, _ = rankine_flow(U, 2 * np.pi, 1.0).velocity(l, 0.0)
print(f"Rankine: semi-comprimento l = {l:.5f}, semi-espessura h = {h:.5f}, u(l,0) = {u_est:.1e}")

# circulacao do vortice por integral de linha
G = 2 * np.pi
u, v = Vortex(G).velocity(2 * np.cos(th[:-1]), 2 * np.sin(th[:-1]))
circ = np.sum((-u * np.sin(th[:-1]) + v * np.cos(th[:-1])) * 2 * (th[1] - th[0]))
print(f"Vortice: circulacao numerica = {circ:.6f} (esperado {G:.6f})")

# Magnus: integra a pressao na superficie
for Gcw in (0.0, 2 * np.pi, 4 * np.pi):
    cp = an.cylinder_cp(th[:-1], U, 1.0, Gcw)
    dth = th[1] - th[0]
    cl = -np.sum(cp * np.sin(th[:-1])) * dth / 2.0
    cd = -np.sum(cp * np.cos(th[:-1])) * dth / 2.0
    print(f"Cilindro Gamma={Gcw:6.3f}: Cl(pressao)={cl:.5f}  Cl(Kutta-Joukowski)={an.cylinder_cl(U,1,Gcw):.5f}  Cd={cd:.1e}")
print("Figura:", rel(out("01_elementares.png")))
