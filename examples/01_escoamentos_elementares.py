"""Casos classicos por superposicao analitica de escoamentos elementares.

Cada caso e gravado em uma subpasta propria de results/01_escoamentos_elementares/
(figura + saida.txt com as verificacoes); o panorama com os seis casos fica na
raiz dessa pasta.
"""
import numpy as np
import matplotlib.pyplot as plt

from _common import out, rel, start
from panelpy import analytic as an
from panelpy.elementary import (Doublet, Source, Uniform, Vortex, cylinder_flow,
                                rankine_flow)
from panelpy.plotting import BLUE, ORANGE, cp_axes

BASE = "01_escoamentos_elementares"
U, R = 1.0, 1.0
X, Y = np.meshgrid(np.linspace(-3.5, 3.5, 501), np.linspace(-2.5, 2.5, 361))
TH = np.linspace(0, 2 * np.pi, 2001)[:-1]
DTH = 2 * np.pi / TH.size


def plot_psi(ax, flow, titulo, levels=np.linspace(-4, 4, 41), cilindro=False, marcas=(),
             psi_corpo=(0.0,), simetrico=False):
    """Linhas de corrente (psi = cte); as linhas do corpo (psi_corpo) em laranja.

    simetrico=True calcula psi so em y >= 0 e usa psi(x,-y) = -psi(x,y), evitando
    o corte de ramo de atan2 no semieixo x < 0 (fonte isolada).
    """
    with np.errstate(all="ignore"):
        if simetrico:
            psi = np.sign(Y) * flow.stream_function(X, np.abs(Y))
        else:
            psi = flow.stream_function(X, Y)
    ax.contour(X, Y, psi, levels=levels, colors="k", linewidths=0.6)
    ax.contour(X, Y, psi, levels=sorted(psi_corpo), colors=ORANGE, linewidths=1.6,
               linestyles="solid")
    if cilindro:
        ax.fill(R * np.cos(TH), R * np.sin(TH), color="0.85", ec="k", zorder=3)
    for (x, y, m) in marcas:
        ax.plot(x, y, m, color="k", ms=7, zorder=4)
    ax.set_title(titulo, fontsize=10)
    ax.set_aspect("equal")
    ax.set_xlim(-3.5, 3.5)
    ax.set_ylim(-2.5, 2.5)


def circulacao(flow, xc, yc, r):
    """Circulacao por integral de linha num circulo de raio r centrado em (xc, yc)."""
    u, v = flow.velocity(xc + r * np.cos(TH), yc + r * np.sin(TH))
    return np.sum(-u * np.sin(TH) + v * np.cos(TH)) * r * DTH


def figura_unica(nome, flow, titulo, **kw):
    fig, ax = plt.subplots(figsize=(7.5, 5.4))
    plot_psi(ax, flow, titulo, **kw)
    fig.tight_layout()
    fig.savefig(out(nome), dpi=130)
    plt.close(fig)
    print("Figura:", rel(out(nome)))


casos = []  # (titulo, flow, kwargs) para o panorama

# ------------------------------------------------------------------ dipolo
start(f"{BASE}/01a_dipolo")
kappa = 2 * np.pi
flow = Doublet(kappa)
titulo = "Dipolo ($\\kappa=2\\pi$)"
print("== Dipolo: phi = kappa cos(theta)/(2 pi r), psi = -kappa sin(theta)/(2 pi r) ==")
for c in (0.5, 1.0, 2.0):
    # linhas de corrente sao circulos tangentes a origem: x^2 + (y + c)^2 = c^2
    x, y = c * np.cos(TH), -c + c * np.sin(TH)
    ok = np.abs(y) > 1e-9
    psi = flow.stream_function(x[ok], y[ok])
    print(f"circulo de raio {c}: psi = {psi.mean():.6f} +- {psi.std():.1e}"
          f"  (esperado kappa/(4 pi c) = {kappa / (4 * np.pi * c):.6f})")
for r in (1.0, 2.0, 4.0):
    u, v = flow.velocity(r, 0.0)
    print(f"|V| em (r={r}, 0) = {np.hypot(u, v):.6f}  (esperado kappa/(2 pi r^2) = {kappa / (2 * np.pi * r**2):.6f})")
print(f"circulacao em torno do dipolo = {circulacao(flow, 0, 0, 1.0):.1e} (esperado 0)")
figura_unica("dipolo.png", flow, titulo, marcas=[(0, 0, "+")])
casos.append((titulo, flow, dict(marcas=[(0, 0, "+")])))

# ------------------------------------------------------------------ vortice
start(f"{BASE}/01b_vortice")
G = 2 * np.pi
flow = Vortex(G)
titulo = "Vortice ($\\Gamma=2\\pi$, anti-horario)"
print("== Vortice livre: V_theta = Gamma/(2 pi r) ==")
for r in (0.5, 1.0, 2.0):
    print(f"circulacao, circulo de raio {r} em torno do vortice: {circulacao(flow, 0, 0, r):.6f} (esperado {G:.6f})")
print(f"circulacao num circulo que NAO envolve o vortice: {circulacao(flow, 2.5, 0, 1.0):.1e}"
      " (esperado 0: escoamento irrotacional fora da singularidade)")
for r in (0.5, 1.0, 2.0):
    u, v = flow.velocity(r, 0.0)
    print(f"V_theta(r={r}) = {v:.6f}  (esperado {G / (2 * np.pi * r):.6f})")
figura_unica("vortice.png", flow, titulo, marcas=[(0, 0, "+")])
casos.append((titulo, flow, dict(marcas=[(0, 0, "+")])))

# ------------------------------------------------------- semi-corpo de Rankine
start(f"{BASE}/01c_semicorpo_rankine")
m = 2 * np.pi
flow = Uniform(U) + Source(m)
titulo = "Semi-corpo de Rankine (uniforme + fonte)"
xs = -m / (2 * np.pi * U)
u, v = flow.velocity(xs, 0.0)
print("== Semi-corpo de Rankine: uniforme U + fonte m na origem ==")
print(f"ponto de estagnacao x_s = -m/(2 pi U) = {xs:.6f};  V(x_s, 0) = ({u:.1e}, {v:.1e})")
psi_s = flow.stream_function(xs, 1e-12)
print(f"psi no ponto de estagnacao = {psi_s:.6f} (linha do corpo: psi = m/2 = {m / 2:.6f})")
# semiespessura assintotica m/(2U): resolve psi(x, y) = m/2 longe a jusante
for xf in (10.0, 100.0, 1000.0):
    lo, hi = 0.0, m / U
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if flow.stream_function(xf, mid) < m / 2 else (lo, mid)
    print(f"semiespessura em x = {xf:7.1f}: {0.5 * (lo + hi):.6f}  (assintota m/(2U) = {m / (2 * U):.6f})")
d = np.linspace(0.2, 4.0, 20)
kw = dict(levels=np.sort(np.r_[-m / 2 - d, np.linspace(-m / 2, m / 2, 17)[1:-1], m / 2 + d]),
          marcas=[(0, 0, "x")], psi_corpo=(-m / 2, m / 2), simetrico=True)
figura_unica("semicorpo_rankine.png", flow, titulo, **kw)
casos.append((titulo, flow, kw))

# ------------------------------------------------------------ oval de Rankine
start(f"{BASE}/01d_oval_rankine")
a = 1.0
flow = rankine_flow(U, m, a)
titulo = "Oval de Rankine (m=2$\\pi$, a=1)"
l = an.rankine_half_length(U, m, a)
h = an.rankine_half_thickness(U, m, a)
u, _ = flow.velocity(np.array([-l, l]), np.zeros(2))
print("== Oval de Rankine: uniforme + fonte em (-a,0) + sumidouro em (+a,0) ==")
print(f"semicomprimento l = a sqrt(1 + m/(pi U a)) = {l:.6f};  u(+-l, 0) = {np.abs(u).max():.1e}")
print(f"semiespessura h = {h:.6f};  psi(0, h) = {flow.stream_function(0.0, h):.1e}")
print(f"vazao liquida (fonte + sumidouro) = {m - m:.1f} -> corpo fechado")
print(f"circulacao em torno do corpo = {circulacao(flow, 0, 0, 2.5):.1e} (esperado 0)")
figura_unica("oval_rankine.png", flow, titulo, marcas=[(-a, 0, "x"), (a, 0, "x")])
casos.append((titulo, flow, dict(marcas=[(-a, 0, "x"), (a, 0, "x")])))


# ------------------------------------------------------------------ cilindros
def caso_cilindro(pasta, Gcw, titulo):
    start(f"{BASE}/{pasta}")
    flow = cylinder_flow(U, R, Gcw)
    print(f"== {titulo.replace('$', '').replace(chr(92), '')} ==")
    psi = flow.stream_function(R * np.cos(TH), R * np.sin(TH))
    print(f"variacao de psi na superficie r=R: {np.ptp(psi):.1e} (esperado 0: e linha de corrente)")
    u, v = flow.velocity(R * np.cos(TH), R * np.sin(TH))
    vr = u * np.cos(TH) + v * np.sin(TH)
    print(f"max |V_r| na superficie: {np.abs(vr).max():.1e} (esperado 0)")
    cp = an.cylinder_cp(TH, U, R, Gcw)
    cl = -np.sum(cp * np.sin(TH)) * DTH / 2
    cd = -np.sum(cp * np.cos(TH)) * DTH / 2
    print(f"Cl (integracao de pressao) = {cl:.6f};  Kutta-Joukowski Gamma/(U R) = {an.cylinder_cl(U, R, Gcw):.6f}")
    print(f"Cd = {cd:.1e} (paradoxo de d'Alembert)")
    s = -Gcw / (4 * np.pi * U * R)
    if abs(s) <= 1:
        th_s = np.degrees([np.arcsin(s), np.pi - np.arcsin(s)])
        print(f"pontos de estagnacao: sin(theta) = -Gamma/(4 pi U R) = {s:.4f} -> theta = "
              f"{np.mod(th_s[0], 360):.2f} e {np.mod(th_s[1], 360):.2f} graus")
    fig, axs = plt.subplots(1, 2, figsize=(13, 4.8))
    plot_psi(axs[0], flow, titulo, cilindro=True)
    axs[1].plot(np.degrees(TH), cp, color=BLUE)
    cp_axes(axs[1], "$\\theta$ (graus)")
    axs[1].set_title("Cp na superficie (exato)")
    fig.tight_layout()
    fig.savefig(out(f"{pasta[4:]}.png"), dpi=130)
    plt.close(fig)
    print("Figura:", rel(out(f"{pasta[4:]}.png")))
    casos.append((titulo, flow, dict(cilindro=True)))


caso_cilindro("01e_cilindro", 0.0, "Cilindro: uniforme + dipolo")
caso_cilindro("01f_cilindro_circulacao", 4 * np.pi,
              "Cilindro com circulacao $\\Gamma=4\\pi U R$")

# ------------------------------------------------------------------ panorama
start(BASE)
print("Panorama dos seis escoamentos elementares; detalhes nas subpastas:")
fig, axs = plt.subplots(2, 3, figsize=(15, 8.4))
for ax, (titulo, flow, kw) in zip(axs.ravel(), casos):
    plot_psi(ax, flow, titulo, **kw)
fig.suptitle("Escoamentos potenciais elementares (linhas de corrente $\\psi$ = cte; laranja: linha do corpo)")
fig.tight_layout()
fig.savefig(out("panorama.png"), dpi=130)
for sub in ("01a_dipolo", "01b_vortice", "01c_semicorpo_rankine", "01d_oval_rankine",
            "01e_cilindro", "01f_cilindro_circulacao"):
    print(f"  {BASE}/{sub}/")
print("Figura:", rel(out("panorama.png")))
