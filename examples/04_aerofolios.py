"""Aerofolios pelo metodo de Hess-Smith (fontes + vortice + condicao de Kutta).

Validacao:
* Karman-Trefftz / Joukowski: solucao EXATA por transformacao conforme.
* NACA 4 digitos: teoria do aerofolio fino (Cl = 2 pi (alpha - alpha_L0)).
Aplicacao: curvas Cl x alpha e efeito solo (metodo das imagens).

Lembrete: escoamento potencial nao tem estol nem arrasto; as curvas Cl x alpha
sao lineares e so representam a realidade em angulos pequenos (regime colado).
"""
import numpy as np
import matplotlib.pyplot as plt

from _common import out
from panelpy import Panels, rotate, solve
from panelpy import analytic as an, geometry as geo
from panelpy.plotting import BLUE, ORANGE, GRAY, cp_axes, streamlines

ALPHA = np.radians(5.0)

# ------------------------------------------------ 1) Karman-Trefftz exato
print("== Aerofolio de Karman-Trefftz (mx=0.1, my=0.1, tau=10 graus), alpha=5 graus ==")
kt = an.Joukowski(0.1, 0.1, tau_deg=10.0)
print(f"Cl exato = {kt.cl(ALPHA):.5f}")
print(f"{'N':>5} {'Cl pressao':>11} {'erro %':>8} {'Cl K-J':>9} {'erro %':>8} {'Cd':>9}")
Ns = [50, 100, 200, 400, 800]
e_p, e_kj = [], []
for N in Ns:
    p = Panels.from_points(*kt.coordinates(N))
    f = solve(p, alpha=ALPHA).forces()
    e_p.append(abs(f["cl"] / kt.cl(ALPHA) - 1) * 100)
    e_kj.append(abs(f["cl_kj"] / kt.cl(ALPHA) - 1) * 100)
    print(f"{N:5d} {f['cl']:11.5f} {e_p[-1]:8.3f} {f['cl_kj']:9.5f} {e_kj[-1]:8.3f} {f['cd']:9.1e}")

print("\n== Joukowski classico (bordo de fuga em cuspide, tau=0) ==")
jk = an.Joukowski(0.1, 0.1, tau_deg=0.0)
for N in (100, 400, 1600):
    f = solve(Panels.from_points(*jk.coordinates(N)), alpha=ALPHA).forces()
    print(f"N={N:5d}  Cl={f['cl']:.5f}  Cl_KJ={f['cl_kj']:.5f}  exato={jk.cl(ALPHA):.5f}")
print("(a cuspide degrada a convergencia do Hess-Smith com paineis de intensidade constante)")

fig, axs = plt.subplots(2, 2, figsize=(14, 9))
xe, cpe = kt.surface_cp(ALPHA, 600)
p = Panels.from_points(*kt.coordinates(120))
s = solve(p, alpha=ALPHA)
axs[0, 0].plot(xe, cpe, color=BLUE, lw=1.2, label="exato (transformacao conforme)")
axs[0, 0].plot(p.xc, s.cp, "o", ms=3.5, mfc="none", color=ORANGE, label="Hess-Smith N=120")
cp_axes(axs[0, 0])
axs[0, 0].legend(fontsize=8)
axs[0, 0].set_title("Karman-Trefftz, $\\alpha$=5$^\\circ$")
axs[0, 1].loglog(Ns, e_p, "o-", color=BLUE, label="Cl por pressao")
axs[0, 1].loglog(Ns, e_kj, "s-", color=ORANGE, label="Cl por Kutta-Joukowski")
axs[0, 1].loglog(Ns, e_p[0] * (np.array(Ns) / Ns[0]) ** -1.0, "--", color=GRAY, label="$O(N^{-1})$")
axs[0, 1].set_xlabel("N paineis")
axs[0, 1].set_ylabel("erro em Cl (%)")
axs[0, 1].grid(alpha=0.3, which="both")
axs[0, 1].legend(fontsize=8)
axs[0, 1].set_title("Convergencia vs. solucao exata")

# ------------------------------------------------ 2) NACA: Cl x alpha
print("\n== NACA 4 digitos (N=200): Cl x alpha ==")
alphas = np.radians(np.arange(-6, 13, 2.0))
cores = {"0012": BLUE, "2412": ORANGE, "4412": "#2a9d55"}
tabela = []
for code, cor in cores.items():
    p = Panels.from_points(*geo.naca4(code, 200))
    cl = np.array([solve(p, alpha=a).forces()["cl"] for a in alphas])
    slope, inter = np.polyfit(alphas, cl, 1)
    a0_pan = -inter / slope
    a0_taf = an.thin_airfoil_alpha0(geo.naca4_camber(code)[1])
    tabela.append((code, slope, np.degrees(a0_pan), np.degrees(a0_taf)))
    axs[1, 0].plot(np.degrees(alphas), cl, "o-", color=cor, ms=4, label=f"NACA {code} (paineis)")
    axs[1, 0].plot(np.degrees(alphas), an.thin_airfoil_cl(alphas, a0_taf), "--", color=cor, lw=1)
print(f"{'perfil':>8} {'dCl/da (1/rad)':>15} {'a_L0 paineis':>13} {'a_L0 aerof. fino':>17}")
for code, sl, a0p, a0t in tabela:
    print(f"{code:>8} {sl:15.4f} {a0p:12.3f}° {a0t:16.3f}°")
print(f"(teoria do aerofolio fino: dCl/da = 2 pi = {2*np.pi:.4f}/rad)")
axs[1, 0].set_xlabel("$\\alpha$ (graus)")
axs[1, 0].set_ylabel("$C_l$")
axs[1, 0].grid(alpha=0.3)
axs[1, 0].legend(fontsize=8, title="tracejado: aerofolio fino", title_fontsize=8)
axs[1, 0].set_title("Curvas de sustentacao (potencial: sem estol)")

# ------------------------------------------------ 3) efeito solo
print("\n== NACA 4412, alpha=4 graus: efeito solo (altura do bordo de fuga h/c) ==")
xa, ya = geo.naca4("4412", 200)
xr, yr = rotate(xa, ya, np.radians(4.0), 0.25, 0.0)
f_inf = solve(Panels.from_points(xr, yr)).forces()
hs = np.array([0.1, 0.15, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0])
cls = []
for h in hs:
    yy = yr - yr[0] + h  # bordo de fuga a altura h do solo (y = 0)
    cls.append(solve(Panels.from_points(xr, yy), ground=0.0).forces()["cl"])
    print(f"h/c = {h:4.2f}: Cl = {cls[-1]:.4f}   (fora do solo: {f_inf['cl']:.4f})")
axs[1, 1].plot(hs, cls, "o-", color=BLUE, label="com solo (imagens)")
axs[1, 1].axhline(f_inf["cl"], color=GRAY, ls="--", label="fora do solo")
axs[1, 1].set_xlabel("altura do bordo de fuga h/c")
axs[1, 1].set_ylabel("$C_l$")
axs[1, 1].grid(alpha=0.3)
axs[1, 1].legend(fontsize=8)
axs[1, 1].set_title("NACA 4412, $\\alpha$=4$^\\circ$: efeito solo")
fig.tight_layout()
fig.savefig(out("04_aerofolios.png"), dpi=130)

# linhas de corrente
fig, axs = plt.subplots(1, 2, figsize=(15, 4.5))
p = Panels.from_points(*geo.naca4("4412", 200))
s = solve(p, alpha=np.radians(8))
streamlines(axs[0], s.velocity, (-0.5, 1.6), (-0.5, 0.6), body=p.polygon(), n=260,
            U=1.0, density=2.0)
axs[0].set_title("NACA 4412, $\\alpha$=8$^\\circ$")
p = Panels.from_points(xr, yr - yr[0] + 0.15)
s = solve(p, ground=0.0)
streamlines(axs[1], s.velocity, (-0.5, 1.6), (-0.05, 0.6), body=p.polygon(), ground=0.0,
            n=260, density=2.0)
axs[1].set_title("NACA 4412, $\\alpha$=4$^\\circ$, h/c=0.15 (efeito solo)")
fig.tight_layout()
fig.savefig(out("04_aerofolios_campo.png"), dpi=130)
print("Figuras:", out("04_aerofolios.png"), out("04_aerofolios_campo.png"))
