"""Cilindro pelo metodo dos paineis (fontes) comparado com a solucao exata,
com e sem circulacao imposta (cilindro girante / efeito Magnus)."""
import numpy as np
import matplotlib.pyplot as plt

from _common import out, rel, start
from panelpy import Panels, solve
from panelpy import analytic as an, geometry as geo
from panelpy.plotting import BLUE, ORANGE, cp_axes, streamlines

start("02_cilindro_paineis")

R, U = 1.0, 1.0
fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))

# Cp x theta
th = np.linspace(0, 2 * np.pi, 400)
for Gcw, cor in ((0.0, BLUE), (2 * np.pi, ORANGE)):
    p = Panels.from_points(*geo.circle(24, R))
    s = solve(p, U, circulation=Gcw) if Gcw else solve(p, U, lifting=False)
    tc = np.mod(np.arctan2(p.yc, p.xc), 2 * np.pi)
    axs[0].plot(np.degrees(th), an.cylinder_cp(th, U, R, Gcw), color=cor, lw=1.2,
                label=f"exato, $\\Gamma$={Gcw:.2f}")
    axs[0].plot(np.degrees(tc), s.cp, "o", ms=4, color=cor, mfc="none",
                label=f"paineis N=24, $\\Gamma$={Gcw:.2f}")
cp_axes(axs[0], "$\\theta$ (graus)")
axs[0].set_title("Cp na superficie")
axs[0].legend(fontsize=8)

# convergencia
print("== Cilindro: convergencia ==")
print(f"{'N':>5} {'erro max Cp':>12} {'Cd':>10} {'soma sigma*L':>13} {'Cl (Gamma=2pi)':>15} {'erro Cl':>9}")
Ns = [8, 16, 32, 64, 128, 256]
err_cl = []
for N in Ns:
    p = Panels.from_points(*geo.circle(N, R))
    s = solve(p, U, lifting=False)
    e = np.abs(s.cp - an.cylinder_cp(np.arctan2(p.yc, p.xc), U, R)).max()
    f = s.forces(2 * R)
    s2 = solve(p, U, circulation=2 * np.pi)
    cl = s2.forces(2 * R)["cl"]
    err_cl.append(abs(cl - an.cylinder_cl(U, R, 2 * np.pi)))
    print(f"{N:5d} {e:12.2e} {f['cd']:10.1e} {np.sum(s.sigma * p.length):13.1e} {cl:15.6f} {err_cl[-1]:9.2e}")
axs[1].loglog(Ns, err_cl, "o-", color=BLUE, label="|Cl - Cl exato|, $\\Gamma=2\\pi$")
axs[1].loglog(Ns, err_cl[0] * (np.array(Ns) / Ns[0]) ** -1.0, "--", color="0.5", label="$O(N^{-1})$")
axs[1].set_xlabel("N paineis")
axs[1].set_title("Convergencia (integracao de pressao)")
axs[1].grid(alpha=0.3, which="both")
axs[1].legend(fontsize=8)

p = Panels.from_points(*geo.circle(128, R))
s = solve(p, U, circulation=2 * np.pi)
streamlines(axs[2], s.velocity, (-3, 3), (-2, 2), body=p.polygon(), n=160)
axs[2].set_title("Paineis: cilindro girante ($\\Gamma=2\\pi$)")
fig.tight_layout()
fig.savefig(out("02_cilindro.png"), dpi=130)
print("Figura:", rel(out("02_cilindro.png")))
