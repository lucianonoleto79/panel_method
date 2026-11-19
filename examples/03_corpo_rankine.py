"""Oval de Rankine: o contorno exato (psi = 0) e discretizado em paineis e a
distribuicao de pressao calculada pelo metodo das fontes e comparada com a
superposicao analitica fonte + sumidouro + uniforme."""
import numpy as np
import matplotlib.pyplot as plt

from _common import out, rel, start
from panelpy import Panels, solve
from panelpy import analytic as an, geometry as geo
from panelpy.elementary import rankine_flow
from panelpy.plotting import BLUE, ORANGE, cp_axes, streamlines

start("03_corpo_rankine")

U, m, a = 1.0, 2 * np.pi, 1.0
l = an.rankine_half_length(U, m, a)
h = an.rankine_half_thickness(U, m, a)

fig, axs = plt.subplots(1, 3, figsize=(17, 4.6))
p = Panels.from_points(*geo.rankine_oval(40, U, m, a))
s = solve(p, U, lifting=False)
phi = np.linspace(0, np.pi, 300)
r = an.rankine_contour_radius(phi, U, m, a)
axs[0].plot(r * np.cos(phi), an.rankine_cp(phi, U, m, a), color=BLUE, label="exato (superposicao)")
up = p.yc > 0
axs[0].plot(p.xc[up], s.cp[up], "o", color=ORANGE, mfc="none", ms=5, label="paineis N=40")
cp_axes(axs[0], "x")
axs[0].legend(fontsize=8)
axs[0].set_title(f"Oval de Rankine: l={l:.4f}, h={h:.4f}")

print("== Oval de Rankine: convergencia ==")
print(f"{'N':>5} {'erro max Cp':>12} {'Cd':>10} {'Cl':>10}")
Ns = [16, 32, 64, 128, 256]
errs = []
for N in Ns:
    p = Panels.from_points(*geo.rankine_oval(N, U, m, a))
    s = solve(p, U, lifting=False)
    phic = np.arctan2(p.yc, p.xc)
    e = np.abs(s.cp - an.rankine_cp(phic, U, m, a)).max()
    f = s.forces(2 * l)
    errs.append(e)
    print(f"{N:5d} {e:12.2e} {f['cd']:10.1e} {f['cl']:10.1e}")
axs[1].loglog(Ns, errs, "o-", color=BLUE, label="erro max de Cp")
axs[1].loglog(Ns, errs[0] * (np.array(Ns) / Ns[0]) ** -1.0, "--", color="0.5", label="$O(N^{-1})$")
axs[1].set_xlabel("N paineis")
axs[1].grid(alpha=0.3, which="both")
axs[1].legend(fontsize=8)
axs[1].set_title("Convergencia")

# linhas de corrente: paineis (cor) x exato (contorno de psi=0, tracejado)
p = Panels.from_points(*geo.rankine_oval(128, U, m, a))
s = solve(p, U, lifting=False)
streamlines(axs[2], s.velocity, (-3.5, 3.5), (-2.5, 2.5), body=p.polygon(), n=160)
axs[2].plot([-a, a], [0, 0], "x", color="k", zorder=5)
axs[2].set_title("Paineis (N=128); x = fonte/sumidouro analiticos")
fig.tight_layout()
fig.savefig(out("03_rankine.png"), dpi=130)
print("Figura:", rel(out("03_rankine.png")))
