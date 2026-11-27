"""Trabalho de CFD: cilindro 2D, Re = 100 -- comparacao Cordel x fvm-cfd x paineis.

Le os resultados ja gerados pelos codigos viscosos (repositorios irmaos)
e resolve aqui o escoamento potencial pelo metodo dos paineis:

  Cordel  (FEM/CBS):        <cordel>/resultados/<saida>/resultados/<saida>.npz e _forcas.npz
                            python cordel2008v4.py cylinder --tmax 200 --saida cylinder_trabalho
  fvm-cfd, malha O:         <fvm>/results/cilindro_O_Re100_NT128_R40/ (campos.npz, cp_parede.csv)
                            python python/cylinder.py --NT 128
  fvm-cfd, malha uniforme   <fvm>/results/cilindro_mac_Re100_h0.04/ (campos.npz, cp_parede.csv)
  + fronteira imersa:       python python/cylinder2d.py --Re 100 --h 0.04
  paineis:                  panelpy.solve (fontes, sem circulacao), N = 64

Sem o campos.npz (nao versionado), usa so o Cp (cp_parede.csv) e o resultados.json.

Grava em results/06_trabalho_cilindro/: comparacao_cp.png, comparacao_forcas.png,
comparacao_vorticidade.png, estudo_malha_fvm.png e saida.txt (tabela).

    python trabalho_cilindro/comparar.py [--cordel ../../cordel] [--fvm ../../fvm-cfd]
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "examples"))
from _common import ROOT, out, rel, start  # noqa: E402  (configura matplotlib e sys.path)

import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Circle  # noqa: E402

from panelpy import Panels, solve  # noqa: E402
from panelpy import analytic as an, geometry as geo  # noqa: E402
from panelpy.plotting import BLUE, GRAY, ORANGE  # noqa: E402

GREEN = "#1f9e6a"
PROJ = os.path.dirname(os.path.dirname(ROOT))          # projetos_cfd/
LIT = dict(St=(0.164, 0.165), Cd=(1.33, 1.35))         # Williamson (1996); Park et al. (1998), Henderson (1995)


def strouhal(t, cl):
    """Frequencia dominante de Cl (FFT com zero-padding): mesma rotina dos dois codigos.
    nan se o escoamento for permanente (amplitude de Cl < 0,02, como no Cordel): o pico seria ruido."""
    if len(t) < 8 or np.ptp(cl) / 2 < 0.02:
        return float("nan")
    ti = np.linspace(t[0], t[-1], 4096)
    c = np.interp(ti, t, cl)
    c = (c - c.mean()) * np.hanning(len(c))
    F = np.abs(np.fft.rfft(c, 8 * len(c)))
    f = np.fft.rfftfreq(8 * len(c), ti[1] - ti[0])
    return float(f[1:][np.argmax(F[1:])])


def espectro(t, cl):
    ti = np.linspace(t[0], t[-1], 4096)
    c = np.interp(ti, t, cl)
    c = (c - c.mean()) * np.hanning(len(c))
    F = np.abs(np.fft.rfft(c, 8 * len(c)))
    return np.fft.rfftfreq(8 * len(c), ti[1] - ti[0]), F / F.max()


def resumo_forcas(t, cd, cl, t_media):
    s = t >= t_media
    return dict(cd_med=cd[s].mean(), cl_med=cl[s].mean(), cl_rms=np.sqrt(np.mean((cl[s] - cl[s].mean()) ** 2)),
                cl_amp=np.ptp(cl[s]) / 2, cd_rms=np.sqrt(np.mean((cd[s] - cd[s].mean()) ** 2)),
                st=strouhal(t[s], cl[s]))


def resumo_cp(theta, cp):
    """theta medido a partir de +x (base); devolve o angulo a partir da estagnacao
    e o arrasto de pressao Cd_p = -1/2 * integral(Cp cos(theta) dtheta) (R/D = 1/2)."""
    th = np.mod(np.asarray(theta), 2 * np.pi)
    k = np.argsort(th)
    thp = np.r_[th[k], th[k][0] + 2 * np.pi]
    cpp = np.r_[np.asarray(cp)[k], np.asarray(cp)[k][0]]
    cd_p = -0.5 * np.trapezoid(cpp * np.cos(thp), thp)
    ang = np.degrees(np.abs(np.pi - np.asarray(theta)))
    o = np.argsort(ang)
    ang, cp = ang[o], np.asarray(cp)[o]
    return dict(ang=ang, cp=cp, cp_estag=cp[0], cp_min=cp.min(), ang_cp_min=ang[np.argmin(cp)], cp_base=cp[-1],
                cd_p=cd_p)


def ler_cordel(pasta, saida):
    base = os.path.join(pasta, "resultados", saida, "resultados", saida)
    if not os.path.exists(base + "_forcas.npz"):
        print(f"AVISO: Cordel sem resultados em {base}_forcas.npz -- rode o caso antes")
        return None
    F = np.load(base + "_forcas.npz")
    C = np.load(base + ".npz") if os.path.exists(base + ".npz") else None
    t, cd, cl = F["t"], F["cd"], F["cl"]
    d = dict(nome="Cordel (FEM)", cor=BLUE, t=t, cd=cd, cl=cl, **resumo_forcas(t, cd, cl, t[-1] / 2))
    if "cp_medio" in F.files:
        d.update(resumo_cp(F["theta"], F["cp_medio"]))
    if C is not None:
        d.update(x=C["x"], y=C["y"], itc=C["itc"], u=C["u"], v=C["v"],
                 custo=f"{int(C['nn'])} nós, {len(C['itc'])} triângulos, Δt ≈ {float(C['dt']):.2e}, "
                       f"{int(C['it'])} passos")
    log = glob.glob(os.path.join(pasta, "resultados", saida, "**", "*.log"), recursive=True)
    d["log"] = log[0] if log else None
    return d


def _ler_cp_csv(arq):
    d = np.loadtxt(arq, delimiter=",", skiprows=1)
    return d[:, 0], d[:, 2]


def _sem_campos(d, r, arq):
    print(f"AVISO: {arq} ausente: {d['nome']} entra so com Cp e o resumo do resultados.json")
    d.update(cd_med=r["Cd_mean"], cl_rms=r["Cl_rms"], cl_amp=r["Cl_amp"],
             st=np.nan if r.get("St") is None else r["St"])


def ler_fvm_o(pasta, nt, rext):
    """fvm-cfd, malha O (python/cylinder.py): historico e campos do campos.npz, Cp do CSV."""
    caso = os.path.join(pasta, "results", f"cilindro_O_Re100_NT{nt}_R{rext:g}")
    if not os.path.exists(os.path.join(caso, "cp_parede.csv")):
        print(f"AVISO: fvm-cfd (malha O) sem cp_parede.csv em {caso} -- rode python/cylinder.py --NT {nt}")
        return None
    r = json.load(open(os.path.join(caso, "resultados.json"), encoding="utf-8"))
    d = dict(nome=f"fvm-cfd, malha O {r['NR']}×{r['NT']}", cor=ORANGE, resumo=r,
             **resumo_cp(*_ler_cp_csv(os.path.join(caso, "cp_parede.csv"))),
             custo=f"{r['NR']}x{r['NT']} = {r['cells']} volumes, Δt = {r['dt']:.3g}, {r['steps']} passos, "
                   f"{r['seconds'] / 60:.0f} min")
    arq = os.path.join(caso, "campos.npz")
    if not os.path.exists(arq):
        _sem_campos(d, r, arq)
        return d
    D = np.load(arq)
    H = D["history"]
    t, cd, cl = H[:, 0] * float(D["dt"]), H[:, 1], H[:, 2]
    d.update(t=t, cd=cd, cl=cl, **resumo_forcas(t, cd, cl, t[-1] / 2))
    sys.path.insert(0, os.path.join(pasta, "python"))
    from polar import PolarGrid  # noqa: E402  (repositorio fvm-cfd)
    g = PolarGrid(R0=0.5, Rout=r["Rext"], NT=r["NT"], aspect=0.7)
    rut = g.rc[:, None] * D["ut"]
    w = ((rut[1:] - rut[:-1]) / g.drc[:, None]
         - (D["ur"][1:-1] - np.roll(D["ur"][1:-1], 1, 1)) / g.dth) / g.rf[1:-1, None]
    Rr, T = np.meshgrid(g.rf[1:-1], np.r_[g.thf, 2 * np.pi], indexing="ij")
    d.update(Xo=Rr * np.cos(T), Yo=Rr * np.sin(T), wo=np.c_[w, w[:, :1]])
    return d


def ler_fvm_mac(pasta, h):
    """fvm-cfd, malha uniforme + fronteira imersa (python/cylinder2d.py)."""
    caso = os.path.join(pasta, "results", f"cilindro_mac_Re100_h{h:g}")
    if not os.path.exists(os.path.join(caso, "cp_parede.csv")):
        print(f"AVISO: fvm-cfd (malha uniforme) sem resultados em {caso} -- rode python/cylinder2d.py")
        return None
    r = json.load(open(os.path.join(caso, "resultados.json"), encoding="utf-8"))
    d = dict(nome=f"fvm-cfd, uniforme D/{r['D_over_h']:g} + FI", cor=GREEN, resumo=r,
             **resumo_cp(*_ler_cp_csv(os.path.join(caso, "cp_parede.csv"))),
             custo=f"{r['nx']}x{r['ny']} = {r['cells']} células, Δt = {r['dt']:.3g}, {r['steps']} passos, "
                   f"{r['seconds'] / 60:.0f} min")
    arq = os.path.join(caso, "campos.npz")
    if not os.path.exists(arq):
        _sem_campos(d, r, arq)
        return d
    D = np.load(arq)
    t, cd, cl = D["t"], D["cd"], D["cl"]
    d.update(t=t, cd=cd, cl=cl, **resumo_forcas(t, cd, cl, r["t_mean"]), xf=D["xf"], yf=D["yf"], w=D["w"])
    return d


def estudos_fvm(pasta):
    """Resumos dos estudos de malha do fvm-cfd: malha O (R_ext = 40 D) e malha uniforme."""
    def ler(nome):
        arq = os.path.join(pasta, "results", nome, "resultados.json")
        return json.load(open(arq, encoding="utf-8")) if os.path.exists(arq) else []
    o = sorted([r for r in ler("estudo_cilindro_O_Re100") if r["Rext"] == 40], key=lambda r: r["cells"])
    m = sorted(ler("estudo_cilindro_mac_Re100"), key=lambda r: r["cells"])
    return o, m


def vorticidade_fem(x, y, itc, u, v):
    xe, ye = x[itc], y[itc]
    b = np.stack([ye[:, 1] - ye[:, 2], ye[:, 2] - ye[:, 0], ye[:, 0] - ye[:, 1]], 1)
    c = np.stack([xe[:, 2] - xe[:, 1], xe[:, 0] - xe[:, 2], xe[:, 1] - xe[:, 0]], 1)
    a2 = (xe[:, 1] - xe[:, 0]) * (ye[:, 2] - ye[:, 0]) - (xe[:, 2] - xe[:, 0]) * (ye[:, 1] - ye[:, 0])
    w_el = ((v[itc] * b).sum(1) - (u[itc] * c).sum(1)) / a2
    pa = np.abs(a2)
    return (np.bincount(itc.ravel(), np.repeat(w_el * pa, 3), len(x))
            / np.bincount(itc.ravel(), np.repeat(pa, 3), len(x)))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--cordel", default=os.path.join(PROJ, "cordel"))
    ap.add_argument("--saida-cordel", default="cylinder_trabalho")
    ap.add_argument("--fvm", default=os.path.join(PROJ, "fvm-cfd"))
    ap.add_argument("--h-fvm", type=float, default=0.04, help="h da malha uniforme (cylinder2d.py)")
    ap.add_argument("--nt-fvm", type=int, default=128, help="NT da malha O (cylinder.py)")
    ap.add_argument("--rext-fvm", type=float, default=40)
    ap.add_argument("--n-paineis", type=int, default=64)
    a = ap.parse_args()
    start("06_trabalho_cilindro")

    # --- escoamento potencial (paineis) --------------------------------------
    R, U = 0.5, 1.0
    p = Panels.from_points(*geo.circle(a.n_paineis, R))
    sol = solve(p, U, lifting=False)
    f = sol.forces(2 * R)
    pan = dict(nome=f"painéis (potencial, N = {a.n_paineis})", cor=GRAY,
               **resumo_cp(np.arctan2(p.yc, p.xc) % (2 * np.pi), sol.cp), cd_med=f["cd"], cl_med=f["cl"])
    print("== Escoamento potencial pelo metodo dos paineis ==")
    print(f"N = {a.n_paineis}: Cd = {f['cd']:.2e}, Cl = {f['cl']:.2e} (d'Alembert); "
          f"Cp estagnacao = {pan['cp_estag']:.4f}, Cp minimo = {pan['cp_min']:.4f} "
          f"em {pan['ang_cp_min']:.1f} graus, Cp base = {pan['cp_base']:.4f}")

    cod = [c for c in (ler_cordel(a.cordel, a.saida_cordel), ler_fvm_o(a.fvm, a.nt_fvm, a.rext_fvm),
                       ler_fvm_mac(a.fvm, a.h_fvm)) if c is not None]

    # --- tabela ---------------------------------------------------------------
    print("\n== Re = 100: medias na metade final de cada rodada ==")
    cols = [("Cd medio", "cd_med", "{:.3f}"), ("Cl medio", "cl_med", "{:+.4f}"), ("Cl' RMS", "cl_rms", "{:.3f}"),
            ("Cl' amplitude", "cl_amp", "{:.3f}"), ("Cd' RMS", "cd_rms", "{:.4f}"), ("St", "st", "{:.4f}"),
            ("Cp estagnacao", "cp_estag", "{:.3f}"), ("Cp minimo", "cp_min", "{:.3f}"),
            ("angulo do Cp min (graus)", "ang_cp_min", "{:.0f}"), ("Cp base (180 graus)", "cp_base", "{:.3f}"),
            ("Cd de pressao (int. Cp)", "cd_p", "{:.3f}")]
    todos = cod + [pan]
    print(f"{'':26s}" + "".join(f"{c['nome'][:30]:>32s}" for c in todos) + f"{'literatura':>16s}")
    for nome, k, fmt in cols:
        lit = {"cd_med": "1,33-1,35", "st": "0,164-0,165"}.get(k, "")
        print(f"{nome:26s}" + "".join(f"{fmt.format(c[k]) if k in c else '-':>32s}" for c in todos)
              + f"{lit:>16s}")
    for c in cod:
        if "custo" in c:
            print(f"custo {c['nome']}: {c['custo']}")
    print("Literatura (Re = 100, citada de memoria -- confira): St = 0,164-0,165 (Williamson 1996), "
          "Cd = 1,33-1,35 (Park, Kwon & Choi 1998; Henderson 1995).")

    # --- Cp -------------------------------------------------------------------
    fig, ax = plt.subplots(figsize=(8, 4.8))
    th = np.linspace(0, 180, 361)
    ax.plot(th, an.cylinder_cp(np.radians(180 - th), U, R), "--", color=GRAY, lw=1.2,
            label="potencial exato: 1 − 4 sen²θ")
    ax.plot(pan["ang"], pan["cp"], "o", ms=4, mfc="none", color=GRAY, label=pan["nome"])
    for c in cod:
        if "cp" in c:
            ax.plot(c["ang"], c["cp"], "-", lw=1.6, color=c["cor"], label=c["nome"] + ", média temporal")
    ax.axhline(0, color="k", lw=0.5)
    ax.set(xlabel="θ a partir do ponto de estagnação (graus)", ylabel="Cp", xlim=(0, 180), ylim=(-3.2, 1.3),
           title="Coeficiente de pressão na parede, Re = 100")
    ax.grid(alpha=0.25); ax.legend(frameon=False, fontsize=8.5, loc="lower left")
    fig.tight_layout(); fig.savefig(out("comparacao_cp.png"), dpi=140); plt.close(fig)

    # --- forcas e espectro ----------------------------------------------------
    serie = [c for c in cod if "t" in c]
    if serie:
        fig, axs = plt.subplots(3, 1, figsize=(10, 9))
        for c in serie:
            axs[0].plot(c["t"], c["cd"], color=c["cor"], lw=1.0, label=c["nome"])
            axs[1].plot(c["t"], c["cl"], color=c["cor"], lw=1.0, label=c["nome"])
            s = c["t"] >= c["t"][-1] / 2
            fr, F = espectro(c["t"][s], c["cl"][s])
            axs[2].plot(fr, F, color=c["cor"], lw=1.4, label=f"{c['nome']}: St = {c['st']:.3f}")
        for ax in axs[:2]:
            ax.axvspan(serie[0]["t"][-1] / 2, serie[0]["t"][-1], color="#e8f1f4", zorder=0)
        axs[0].set(ylabel="Cd", ylim=(0.8, 2.0),
                   title=f"Arrasto (faixa sombreada: janela das médias de {serie[0]['nome']})")
        axs[1].set(ylabel="Cl", xlabel="t U / D", ylim=(-0.6, 0.6), title="Sustentação")
        axs[2].axvspan(*LIT["St"], color="#cdd9de", label="literatura: St = 0,164–0,165")
        axs[2].set(xlim=(0, 0.6), xlabel="f D / U", ylabel="|FFT(Cl)| (normalizado)",
                   title="Espectro de Cl")
        for ax in axs:
            ax.grid(alpha=0.25); ax.legend(frameon=False, fontsize=8.5)
        fig.tight_layout(); fig.savefig(out("comparacao_forcas.png"), dpi=140); plt.close(fig)

    # --- vorticidade ----------------------------------------------------------
    viz = [c for c in cod if "w" in c or "itc" in c or "wo" in c]
    if viz:
        from matplotlib.tri import Triangulation
        fig, axs = plt.subplots(len(viz), 1, figsize=(11, 3.6 * len(viz)), squeeze=False)
        lim, lev = 3.0, np.linspace(-3.0, 3.0, 41)
        for ax, c in zip(axs[:, 0], viz):
            if "itc" in c:
                w = vorticidade_fem(c["x"], c["y"], c["itc"], c["u"], c["v"])
                cs = ax.tricontourf(Triangulation(c["x"], c["y"], c["itc"]), np.clip(w, -lim, lim),
                                    levels=lev, cmap="RdBu_r", extend="both")
            elif "wo" in c:
                cs = ax.contourf(c["Xo"], c["Yo"], np.clip(c["wo"], -lim, lim), levels=lev,
                                 cmap="RdBu_r", extend="both")
            else:
                cs = ax.contourf(c["xf"][1:-1], c["yf"][1:-1], np.clip(c["w"], -lim, lim).T,
                                 levels=lev, cmap="RdBu_r", extend="both")
            ax.add_patch(Circle((0, 0), R, fc="white", ec="#0b0b0b"))
            ax.set(xlim=(-2, 15), ylim=(-3, 3), aspect="equal", title=f"{c['nome']}: vorticidade em t = {c['t'][-1]:.0f}")
            fig.colorbar(cs, ax=ax, shrink=0.85, label="ω D/U")
        fig.tight_layout(); fig.savefig(out("comparacao_vorticidade.png"), dpi=130); plt.close(fig)

    # --- estudo de malha do fvm-cfd (malha O x malha uniforme) ----------------
    O, M = estudos_fvm(a.fvm)
    if O or M:
        print("\n== fvm-cfd: estudo de malha (Re = 100) ==")
        print(f"{'abordagem':>16} {'malha':>10} {'volumes':>8} {'Cd':>7} {'Cl amp':>7} {'St':>7} {'tempo':>11}")
        for r in O:
            m = f"{r['NR']}x{r['NT']}"
            print(f"{'malha O':>16} {m:>10} {r['cells']:8d} {r['Cd_mean']:7.3f} "
                  f"{r['Cl_amp']:7.3f} {r['St']:7.4f} {r['seconds'] / 60:7.1f} min")
        for r in M:
            m = f"D/{r['D_over_h']:g}"
            print(f"{'uniforme + FI':>16} {m:>10} {r['cells']:8d} {r['cd_mean']:7.3f} "
                  f"{r['cl_amp']:7.3f} {r['St']:7.4f} {r['seconds'] / 60:7.1f} min")
        fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
        for lista, kcd, cor, lab in ((O, "Cd_mean", ORANGE, "malha O (R_ext = 40 D)"),
                                     (M, "cd_mean", GREEN, "malha uniforme + FI (bloqueio 7 %)")):
            if lista:
                n = [r["cells"] for r in lista]
                axs[0].plot(n, [r[kcd] for r in lista], "o-", color=cor, label=lab)
                axs[1].plot(n, [r["St"] for r in lista], "o-", color=cor, label=lab)
        axs[0].axhspan(*LIT["Cd"], color="#cdd9de")
        axs[0].set(xscale="log", xlabel="volumes", ylabel="Cd médio", title="Arrasto × malha (faixa: literatura)")
        axs[1].axhspan(*LIT["St"], color="#cdd9de")
        axs[1].set(xscale="log", xlabel="volumes", ylabel="St", title="Strouhal × malha (faixa: literatura)")
        from matplotlib.ticker import NullFormatter
        for ax in axs:
            ax.xaxis.set_minor_formatter(NullFormatter())
            ax.grid(alpha=0.25, which="both")
            ax.legend(frameon=False, fontsize=8.5)
        fig.tight_layout(); fig.savefig(out("estudo_malha_fvm.png"), dpi=140); plt.close(fig)

    print("\nFiguras: " + ", ".join(rel(out(n)) for n in
                                     ("comparacao_cp.png", "comparacao_forcas.png", "comparacao_vorticidade.png")))


if __name__ == "__main__":
    main()
