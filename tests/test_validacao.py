"""Testes de validacao contra solucoes analiticas da teoria potencial."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from panelpy import Panels, rotate, solve  # noqa: E402
from panelpy import analytic as an, geometry as geo  # noqa: E402
from panelpy.ahmed import ahmed_forces, ahmed_profile  # noqa: E402
from panelpy.elementary import (Doublet, Source, Uniform, Vortex,  # noqa: E402
                                cylinder_flow, rankine_flow)


# ------------------------------------------------------------ elementares
def test_cilindro_superposicao_psi_zero_na_superficie():
    th = np.linspace(0, 2 * np.pi, 361)
    psi = cylinder_flow(1.0, 1.5).stream_function(1.5 * np.cos(th), 1.5 * np.sin(th))
    assert np.abs(psi).max() < 1e-12


def test_velocidade_e_gradiente_do_potencial():
    flow = Uniform(1.3, 0.2) + Source(2.0, -1, 0.3) + Vortex(1.5, 0.5, -0.4) + Doublet(0.7, 0.2, 1.0)
    x, y, h = 2.1, -1.7, 1e-6
    u, v = flow.velocity(x, y)
    dphix = (flow.potential(x + h, y) - flow.potential(x - h, y)) / (2 * h)
    dphiy = (flow.potential(x, y + h) - flow.potential(x, y - h)) / (2 * h)
    dpsix = (flow.stream_function(x + h, y) - flow.stream_function(x - h, y)) / (2 * h)
    dpsiy = (flow.stream_function(x, y + h) - flow.stream_function(x, y - h)) / (2 * h)
    assert np.allclose([u, v], [dphix, dphiy], atol=1e-6)
    assert np.allclose([u, v], [dpsiy, -dpsix], atol=1e-6)


def test_circulacao_do_vortice():
    G = 3.0
    th = np.linspace(0, 2 * np.pi, 2001)[:-1]
    u, v = Vortex(G).velocity(np.cos(th), np.sin(th))
    circ = np.sum(-u * np.sin(th) + v * np.cos(th)) * (2 * np.pi / th.size)
    assert circ == pytest.approx(G, rel=1e-10)


def test_rankine_pontos_de_estagnacao():
    U, m, a = 1.0, 2 * np.pi, 1.0
    l = an.rankine_half_length(U, m, a)
    u, v = rankine_flow(U, m, a).velocity(np.array([-l, l]), np.zeros(2))
    assert np.abs(u).max() < 1e-12
    h = an.rankine_half_thickness(U, m, a)
    assert rankine_flow(U, m, a).stream_function(0.0, h) == pytest.approx(0.0, abs=1e-12)


@pytest.mark.parametrize("G", [0.0, 2.0, 4 * np.pi])
def test_magnus_integracao_de_pressao(G):
    th = np.linspace(0, 2 * np.pi, 2001)[:-1]
    cp = an.cylinder_cp(th, 1.0, 1.0, G)
    cl = -np.sum(cp * np.sin(th)) * (2 * np.pi / th.size) / 2
    assert cl == pytest.approx(an.cylinder_cl(1.0, 1.0, G), abs=1e-10)


# ------------------------------------------------------------ paineis
@pytest.mark.parametrize("N", [16, 64, 128])
def test_cilindro_paineis_cp(N):
    p = Panels.from_points(*geo.circle(N))
    s = solve(p, lifting=False)
    assert np.abs(s.cp - an.cylinder_cp(np.arctan2(p.yc, p.xc))).max() < 1e-10
    assert abs(np.sum(s.sigma * p.length)) < 1e-12  # corpo fechado: vazao liquida nula
    assert abs(s.forces(2.0)["cd"]) < 1e-12        # d'Alembert


def test_cilindro_girante_paineis():
    p = Panels.from_points(*geo.circle(256))
    f = solve(p, circulation=2 * np.pi).forces(2.0)
    assert f["cl"] == pytest.approx(an.cylinder_cl(1.0, 1.0, 2 * np.pi), rel=5e-3)
    assert abs(f["cd"]) < 1e-10


def test_rankine_paineis():
    U, m, a = 1.0, 2 * np.pi, 1.0
    p = Panels.from_points(*geo.rankine_oval(256, U, m, a))
    s = solve(p, U, lifting=False)
    err = np.abs(s.cp - an.rankine_cp(np.arctan2(p.yc, p.xc), U, m, a)).max()
    assert err < 1e-2


def test_karman_trefftz_cl():
    kt = an.Joukowski(0.1, 0.1, tau_deg=10.0)
    a = np.radians(5)
    f = solve(Panels.from_points(*kt.coordinates(400)), alpha=a).forces()
    assert f["cl"] == pytest.approx(kt.cl(a), rel=1e-2)
    assert f["cl_kj"] == pytest.approx(kt.cl(a), rel=2e-3)
    assert abs(f["cd"]) < 2e-3


def test_joukowski_tau0_equivale_a_transformacao_classica():
    j = an.Joukowski(0.1, 0.05)
    zeta = j._zeta(np.linspace(0.1, 6.0, 50))
    assert np.allclose(j._map(zeta), zeta + 1 / zeta)


def test_naca_simetrico_e_aerofolio_fino():
    p = Panels.from_points(*geo.naca4("0012", 200))
    assert abs(solve(p, alpha=0.0).forces()["cl"]) < 1e-10
    p = Panels.from_points(*geo.naca4("2412", 200))
    al = np.radians([0.0, 4.0])
    cl = [solve(p, alpha=a).forces()["cl"] for a in al]
    slope = (cl[1] - cl[0]) / (al[1] - al[0])
    a0 = -cl[0] / slope
    a0_taf = an.thin_airfoil_alpha0(geo.naca4_camber("2412")[1])
    assert np.degrees(abs(a0 - a0_taf)) < 0.15
    assert 2 * np.pi < slope < 1.15 * 2 * np.pi  # espessura aumenta ligeiramente dCl/da


def test_efeito_solo_parede_impermeavel():
    x, y = geo.naca4("4412", 120)
    x, y = rotate(x, y, np.radians(4), 0.25, 0.0)
    s = solve(Panels.from_points(x, y + 0.3), ground=0.0)
    xs = np.linspace(-3, 4, 71)
    _, v = s.velocity(xs, np.zeros_like(xs))
    assert np.abs(v).max() < 1e-12


# ------------------------------------------------------------ Ahmed
def test_ahmed_fechado_dalembert():
    x, y, info = ahmed_profile(25, 1200, separation=None)
    f = solve(Panels.from_points(x, y), lifting=False).forces(1.044)
    assert abs(f["cl"]) < 1e-2 and abs(f["cd"]) < 2e-3


def test_ahmed_esteira_convergencia_malha():
    cls = []
    for n in (400, 800):
        x, y, info = ahmed_profile(25, n, separation="slant_end")
        cls.append(ahmed_forces(solve(Panels.from_points(x, y), ground=-0.05), info)["cl_body"])
    assert cls[0] == pytest.approx(cls[1], rel=5e-3)
