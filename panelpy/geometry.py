"""Geradores de contorno. Todos retornam (x, y) fechados (primeiro no == ultimo),
no sentido HORARIO, comecando no "bordo de fuga" (onde a condicao de Kutta
e aplicada pelo solver de Hess-Smith)."""
from __future__ import annotations

import numpy as np

from .analytic import Joukowski, rankine_contour_radius


def circle(n=64, R=1.0, xc=0.0, yc=0.0):
    """Circulo com n paineis, comecando em (R, 0) e descendo (horario)."""
    th = np.linspace(0.0, -2.0 * np.pi, n + 1)
    x, y = xc + R * np.cos(th), yc + R * np.sin(th)
    x[-1], y[-1] = x[0], y[0]
    return x, y


def rankine_oval(n=80, U=1.0, m=1.0, a=1.0):
    """Contorno psi = 0 do oval de Rankine, com espacamento uniforme no angulo polar."""
    phi = np.linspace(0.0, -2.0 * np.pi, n + 1)
    r = rankine_contour_radius(phi, U, m, a)
    x, y = r * np.cos(phi), r * np.sin(phi)
    x[-1], y[-1] = x[0], y[0]
    return x, y


# --------------------------------------------------------------------------
# NACA 4 digitos
# --------------------------------------------------------------------------
def naca4_camber(code: str):
    """Retorna (yc(x), dyc/dx(x)) da linha media de um perfil NACA 4 digitos."""
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0

    def yc(x):
        x = np.asarray(x, float)
        if m == 0:
            return np.zeros_like(x)
        return np.where(x < p, m / p**2 * (2 * p * x - x**2),
                        m / (1 - p) ** 2 * ((1 - 2 * p) + 2 * p * x - x**2))

    def dyc(x):
        x = np.asarray(x, float)
        if m == 0:
            return np.zeros_like(x)
        return np.where(x < p, 2 * m / p**2 * (p - x), 2 * m / (1 - p) ** 2 * (p - x))

    return yc, dyc


def naca4(code="0012", n=160, closed_te=True):
    """Perfil NACA 4 digitos com n paineis (n par), espacamento cosseno, corda 1.

    closed_te=True usa o coeficiente -0.1036 (bordo de fuga fechado), adequado a
    condicao de Kutta; com False usa o -0.1015 original (bordo de fuga aberto,
    que e entao fechado por um painel curto).
    """
    if n % 2:
        raise ValueError("n deve ser par")
    t = int(code[2:]) / 100.0
    a4 = -0.1036 if closed_te else -0.1015
    beta = np.linspace(0.0, np.pi, n // 2 + 1)
    xs = 0.5 * (1.0 - np.cos(beta))  # 0 (BA) -> 1 (BF)
    yt = 5 * t * (0.2969 * np.sqrt(xs) - 0.1260 * xs - 0.3516 * xs**2
                  + 0.2843 * xs**3 + a4 * xs**4)
    yc, dyc = naca4_camber(code)
    th = np.arctan(dyc(xs))
    xu, yu = xs - yt * np.sin(th), yc(xs) + yt * np.cos(th)
    xl, yl = xs + yt * np.sin(th), yc(xs) - yt * np.cos(th)
    # horario: BF -> intradorso -> BA -> extradorso -> BF
    x = np.concatenate([xl[::-1], xu[1:]])
    y = np.concatenate([yl[::-1], yu[1:]])
    if closed_te:
        x[0] = x[-1] = 1.0
        y[0] = y[-1] = 0.5 * (y[0] + y[-1])
    return x, y


def joukowski(n=160, mx=0.1, my=0.0, tau_deg=0.0):
    """Aerofolio de Joukowski (tau=0) ou Karman-Trefftz, corda 1, BA em x = 0."""
    return Joukowski(mx, my, tau_deg).coordinates(n)


# --------------------------------------------------------------------------
# Utilitarios para montar contornos por segmentos
# --------------------------------------------------------------------------
def _line(p0, p1, n, cluster=True):
    """n pontos de p0 (incluso) a p1 (excluso), com agrupamento cosseno nas pontas."""
    s = np.linspace(0.0, 1.0, n + 1)
    if cluster:
        s = 0.5 * (1.0 - np.cos(np.pi * s))
    p0, p1 = np.asarray(p0, float), np.asarray(p1, float)
    return (p0[None, :] + s[:, None] * (p1 - p0)[None, :])[:-1]


def _arc(center, R, a0, a1, n):
    """n pontos num arco de a0 (incluso) a a1 (excluso), angulos em rad."""
    a = np.linspace(a0, a1, n + 1)[:-1]
    return np.column_stack([center[0] + R * np.cos(a), center[1] + R * np.sin(a)])
