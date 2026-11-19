"""Solvers de paineis 2D baseados exclusivamente em escoamento potencial.

* ``solve(..., lifting=False)`` -> metodo de fontes (Hess & Smith, 1967 para
  corpos sem sustentacao): N fontes constantes, condicao de impenetrabilidade.
* ``solve(..., lifting=True)``  -> metodo de Hess-Smith: N fontes constantes +
  um vortice de intensidade uniforme gamma em todos os paineis + condicao de
  Kutta (velocidades tangenciais iguais em modulo nos dois paineis adjacentes
  ao bordo de fuga, que e o primeiro/ultimo no do contorno).
* ``circulation=...`` impoe a circulacao (horaria) em vez da condicao de Kutta
  (ex.: cilindro girante / efeito Magnus).
* ``ground=h`` adiciona o efeito solo pelo metodo das imagens (parede em y = h).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .panels import Panels, influence


@dataclass
class Solution:
    panels: Panels
    U: float
    alpha: float
    sigma: np.ndarray
    gamma: float
    y_wall: float | None
    vt: np.ndarray = field(init=False)
    cp: np.ndarray = field(init=False)

    @property
    def circulation_cw(self):
        """Circulacao total no sentido horario (positiva gera sustentacao para +x)."""
        return -self.gamma * self.panels.perimeter

    # ------------------------------------------------------------------ campo
    def velocity(self, x, y):
        """Velocidade (u, v) em pontos arbitrarios do campo (fora do corpo)."""
        x = np.asarray(x, float)
        shape = x.shape
        Us, Vs, Uv, Vv = influence(self.panels, x.ravel(), np.asarray(y, float).ravel(),
                                   y_wall=self.y_wall)
        u = self.U * np.cos(self.alpha) + Us @ self.sigma + self.gamma * Uv.sum(axis=1)
        v = self.U * np.sin(self.alpha) + Vs @ self.sigma + self.gamma * Vv.sum(axis=1)
        return u.reshape(shape), v.reshape(shape)

    # ----------------------------------------------------------------- forcas
    def forces(self, ref_length=1.0, x_ref=0.25, y_ref=0.0):
        """Coeficientes por integracao da pressao nos paineis.

        Retorna dict com cl, cd (deve ser ~0: paradoxo de d'Alembert), cm
        (positivo cabrando, em torno de (x_ref, y_ref)), cx, cy e cl_kj
        (Kutta-Joukowski, valido sem efeito solo).
        """
        p = self.panels
        fx = -np.sum(self.cp * p.length * p.nx) / ref_length
        fy = -np.sum(self.cp * p.length * p.ny) / ref_length
        a = self.alpha
        cl = -fx * np.sin(a) + fy * np.cos(a)
        cd = fx * np.cos(a) + fy * np.sin(a)
        dfx = -self.cp * p.length * p.nx
        dfy = -self.cp * p.length * p.ny
        # momento cabrador (nariz para cima = horario) = -(r x F)_z
        cm = np.sum(-(p.xc - x_ref) * dfy + (p.yc - y_ref) * dfx) / ref_length**2
        cl_kj = 2.0 * self.circulation_cw / (self.U * ref_length)
        return dict(cl=cl, cd=cd, cm=cm, cx=fx, cy=fy, cl_kj=cl_kj)


def solve(panels: Panels, U=1.0, alpha=0.0, lifting=True, circulation=None, ground=None):
    """Resolve o escoamento potencial em torno do corpo descrito por ``panels``.

    Parametros
    ----------
    alpha : angulo de ataque em radianos (deve ser 0 com efeito solo; gire a
        geometria nesse caso, pois o escoamento livre e paralelo ao solo).
    lifting : True usa Hess-Smith com condicao de Kutta no primeiro/ultimo painel.
    circulation : circulacao horaria imposta (ignora a condicao de Kutta).
    ground : coordenada y do solo (metodo das imagens) ou None.
    """
    if ground is not None and abs(alpha) > 1e-12:
        raise ValueError("Com efeito solo use alpha=0 e gire a geometria.")
    N = panels.n
    Us, Vs, Uv, Vv = influence(panels, panels.xc, panels.yc, self_points=True, y_wall=ground)
    nx, ny, tx, ty = panels.nx[:, None], panels.ny[:, None], panels.tx[:, None], panels.ty[:, None]
    An = Us * nx + Vs * ny        # normal  devida as fontes
    At = Us * tx + Vs * ty        # tangencial devida as fontes
    Bn = (Uv * nx + Vv * ny).sum(axis=1)  # normal devida ao vortice uniforme
    Bt = (Uv * tx + Vv * ty).sum(axis=1)
    ux, uy = U * np.cos(alpha), U * np.sin(alpha)
    Vn_inf = ux * panels.nx + uy * panels.ny
    Vt_inf = ux * panels.tx + uy * panels.ty

    if circulation is not None:
        gamma = -circulation / panels.perimeter
        sigma = np.linalg.solve(An, -Vn_inf - gamma * Bn)
    elif lifting:
        A = np.zeros((N + 1, N + 1))
        b = np.zeros(N + 1)
        A[:N, :N] = An
        A[:N, N] = Bn
        b[:N] = -Vn_inf
        # Kutta: Vt(1) + Vt(N) = 0
        A[N, :N] = At[0] + At[-1]
        A[N, N] = Bt[0] + Bt[-1]
        b[N] = -(Vt_inf[0] + Vt_inf[-1])
        x = np.linalg.solve(A, b)
        sigma, gamma = x[:N], x[N]
    else:
        gamma = 0.0
        sigma = np.linalg.solve(An, -Vn_inf)

    sol = Solution(panels, U, alpha, sigma, gamma, ground)
    sol.vt = Vt_inf + At @ sigma + gamma * Bt
    sol.cp = 1.0 - (sol.vt / U) ** 2
    return sol


def rotate(x, y, angle, x0=0.0, y0=0.0):
    """Gira pontos de ``angle`` rad no sentido horario (nariz para cima com angle>0
    quando o corpo aponta para -x), em torno de (x0, y0). Usado para simular
    angulo de ataque no efeito solo."""
    c, s = np.cos(angle), np.sin(angle)
    dx, dy = np.asarray(x) - x0, np.asarray(y) - y0
    return x0 + c * dx + s * dy, y0 - s * dx + c * dy
