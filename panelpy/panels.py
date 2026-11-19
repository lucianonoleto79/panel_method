"""Discretizacao em paineis retos e coeficientes de influencia analiticos.

Cada painel j vai do no (x_j, y_j) ao no (x_{j+1}, y_{j+1}). Os nos devem
percorrer o contorno no sentido HORARIO (convencao de Hess-Smith/Anderson),
de modo que a normal n = (-t_y, t_x) aponta para FORA do corpo.

Para um painel de comprimento L, no sistema local (xi ao longo do painel,
eta na direcao da normal), uma distribuicao de fontes de intensidade
unitaria e constante induz (Katz & Plotkin, cap. 10):

    u_xi  = (1 / 2 pi) ln(r1 / r2)
    u_eta = (1 / 2 pi) (theta2 - theta1)

com r1, r2 as distancias as extremidades e theta_k = atan2(eta, xi - xi_k).
Um vortice de intensidade unitaria (anti-horario) induz a mesma velocidade
girada de +90 graus:  (u_xi, u_eta)_vortice = (-u_eta, u_xi)_fonte.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Panels:
    """Geometria dos paineis (arrays de tamanho N)."""

    xa: np.ndarray  # no inicial
    ya: np.ndarray
    xb: np.ndarray  # no final
    yb: np.ndarray

    def __post_init__(self):
        dx, dy = self.xb - self.xa, self.yb - self.ya
        self.length = np.hypot(dx, dy)
        if np.any(self.length <= 0):
            raise ValueError("Ha paineis de comprimento nulo (nos repetidos).")
        self.tx, self.ty = dx / self.length, dy / self.length
        self.nx, self.ny = -self.ty, self.tx
        self.xc = 0.5 * (self.xa + self.xb)  # pontos de controle (colocacao)
        self.yc = 0.5 * (self.ya + self.yb)

    @classmethod
    def from_points(cls, x, y, check_orientation=True):
        """Cria paineis a partir de nos (o contorno e fechado automaticamente)."""
        x, y = np.asarray(x, float), np.asarray(y, float)
        if np.hypot(x[0] - x[-1], y[0] - y[-1]) > 1e-12:
            x, y = np.append(x, x[0]), np.append(y, y[0])
        if check_orientation:
            area2 = np.sum(x[:-1] * y[1:] - x[1:] * y[:-1])  # > 0 se anti-horario
            if area2 > 0:
                raise ValueError("Os nos devem estar em sentido horario.")
        return cls(x[:-1], y[:-1], x[1:], y[1:])

    @property
    def n(self):
        return self.xa.size

    @property
    def perimeter(self):
        return self.length.sum()

    def mirrored(self, y_wall):
        """Paineis espelhados em relacao a reta y = y_wall (para ilustracao)."""
        return Panels(self.xb, 2 * y_wall - self.yb, self.xa, 2 * y_wall - self.ya)

    def polygon(self):
        return np.append(self.xa, self.xa[0]), np.append(self.ya, self.ya[0])


def _source_local(panels: Panels, X, Y):
    """Velocidade local (u_xi, u_eta) de fontes unitarias. X, Y: (M,1); retorno (M,N)."""
    dx = X - panels.xa[None, :]
    dy = Y - panels.ya[None, :]
    xi = dx * panels.tx + dy * panels.ty
    eta = -dx * panels.ty + dy * panels.tx
    L = panels.length[None, :]
    r1sq = xi**2 + eta**2
    r2sq = (xi - L) ** 2 + eta**2
    u_xi = np.log(r1sq / r2sq) / (4.0 * np.pi)
    u_eta = (np.arctan2(eta, xi - L) - np.arctan2(eta, xi)) / (2.0 * np.pi)
    return u_xi, u_eta


def influence(panels: Panels, X, Y, self_points=False, y_wall=None):
    """Velocidades induzidas por fontes e vortices unitarios de cada painel.

    Parametros
    ----------
    X, Y : arrays 1D com os M pontos de avaliacao.
    self_points : True quando (X, Y) sao os pontos de controle dos proprios
        paineis; aplica o limite pelo lado externo (u_eta = 1/2 para fonte).
    y_wall : se dado, inclui as imagens espelhadas em y = y_wall (efeito solo).

    Retorna
    -------
    (Us, Vs, Uv, Vv) : matrizes (M, N) com as componentes globais (x, y) da
        velocidade devida a fontes unitarias (Us, Vs) e vortices unitarios
        anti-horarios (Uv, Vv) em cada painel.
    """
    X = np.asarray(X, float).reshape(-1, 1)
    Y = np.asarray(Y, float).reshape(-1, 1)
    us_xi, us_eta = _source_local(panels, X, Y)
    if self_points:
        idx = np.arange(panels.n)
        us_xi[idx, idx] = 0.0
        us_eta[idx, idx] = 0.5
    tx, ty, nx, ny = panels.tx, panels.ty, panels.nx, panels.ny
    Us = us_xi * tx + us_eta * nx
    Vs = us_xi * ty + us_eta * ny
    uv_xi, uv_eta = -us_eta, us_xi
    Uv = uv_xi * tx + uv_eta * nx
    Vv = uv_xi * ty + uv_eta * ny
    if y_wall is not None:
        # Imagem: V_img(P) = M V(P'), com P' o reflexo de P e M = diag(1, -1).
        # Fonte-imagem tem a mesma intensidade e vortice-imagem o sinal oposto;
        # ambos sao obtidos automaticamente por essa reflexao do campo.
        Yr = 2.0 * y_wall - Y
        si_xi, si_eta = _source_local(panels, X, Yr)
        Us += si_xi * tx + si_eta * nx
        Vs -= si_xi * ty + si_eta * ny
        vi_xi, vi_eta = -si_eta, si_xi
        Uv += vi_xi * tx + vi_eta * nx
        Vv -= vi_xi * ty + vi_eta * ny
    return Us, Vs, Uv, Vv
