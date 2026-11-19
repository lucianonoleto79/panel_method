"""Escoamentos potenciais elementares 2D (solucoes analiticas exatas).

Convencoes
----------
* Velocidades: u = d(phi)/dx = d(psi)/dy ,  v = d(phi)/dy = -d(psi)/dx
* Circulacao ``gamma`` dos vortices: positiva no sentido ANTI-horario.
* Todas as funcoes aceitam escalares ou arrays numpy (broadcast).

Os escoamentos podem ser somados (superposicao), pois a equacao de Laplace
e linear:  ``fluxo = Uniforme(1) + Dipolo(2*pi)``  -> cilindro de raio 1.
"""
from __future__ import annotations

import numpy as np

TWO_PI = 2.0 * np.pi


class Flow:
    """Classe-base: um escoamento potencial com phi, psi e velocidade."""

    def velocity(self, x, y):
        raise NotImplementedError

    def potential(self, x, y):
        raise NotImplementedError

    def stream_function(self, x, y):
        raise NotImplementedError

    def __add__(self, other: "Flow") -> "Superposition":
        return Superposition([self, other])


class Superposition(Flow):
    """Soma de varios escoamentos elementares."""

    def __init__(self, flows):
        self.flows = []
        for f in flows:
            self.flows.extend(f.flows if isinstance(f, Superposition) else [f])

    def velocity(self, x, y):
        u = np.zeros(np.broadcast(x, y).shape)
        v = np.zeros_like(u)
        for f in self.flows:
            du, dv = f.velocity(x, y)
            u = u + du
            v = v + dv
        return u, v

    def potential(self, x, y):
        return sum(f.potential(x, y) for f in self.flows)

    def stream_function(self, x, y):
        return sum(f.stream_function(x, y) for f in self.flows)

    def __add__(self, other):
        return Superposition([self, other])


class Uniform(Flow):
    """Escoamento uniforme de modulo U e angulo de ataque alpha (rad)."""

    def __init__(self, U=1.0, alpha=0.0):
        self.U, self.alpha = U, alpha

    def velocity(self, x, y):
        shape = np.broadcast(x, y).shape
        return (np.full(shape, self.U * np.cos(self.alpha)),
                np.full(shape, self.U * np.sin(self.alpha)))

    def potential(self, x, y):
        return self.U * (x * np.cos(self.alpha) + y * np.sin(self.alpha))

    def stream_function(self, x, y):
        return self.U * (y * np.cos(self.alpha) - x * np.sin(self.alpha))


class Source(Flow):
    """Fonte (m > 0) ou sumidouro (m < 0) de vazao m por unidade de envergadura."""

    def __init__(self, m, x0=0.0, y0=0.0):
        self.m, self.x0, self.y0 = m, x0, y0

    def velocity(self, x, y):
        dx, dy = x - self.x0, y - self.y0
        r2 = dx**2 + dy**2
        return self.m / TWO_PI * dx / r2, self.m / TWO_PI * dy / r2

    def potential(self, x, y):
        return self.m / TWO_PI * 0.5 * np.log((x - self.x0) ** 2 + (y - self.y0) ** 2)

    def stream_function(self, x, y):
        return self.m / TWO_PI * np.arctan2(y - self.y0, x - self.x0)


class Vortex(Flow):
    """Vortice livre (potencial) de circulacao gamma (anti-horario > 0)."""

    def __init__(self, gamma, x0=0.0, y0=0.0):
        self.gamma, self.x0, self.y0 = gamma, x0, y0

    def velocity(self, x, y):
        dx, dy = x - self.x0, y - self.y0
        r2 = dx**2 + dy**2
        return -self.gamma / TWO_PI * dy / r2, self.gamma / TWO_PI * dx / r2

    def potential(self, x, y):
        return self.gamma / TWO_PI * np.arctan2(y - self.y0, x - self.x0)

    def stream_function(self, x, y):
        return -self.gamma / TWO_PI * 0.5 * np.log((x - self.x0) ** 2 + (y - self.y0) ** 2)


class Doublet(Flow):
    """Dipolo de intensidade kappa com eixo na direcao -x (convencao de Anderson).

    phi = kappa cos(theta) / (2 pi r) ,  psi = -kappa sin(theta) / (2 pi r)
    Somado a um uniforme U, gera o cilindro de raio R = sqrt(kappa / (2 pi U)).
    """

    def __init__(self, kappa, x0=0.0, y0=0.0):
        self.kappa, self.x0, self.y0 = kappa, x0, y0

    def velocity(self, x, y):
        dx, dy = x - self.x0, y - self.y0
        r4 = (dx**2 + dy**2) ** 2
        c = self.kappa / TWO_PI
        return c * (dy**2 - dx**2) / r4, -c * 2.0 * dx * dy / r4

    def potential(self, x, y):
        dx, dy = x - self.x0, y - self.y0
        return self.kappa / TWO_PI * dx / (dx**2 + dy**2)

    def stream_function(self, x, y):
        dx, dy = x - self.x0, y - self.y0
        return -self.kappa / TWO_PI * dy / (dx**2 + dy**2)


# --------------------------------------------------------------------------
# Composicoes classicas
# --------------------------------------------------------------------------
def cylinder_flow(U=1.0, R=1.0, gamma_cw=0.0) -> Superposition:
    """Cilindro (uniforme + dipolo) com circulacao horaria opcional gamma_cw."""
    flow = Uniform(U) + Doublet(TWO_PI * U * R**2)
    if gamma_cw != 0.0:
        flow = flow + Vortex(-gamma_cw)
    return flow


def rankine_flow(U=1.0, m=1.0, a=1.0) -> Superposition:
    """Corpo (oval) de Rankine: fonte em (-a,0), sumidouro em (+a,0) e uniforme U."""
    return Uniform(U) + Source(m, -a, 0.0) + Source(-m, a, 0.0)


def pressure_coefficient(u, v, U=1.0):
    """Cp = 1 - |V|^2 / U^2 (Bernoulli, escoamento incompressivel)."""
    return 1.0 - (u**2 + v**2) / U**2
