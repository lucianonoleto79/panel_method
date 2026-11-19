"""Solucoes analiticas de referencia (teoria potencial) usadas na validacao."""
from __future__ import annotations

import numpy as np

from .elementary import rankine_flow


# --------------------------------------------------------------------------
# Cilindro circular
# --------------------------------------------------------------------------
def cylinder_cp(theta, U=1.0, R=1.0, gamma_cw=0.0):
    """Cp exato na superficie do cilindro (theta medido a partir de +x, anti-horario).

    V_theta = -2 U sin(theta) - Gamma_cw / (2 pi R)
    """
    vt = -2.0 * U * np.sin(theta) - gamma_cw / (2.0 * np.pi * R)
    return 1.0 - (vt / U) ** 2


def cylinder_cl(U=1.0, R=1.0, gamma_cw=0.0):
    """Cl por unidade de diametro (Kutta-Joukowski): L'/(q 2R) = Gamma/(U R)."""
    return gamma_cw / (U * R)


# --------------------------------------------------------------------------
# Oval de Rankine
# --------------------------------------------------------------------------
def rankine_half_length(U=1.0, m=1.0, a=1.0):
    """Semi-comprimento l = a sqrt(1 + m / (pi U a)) (pontos de estagnacao +-l)."""
    return a * np.sqrt(1.0 + m / (np.pi * U * a))


def rankine_half_thickness(U=1.0, m=1.0, a=1.0, tol=1e-14):
    """Semi-espessura h: raiz de h = m/(2 pi U) (pi - 2 atan(h/a)) (bissecao)."""
    f = lambda h: h - m / (2.0 * np.pi * U) * (np.pi - 2.0 * np.arctan(h / a))
    lo, hi = 0.0, m / (2.0 * U) + a
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)


def rankine_contour_radius(phi, U=1.0, m=1.0, a=1.0, tol=1e-13):
    """Raio r(phi) da linha psi = 0 (contorno do corpo) ao longo do raio de angulo phi.

    Usa bissecao: dentro do corpo psi tem sinal oposto ao de fora.
    """
    flow = rankine_flow(U, m, a)
    phi = np.atleast_1d(np.asarray(phi, dtype=float))
    out = np.empty_like(phi)
    l = rankine_half_length(U, m, a)
    for k, p in enumerate(phi):
        s = np.sin(p)
        if abs(s) < 1e-12:
            out[k] = l
            continue
        sign = 1.0 if s > 0 else -1.0  # simetria: psi(x,-y) = -psi(x,y)
        lo, hi = 1e-9, 2.0 * l
        while hi - lo > tol:
            mid = 0.5 * (lo + hi)
            psi = sign * flow.stream_function(mid * np.cos(p), mid * s)
            lo, hi = (mid, hi) if psi < 0 else (lo, mid)
        out[k] = 0.5 * (lo + hi)
    return out


def rankine_cp(phi, U=1.0, m=1.0, a=1.0):
    """Cp exato na superficie do oval de Rankine no angulo polar phi."""
    r = rankine_contour_radius(phi, U, m, a)
    u, v = rankine_flow(U, m, a).velocity(r * np.cos(phi), r * np.sin(phi))
    return 1.0 - (u**2 + v**2) / U**2


# --------------------------------------------------------------------------
# Aerofolios de Joukowski e Karman-Trefftz (transformacao conforme)
# --------------------------------------------------------------------------
class Joukowski:
    """Aerofolio de Karman-Trefftz exato (Joukowski quando tau = 0).

    Transformacao:  z = k (1 + w) / (1 - w),  w = ((zeta - 1)/(zeta + 1))^k,
    com k = 2 - tau/pi  (tau = angulo do bordo de fuga). Para k = 2 ela se
    reduz a z = zeta + 1/zeta (Joukowski, bordo de fuga em cuspide).

    O circulo no plano zeta tem centro mu = (-mx, my) e passa por zeta = 1
    (bordo de fuga). mx controla a espessura e my o arqueamento.
    """

    def __init__(self, mx=0.1, my=0.0, tau_deg=0.0):
        self.mu = complex(-mx, my)
        self.a = abs(1.0 - self.mu)
        self.beta = np.arcsin(my / self.a)
        self.k = 2.0 - np.radians(tau_deg) / np.pi
        t = np.linspace(0.0, 2.0 * np.pi, 20001)[1:-1]
        z = self._map(self._zeta(t))
        self.x_te = self.k
        self.x_le = z.real.min()
        self.chord = self.x_te - self.x_le

    def _zeta(self, t):
        """Pontos do circulo; t em (0, 2 pi) percorre no sentido horario a partir do BF."""
        return self.mu + self.a * np.exp(1j * (-self.beta - t))

    def _w(self, zeta):
        r = (zeta - 1.0) / (zeta + 1.0)
        # ramo continuo ao longo do contorno (evita o corte do ramo principal)
        ang = np.unwrap(np.angle(r))
        return np.abs(r) ** self.k * np.exp(1j * self.k * ang)

    def _map(self, zeta):
        w = self._w(zeta)
        return self.k * (1.0 + w) / (1.0 - w)

    def _dzdzeta(self, zeta):
        w = self._w(zeta)
        return 4.0 * self.k**2 * w / ((1.0 - w) ** 2 * (zeta**2 - 1.0))

    def coordinates(self, n=200, dense=20000):
        """n paineis, sentido horario a partir do BF, corda 1, BA em x = 0.

        Os nos sao redistribuidos com espacamento cosseno no comprimento de arco
        de cada superficie (como nos perfis NACA)."""
        t = np.linspace(0.0, 2.0 * np.pi, dense + 1)[1:-1]
        z = self._map(self._zeta(t))
        x = np.r_[self.x_te, z.real, self.x_te]
        y = np.r_[0.0, z.imag, 0.0]
        x, y = cosine_resample(x, y, n)
        x = (x - self.x_le) / self.chord
        y = y / self.chord
        x[0] = x[-1] = 1.0
        y[0] = y[-1] = 0.0
        return x, y

    def circulation(self, alpha, U=1.0):
        """Circulacao horaria imposta pela condicao de Kutta: 4 pi U a sin(alpha + beta)."""
        return 4.0 * np.pi * U * self.a * np.sin(alpha + self.beta)

    def cl(self, alpha):
        """Cl exato (Kutta-Joukowski) baseado na corda."""
        return 2.0 * self.circulation(alpha) / self.chord

    def surface_cp(self, alpha, n=400, U=1.0):
        """(x/c, Cp) exatos ao longo da superficie (exclui o bordo de fuga)."""
        t = np.linspace(0.0, 2.0 * np.pi, n + 1)[1:-1]
        zeta = self._zeta(t)
        zp = zeta - self.mu
        G = self.circulation(alpha, U)
        dw = U * (np.exp(-1j * alpha) - self.a**2 * np.exp(1j * alpha) / zp**2) \
            + 1j * G / (2.0 * np.pi * zp)
        V = np.abs(dw / self._dzdzeta(zeta))
        z = self._map(zeta)
        return (z.real - self.x_le) / self.chord, 1.0 - (V / U) ** 2


def cosine_resample(x, y, n):
    """Redistribui um contorno fechado (BF -> intradorso -> BA -> extradorso -> BF)
    em n paineis com espacamento cosseno no comprimento de arco de cada lado."""
    ile = int(np.argmin(x))
    out = []
    for xs, ys in ((x[: ile + 1], y[: ile + 1]), (x[ile:], y[ile:])):
        s = np.r_[0.0, np.cumsum(np.hypot(np.diff(xs), np.diff(ys)))]
        si = s[-1] * 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, n // 2 + 1)))
        out.append((np.interp(si, s, xs), np.interp(si, s, ys)))
    return (np.r_[out[0][0], out[1][0][1:]], np.r_[out[0][1], out[1][1][1:]])


# --------------------------------------------------------------------------
# Teoria do aerofolio fino
# --------------------------------------------------------------------------
def thin_airfoil_alpha0(camber_slope, n=4000):
    """Angulo de sustentacao nula (rad) pela teoria do aerofolio fino.

    alpha_L0 = -(1/pi) * int_0^pi dyc/dx (cos(theta) - 1) d(theta),  x = (1 - cos theta)/2
    ``camber_slope`` e uma funcao dyc/dx(x) com x em [0, 1] (corda unitaria).
    """
    th = (np.arange(n) + 0.5) * np.pi / n  # ponto medio
    x = 0.5 * (1.0 - np.cos(th))
    return -np.sum(camber_slope(x) * (np.cos(th) - 1.0)) * (np.pi / n) / np.pi


def thin_airfoil_cl(alpha, alpha0=0.0):
    """Cl = 2 pi (alpha - alpha_L0)."""
    return 2.0 * np.pi * (alpha - alpha0)
