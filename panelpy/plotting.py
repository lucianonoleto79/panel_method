"""Funcoes auxiliares de visualizacao (matplotlib)."""
from __future__ import annotations

import numpy as np
from matplotlib.path import Path

BLUE, ORANGE, GRAY = "#1f5fa8", "#d9730d", "#555555"


def field_grid(xlim, ylim, n=200):
    x = np.linspace(*xlim, n)
    y = np.linspace(*ylim, int(n * (ylim[1] - ylim[0]) / (xlim[1] - xlim[0])) or n)
    return np.meshgrid(x, y)


def mask_inside(X, Y, xb, yb, pad=0.0):
    """Mascara booleana dos pontos dentro do poligono (xb, yb)."""
    path = Path(np.column_stack([xb, yb]))
    pts = np.column_stack([X.ravel(), Y.ravel()])
    return path.contains_points(pts, radius=pad).reshape(X.shape) | \
        path.contains_points(pts, radius=-pad).reshape(X.shape)


def streamlines(ax, velocity, xlim, ylim, body=None, ground=None, n=220, U=1.0,
                cp_levels=np.linspace(-2.0, 1.0, 31), density=1.6, colorbar=True):
    """Linhas de corrente sobre o mapa de Cp para uma funcao velocity(X, Y)."""
    X, Y = field_grid(xlim, ylim, n)
    u, v = velocity(X, Y)
    bad = np.zeros_like(X, bool)
    if body is not None:
        bad |= mask_inside(X, Y, body[0], body[1], pad=1e-3 * (xlim[1] - xlim[0]))
    if ground is not None:
        bad |= Y < ground
    u = np.ma.array(u, mask=bad)
    v = np.ma.array(v, mask=bad)
    cp = 1.0 - (u**2 + v**2) / U**2
    cf = ax.contourf(X, Y, np.clip(cp, cp_levels[0], cp_levels[-1]), levels=cp_levels,
                     cmap="RdBu_r", extend="both")
    ax.streamplot(X, Y, u.filled(np.nan), v.filled(np.nan), color="k", linewidth=0.6,
                  density=density, arrowstyle="-", broken_streamlines=False)
    if body is not None:
        ax.fill(body[0], body[1], color="0.85", ec="k", lw=1.0, zorder=3)
    if ground is not None:
        ax.axhspan(ylim[0], ground, color="0.6", zorder=3)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal")
    if colorbar:
        ax.figure.colorbar(cf, ax=ax, shrink=0.8, label="$C_p$")
    return cf


def cp_axes(ax, xlabel="x/c"):
    ax.invert_yaxis()
    ax.set_xlabel(xlabel)
    ax.set_ylabel("$C_p$")
    ax.grid(alpha=0.3)
