"""Corpo de Ahmed em 2D (perfil lateral, plano de simetria) por escoamento potencial.

AVISO: o corpo de Ahmed e um corpo rombudo 3D cujo escoamento e dominado por
separacao, vortices longitudinais nas arestas laterais da rampa e esteira.
Nada disso existe em escoamento potencial 2D. Aqui o objetivo e apenas obter
estimativas de sustentacao sob hipoteses explicitas:

1. ``separation=None`` -> corpo fechado. Sem circulacao, a forca e nula fora
   do solo (d'Alembert) e o efeito solo gera apenas succao para baixo. Com a
   condicao de Kutta em uma unica aresta traseira o resultado depende
   completamente da aresta escolhida (mostrado nos exemplos como contraexemplo).
2. ``separation='slant_end'`` ou ``'roof_edge'`` -> modelo de "corpo
   equivalente": a regiao de agua morta atras do corpo e substituida por uma
   carenagem de esteira (dois segmentos retos) que fecha em uma ponta T, onde se
   aplica a condicao de Kutta. 'slant_end' supoe escoamento colado na rampa
   (separacao na aresta inferior da rampa, tipico de angulos < ~30 graus) e
   'roof_edge' supoe separacao ja na aresta teto/rampa (tipico de > ~30 graus).
   O valor ~30 graus e o regime critico que recordo da literatura (Ahmed et al.,
   1984); confira no artigo original.
"""
from __future__ import annotations

import numpy as np

from .geometry import _arc, _line

AHMED = dict(L=1.044, H=0.288, R=0.100, slant_length=0.222, clearance=0.050)
"""Dimensoes (m) do corpo de Ahmed como as recordo de Ahmed, Ramm & Faltin (1984)
("Some salient features of the time-averaged ground vehicle wake", SAE 840300).
Confira no artigo original antes de comparacoes quantitativas."""


def _bezier_tangent(P0, t0, P2, end_tangent=None):
    """Gerador de pontos de uma Bezier quadratica de P0 a P2 tangente a t0 em P0
    (ou tangente a end_tangent em P2). O ponto de controle fica no meio em x."""
    P0, P2 = np.asarray(P0, float), np.asarray(P2, float)
    if t0 is not None:
        t = np.asarray(t0, float)
        P1 = P0 + t * (0.5 * (P2[0] - P0[0]) / t[0])
    else:
        t = np.asarray(end_tangent, float)
        P1 = P2 - t * (0.5 * (P2[0] - P0[0]) / t[0])

    def gen(k):
        s = 0.5 * (1.0 - np.cos(np.pi * np.linspace(0.0, 1.0, k + 1)))[:-1, None]
        return (1 - s) ** 2 * P0 + 2 * s * (1 - s) * P1 + s**2 * P2

    return gen


def _curve_length(gen, k=2000):
    pts = gen(k)
    return float(np.sum(np.hypot(np.diff(pts[:, 0]), np.diff(pts[:, 1]))))


def _discretize(segs, n):
    total = sum(length for length, _ in segs)
    counts = [max(4, int(round(n * length / total))) for length, _ in segs]
    return np.vstack([gen(k) for (_, gen), k in zip(segs, counts)]), counts


def ahmed_profile(slant_deg=25.0, n=300, separation="slant_end", wake_length=None,
                  wake_tip_y=0.0, kutta_at="top", dims=AHMED):
    """Perfil lateral 2D do corpo de Ahmed (fundo em y = 0, solo em y = -clearance).

    Parametros
    ----------
    separation : None, 'slant_end' ou 'roof_edge' (ver docstring do modulo).
    wake_length : comprimento da carenagem de esteira alem da base (m). Se None:
        com 'slant_end' e rampa > 0 usa o prolongamento tangente da rampa ate
        y = wake_tip_y; caso contrario usa a altura H.
    wake_tip_y : altura da ponta T da esteira (m).
    kutta_at : apenas para separation=None: 'top' (aresta rampa/base) ou 'bottom'.

    Retorna
    -------
    x, y : nos fechados, sentido horario, comecando no ponto de Kutta.
    info : dict com 'wet' (mascara dos paineis fisicos do corpo), pontos de
        separacao 'P_top', 'P_bot', ponta 'T', e dimensoes.
    """
    L, H, R, ls = dims["L"], dims["H"], dims["R"], dims["slant_length"]
    phi = np.radians(slant_deg)
    if slant_deg > 0:
        xs0, ys1 = L - ls * np.cos(phi), H - ls * np.sin(phi)
    else:  # traseira reta (square-back)
        xs0, ys1 = L, H
    B, S, E = (L, 0.0), (L, ys1), (xs0, H)
    body = [
        (L - R, lambda k: _line(B, (R, 0.0), k)),
        (np.pi * R / 2, lambda k: _arc((R, R), R, 1.5 * np.pi, np.pi, k)),
        (H - 2 * R, lambda k: _line((0.0, R), (0.0, H - R), k)),
        (np.pi * R / 2, lambda k: _arc((R, H - R), R, np.pi, 0.5 * np.pi, k)),
        (xs0 - R, lambda k: _line((R, H), E, k)),
    ]
    info = dict(dims=dims, slant_deg=slant_deg, separation=separation, B=B, S=S, E=E)

    if separation is None:
        if slant_deg > 0:
            body.append((ls, lambda k: _line(E, S, k)))
        body.append((ys1, lambda k: _line(S, B, k)))
        pts, _ = _discretize(body, n)
        if kutta_at == "top":
            i0 = int(np.argmin(np.hypot(pts[:, 0] - S[0], pts[:, 1] - S[1])))
            pts = np.roll(pts, -i0, axis=0)
        elif kutta_at != "bottom":
            raise ValueError("kutta_at deve ser 'top' ou 'bottom'")
        info["wet"] = np.ones(len(pts), bool)
        info.update(P_top=S, P_bot=B, T=None)
    else:
        if separation == "slant_end":
            if slant_deg > 0:
                body.append((ls, lambda k: _line(E, S, k)))
            P = S
        elif separation == "roof_edge":
            P = E
        else:
            raise ValueError("separation deve ser None, 'slant_end' ou 'roof_edge'")
        if wake_length is None:
            if separation == "slant_end" and slant_deg > 0:
                wake_length = (ys1 - wake_tip_y) / np.tan(phi)
            else:
                wake_length = H
        T = (L + wake_length, wake_tip_y)
        # A fronteira da esteira sai TANGENTE a superficie no ponto de separacao
        # (sem quina -> sem singularidade de velocidade): curva de Bezier quadratica.
        t_top = (np.cos(phi), -np.sin(phi)) if P is S and slant_deg > 0 else (1.0, 0.0)
        top_curve = _bezier_tangent(P, t_top, T)
        bot_curve = _bezier_tangent(T, None, B, end_tangent=(1.0, 0.0))
        segs = [(_curve_length(bot_curve), bot_curve)] + body + \
               [(_curve_length(top_curve), top_curve)]
        pts, counts = _discretize(segs, n)
        wet = np.ones(len(pts), bool)
        wet[: counts[0]] = False
        wet[-counts[-1]:] = False
        info.update(wet=wet, P_top=P, P_bot=B, T=T, wake_length=wake_length)
    pts = np.vstack([pts, pts[:1]])
    return pts[:, 0], pts[:, 1], info


def ahmed_forces(sol, info, ref_length=None):
    """Forcas no corpo de Ahmed 2D (por unidade de envergadura, adimensionais com q*ref).

    Retorna dict com:
    * cl_eq, cd_eq : integracao de pressao em todo o contorno (corpo equivalente,
      inclui a carenagem de esteira). cd_eq ~ 0 por d'Alembert.
    * cl_body, cd_body : so nas superficies fisicas molhadas + regiao de agua
      morta (base, e rampa se 'roof_edge') a pressao constante Cp_base, tomada
      como a media do Cp nos dois paineis molhados adjacentes aos pontos de
      separacao (hipotese classica de esteira a pressao constante).
    * cp_base
    """
    p = sol.panels
    ref = ref_length or info["dims"]["L"]
    dfx, dfy = -sol.cp * p.length * p.nx, -sol.cp * p.length * p.ny
    out = dict(cl_eq=dfy.sum() / ref, cd_eq=dfx.sum() / ref)
    wet = info["wet"]
    if info.get("T") is None:
        out.update(cl_body=out["cl_eq"], cd_body=out["cd_eq"], cp_base=np.nan)
        return out
    idx = np.flatnonzero(wet)
    cp_b = 0.5 * (sol.cp[idx[0]] + sol.cp[idx[-1]])
    (xP, yP), (xB, yB) = info["P_top"], info["P_bot"]
    # integral de n dl no caminho horario P -> B (regiao de agua morta)
    nlx, nly = -(yB - yP), (xB - xP)
    fx = dfx[wet].sum() - cp_b * nlx
    fy = dfy[wet].sum() - cp_b * nly
    out.update(cl_body=fy / ref, cd_body=fx / ref, cp_base=cp_b)
    return out
