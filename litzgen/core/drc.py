"""Geometric design-rule check for generated coils.

All strands of a coil normally share one net, so KiCad's own DRC cannot see
a strand-to-strand short (which would silently destroy the Litz action).
This checker treats every strand as its own conductor and verifies:

* copper clearance between different strands on every layer
  (tracks, arcs and via pads, modelled as capsules / discs);
* hole-to-hole clearance between all vias.

It is pure numpy/scipy and runs in well under a second for the paper coil.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from .neighbors import query_pairs

from .geometry import CoilGeometry


@dataclass
class Violation:
    kind: str
    layer: str
    pos: Tuple[float, float]
    strands: Tuple[int, int]
    actual: float
    required: float

    def __str__(self) -> str:
        return "%s on %s at (%.3f, %.3f) mm: strands %d/%d, %.4f < %.4f mm" % (
            self.kind, self.layer, self.pos[0], self.pos[1], self.strands[0], self.strands[1],
            self.actual, self.required)


def _arc_points(start, mid, end, max_step_deg=0.25):
    (x1, y1), (x2, y2), (x3, y3) = start, mid, end
    d = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2))
    if abs(d) < 1e-12:
        return [start, end]
    ux = ((x1 ** 2 + y1 ** 2) * (y2 - y3) + (x2 ** 2 + y2 ** 2) * (y3 - y1) + (x3 ** 2 + y3 ** 2) * (y1 - y2)) / d
    uy = ((x1 ** 2 + y1 ** 2) * (x3 - x2) + (x2 ** 2 + y2 ** 2) * (x1 - x3) + (x3 ** 2 + y3 ** 2) * (x2 - x1)) / d
    r = math.hypot(x1 - ux, y1 - uy)
    a1 = math.atan2(y1 - uy, x1 - ux)
    a2 = math.atan2(y2 - uy, x2 - ux)
    a3 = math.atan2(y3 - uy, x3 - ux)

    def ccw(a, b):
        return (b - a) % (2 * math.pi)
    if ccw(a1, a2) <= ccw(a1, a3):
        sweep = ccw(a1, a3)
    else:
        sweep = -ccw(a3, a1)
    n = max(1, int(math.ceil(abs(math.degrees(sweep)) / max_step_deg)))
    pts = [(ux + r * math.cos(a1 + sweep * k / n), uy + r * math.sin(a1 + sweep * k / n)) for k in range(n + 1)]
    pts[0] = start
    pts[-1] = end
    return pts


def _pt_seg(P, A, B):
    d = B - A
    L2 = np.einsum("ij,ij->i", d, d)
    t = np.where(L2 > 1e-18, np.einsum("ij,ij->i", P - A, d) / np.where(L2 > 1e-18, L2, 1.0), 0.0)
    t = np.clip(t, 0.0, 1.0)
    C = A + d * t[:, None]
    return np.linalg.norm(P - C, axis=1), C


def seg_seg_distance(P0, P1, Q0, Q1):
    """Vectorised minimum distance between 2D segments P0-P1 and Q0-Q1 (arrays N x 2).

    Returns (distance, midpoint of the closest pair). Zero-length segments are points.
    """
    cands = []
    d, c = _pt_seg(P0, Q0, Q1); cands.append((d, 0.5 * (P0 + c)))
    d, c = _pt_seg(P1, Q0, Q1); cands.append((d, 0.5 * (P1 + c)))
    d, c = _pt_seg(Q0, P0, P1); cands.append((d, 0.5 * (Q0 + c)))
    d, c = _pt_seg(Q1, P0, P1); cands.append((d, 0.5 * (Q1 + c)))
    D = np.stack([x[0] for x in cands])
    k = np.argmin(D, axis=0)
    idx = np.arange(D.shape[1])
    dist = D[k, idx]
    at = np.stack([x[1] for x in cands])[k, idx]

    def orient(a, b, c):
        return (b[:, 0] - a[:, 0]) * (c[:, 1] - a[:, 1]) - (b[:, 1] - a[:, 1]) * (c[:, 0] - a[:, 0])
    o1 = orient(P0, P1, Q0); o2 = orient(P0, P1, Q1)
    o3 = orient(Q0, Q1, P0); o4 = orient(Q0, Q1, P1)
    cross = (o1 * o2 < 0) & (o3 * o4 < 0)
    dist = np.where(cross, 0.0, dist)
    return dist, at


def _layer_primitives(geo: CoilGeometry):
    p = geo.params
    R = p.stackup.n_layers
    prims = {l: {"a": [], "b": [], "r": [], "s": []} for l in range(R)}
    for t in geo.tracks:
        d = prims[t.layer]
        d["a"].append(t.start); d["b"].append(t.end); d["r"].append(t.width / 2); d["s"].append(t.strand)
    for a in geo.arcs:
        pts = _arc_points(a.start, a.mid, a.end)
        d = prims[a.layer]
        for i in range(len(pts) - 1):
            d["a"].append(pts[i]); d["b"].append(pts[i + 1]); d["r"].append(a.width / 2); d["s"].append(a.strand)
    for v in geo.vias:
        for l in range(min(v.la, v.lb), max(v.la, v.lb) + 1):
            d = prims[l]
            d["a"].append(v.pos); d["b"].append(v.pos); d["r"].append(v.pad / 2); d["s"].append(v.strand)
    out = {}
    for l, d in prims.items():
        if d["a"]:
            out[l] = (np.array(d["a"], float), np.array(d["b"], float), np.array(d["r"], float), np.array(d["s"], int))
    return out


def check(geo: CoilGeometry, clearance: float = None, max_report: int = 50) -> List[Violation]:
    p = geo.params
    clr = p.clearance if clearance is None else clearance
    tol = 1e-4
    names = p.stackup.kicad_layers
    viol: List[Violation] = []
    for l, (A, B, Rr, S) in _layer_primitives(geo).items():
        M = 0.5 * (A + B)
        half = 0.5 * np.linalg.norm(B - A, axis=1)
        reach = 2 * half.max() + 2 * Rr.max() + clr + tol
        pairs = query_pairs(M, reach)
        if len(pairs) == 0:
            continue
        i, j = pairs[:, 0], pairs[:, 1]
        keep = S[i] != S[j]
        i, j = i[keep], j[keep]
        # cheap prefilter on midpoint distance
        dm = np.linalg.norm(M[i] - M[j], axis=1)
        keep = dm <= half[i] + half[j] + Rr[i] + Rr[j] + clr + tol
        i, j = i[keep], j[keep]
        if len(i) == 0:
            continue
        dist, at = seg_seg_distance(A[i], B[i], A[j], B[j])
        gap = dist - Rr[i] - Rr[j]
        bad = np.nonzero(gap < clr - tol)[0]
        for k in bad[:max_report]:
            viol.append(Violation("clearance", names[l], tuple(at[k]), (int(S[i[k]]), int(S[j[k]])),
                                  float(gap[k]), clr))
    # hole to hole (edge-to-edge)
    if geo.vias:
        P = np.array([v.pos for v in geo.vias])
        drill = np.array([v.drill for v in geo.vias])
        need = drill.max() + p.hole_clearance
        pairs = query_pairs(P, need + tol)
        for i, j in pairs[:max_report]:
            d = float(np.linalg.norm(P[i] - P[j])) - 0.5 * (drill[i] + drill[j])
            if d < p.hole_clearance - tol:
                viol.append(Violation("hole-to-hole", "drill", tuple(P[i]),
                                      (geo.vias[i].strand, geo.vias[j].strand), d, p.hole_clearance))
    return viol
