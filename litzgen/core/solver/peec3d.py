"""Strand-level 3D partial-inductance (PEEC) model.

Each strand is discretised into straight filament segments (horizontal
trace pieces and vertical via barrels).  The strand inductance matrix is

    L[s, t] = sum_{i in s} sum_{j in t} Lp[i, j],

    Lp[i, j] = mu0/(4 pi) * integral integral  dl_i . dl_j / R      (Neumann)

evaluated with a three-tier scheme:

* far pairs      midpoint rule (relative error < ~1e-3 beyond 5 lengths)
* mid pairs      4 x 4 Gauss-Legendre along both segments
* near pairs     8 x 8 Gauss-Legendre, plus 2 x 2 points across the trace
                 width for different strands (a width-averaged, GMD-like
                 kernel) or the reduced kernel sqrt(R^2 + a^2) for pieces of
                 the same strand
* self terms     closed form of the reduced kernel:
                 mu0/(2 pi) [ l asinh(l/a) - sqrt(l^2 + a^2) + a ],
                 a = 0.2235 (w + t) (rectangle GMD) or the barrel radius.
                 For l >> a this reproduces Grover's bar formula
                 mu0 l/(2 pi) [ln(2l/(w+t)) + 0.5].

Only the S x S strand matrix is accumulated, so memory stays O(S^2 + N).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Tuple

import numpy as np
from ..neighbors import query_pairs

from ..geometry import CoilGeometry, Dwell, Jog, ViaHop, lateral_offset, phi_of
from ..params import MU0

K = MU0 / (4.0 * math.pi)


@dataclass
class SegmentMesh:
    p0: np.ndarray        # N x 3 (m)
    p1: np.ndarray        # N x 3 (m)
    strand: np.ndarray    # N
    width: np.ndarray     # N (m) lateral width (0 for vias)
    a_eq: np.ndarray      # N (m) equivalent radius for the reduced kernel
    is_via: np.ndarray    # N bool

    @property
    def n(self) -> int:
        return len(self.strand)

    @property
    def d(self) -> np.ndarray:
        return self.p1 - self.p0

    @property
    def length(self) -> np.ndarray:
        return np.linalg.norm(self.p1 - self.p0, axis=1)

    @property
    def mid(self) -> np.ndarray:
        return 0.5 * (self.p0 + self.p1)

    def lateral(self) -> np.ndarray:
        """Unit vector across the trace width (in-plane, normal to the segment)."""
        d = self.d
        lat = np.stack([-d[:, 1], d[:, 0], np.zeros(len(d))], axis=1)
        n = np.linalg.norm(lat, axis=1)
        lat[n > 0] /= n[n > 0, None]
        return lat


def build_mesh(geo: CoilGeometry, max_seg_len: float = 3.0, min_seg_len: float = 0.25) -> SegmentMesh:
    """Discretise all strand paths into straight 3D segments (inputs in mm, output in m)."""
    p = geo.params
    z = p.stackup.layer_z()
    t_cu = p.stackup.copper_thickness
    P0, P1, S, W, A, V = [], [], [], [], [], []

    def xyz(theta, u, layer):
        r = p.radius(theta) + lateral_offset(p, u)
        ph = phi_of(p, theta)
        return (r * math.cos(ph), r * math.sin(ph), z[layer])

    for s, pcs in enumerate(geo.pieces):
        pts: List[Tuple[float, float, float]] = []   # running horizontal polyline
        layer_now = None

        def flush():
            # merge too-short pieces, then emit horizontal segments
            if len(pts) < 2:
                return
            keep = [pts[0]]
            for q in pts[1:-1]:
                if math.dist(q, keep[-1]) >= min_seg_len:
                    keep.append(q)
            if math.dist(pts[-1], keep[-1]) < min_seg_len and len(keep) > 1:
                keep[-1] = pts[-1]
            else:
                keep.append(pts[-1])
            wdt = p.strand_width
            a = 0.2235 * (wdt + t_cu)
            for q0, q1 in zip(keep[:-1], keep[1:]):
                if math.dist(q0, q1) < 1e-9:
                    continue
                P0.append(q0); P1.append(q1); S.append(s); W.append(wdt); A.append(a); V.append(False)

        for pc in pcs:
            if isinstance(pc, (Dwell, Jog)):
                ta, tb = pc.ta, pc.tb
                ua, ub = (pc.u, pc.u) if isinstance(pc, Dwell) else (pc.ua, pc.ub)
                rr = p.radius(0.5 * (ta + tb)) + lateral_offset(p, 0.5 * (ua + ub))
                n = max(1, int(math.ceil(abs(tb - ta) * rr / max_seg_len)))
                for k in range(n + 1):
                    f = k / n
                    q = xyz(ta + f * (tb - ta), ua + f * (ub - ua), pc.layer)
                    if pts and math.dist(pts[-1], q) < 1e-9:
                        continue
                    pts.append(q)
                layer_now = pc.layer
            elif isinstance(pc, ViaHop):
                flush()
                q0 = xyz(pc.t, pc.u, pc.la)
                q1 = xyz(pc.t, pc.u, pc.lb)
                P0.append(q0); P1.append(q1); S.append(s)
                W.append(0.0)
                A.append(max(p.via.drill / 2.0, 1e-3)); V.append(True)
                pts = [q1]
        flush()
    mm = 1e-3
    return SegmentMesh(np.array(P0) * mm, np.array(P1) * mm, np.array(S, int),
                       np.array(W) * mm, np.array(A) * mm, np.array(V, bool))


# ----------------------------------------------------------------------------- kernels
def self_partial(l: np.ndarray, a: np.ndarray) -> np.ndarray:
    """Self partial inductance with the reduced kernel (H)."""
    return 2.0 * K * (l * np.arcsinh(l / a) - np.sqrt(l * l + a * a) + a)


def _gl(n: int):
    x, w = np.polynomial.legendre.leggauss(n)
    return 0.5 * (x + 1.0), 0.5 * w


def pair_quadrature(mesh: SegmentMesh, i: np.ndarray, j: np.ndarray, n: int, lateral: bool,
                    reduced: bool, chunk: int = 4000) -> np.ndarray:
    """Neumann integral for pairs (i, j) with n-point GL along each segment.

    lateral: additionally average over 2 points across each trace width.
    reduced: use the reduced kernel sqrt(R^2 + a_i a_j) (same-conductor pairs).
    """
    xs, ws = _gl(n)
    d = mesh.d
    lat = mesh.lateral()
    if lateral:
        ly = np.array([-0.5, 0.5]) / math.sqrt(3.0)
        lw = np.array([0.5, 0.5])
    else:
        ly = np.array([0.0])
        lw = np.array([1.0])
    out = np.empty(len(i))
    for c0 in range(0, len(i), chunk):
        ii = i[c0:c0 + chunk]
        jj = j[c0:c0 + chunk]
        dot = np.einsum("ij,ij->i", d[ii], d[jj])
        # points: pair x (n * nl)
        Pi = (mesh.p0[ii][:, None, None, :] + d[ii][:, None, None, :] * xs[None, :, None, None]
              + (mesh.width[ii][:, None, None, None] * ly[None, None, :, None]) * lat[ii][:, None, None, :])
        Pj = (mesh.p0[jj][:, None, None, :] + d[jj][:, None, None, :] * xs[None, :, None, None]
              + (mesh.width[jj][:, None, None, None] * ly[None, None, :, None]) * lat[jj][:, None, None, :])
        Pi = Pi.reshape(len(ii), -1, 3)
        Pj = Pj.reshape(len(jj), -1, 3)
        wi = (ws[:, None] * lw[None, :]).reshape(-1)
        diff = Pi[:, :, None, :] - Pj[:, None, :, :]
        R2 = np.einsum("pabk,pabk->pab", diff, diff)
        if reduced:
            R2 = R2 + (mesh.a_eq[ii] * mesh.a_eq[jj])[:, None, None]
        R2 = np.maximum(R2, 1e-24)
        integ = np.einsum("a,b,pab->p", wi, wi, 1.0 / np.sqrt(R2))
        out[c0:c0 + chunk] = K * dot * integ
    return out


def strand_inductance(mesh: SegmentMesh, n_strands: int, near_factor: float = 1.5,
                      mid_factor: float = 5.0, block: int = 400) -> np.ndarray:
    """S x S strand inductance matrix (H)."""
    N = mesh.n
    mid = mesh.mid
    d = mesh.d
    l = mesh.length
    Onehot = np.zeros((N, n_strands))
    Onehot[np.arange(N), mesh.strand] = 1.0
    L = np.zeros((n_strands, n_strands))
    sq = np.einsum("ij,ij->i", mid, mid)
    # 1) midpoint rule for all pairs (i != j)
    for a0 in range(0, N, block):
        a1 = min(N, a0 + block)
        D2 = sq[a0:a1, None] + sq[None, :] - 2.0 * (mid[a0:a1] @ mid.T)
        D = np.sqrt(np.maximum(D2, 0.0))
        G = d[a0:a1] @ d.T
        with np.errstate(divide="ignore", invalid="ignore"):
            G = np.where(D > 0, G / D, 0.0)
        idx = np.arange(a0, a1)
        G[idx - a0, idx] = 0.0
        L += Onehot[a0:a1].T @ (K * G) @ Onehot
    # 2) corrections for mid / near pairs
    lmax = max(l.max(), 3.0 * mesh.width.max())
    pairs = query_pairs(mid, mid_factor * lmax)
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        Dm = np.linalg.norm(mid[i] - mid[j], axis=1)
        # pairs within a few trace widths need width averaging even when the segments are short
        Lp = np.maximum(np.maximum(l[i], l[j]), 3.0 * np.maximum(mesh.width[i], mesh.width[j]))
        sel = Dm < mid_factor * Lp
        i, j, Dm, Lp = i[sel], j[sel], Dm[sel], Lp[sel]
        mp = K * np.einsum("ij,ij->i", d[i], d[j]) / Dm
        exact = np.empty(len(i))
        near = Dm < near_factor * Lp
        same = mesh.strand[i] == mesh.strand[j]
        trace_pair = (~mesh.is_via[i]) & (~mesh.is_via[j])
        m = ~near
        if m.any():
            exact[m] = pair_quadrature(mesh, i[m], j[m], 4, lateral=False, reduced=False)
        m = near & same
        if m.any():
            exact[m] = pair_quadrature(mesh, i[m], j[m], 8, lateral=False, reduced=True)
        m = near & ~same & trace_pair
        if m.any():
            exact[m] = pair_quadrature(mesh, i[m], j[m], 8, lateral=True, reduced=False, chunk=1500)
        m = near & ~same & ~trace_pair
        if m.any():
            exact[m] = pair_quadrature(mesh, i[m], j[m], 8, lateral=False, reduced=False)
        corr = exact - mp
        si, sj = mesh.strand[i], mesh.strand[j]
        np.add.at(L, (si, sj), corr)
        np.add.at(L, (sj, si), corr)
    # 3) self terms
    ls = self_partial(l, mesh.a_eq)
    np.add.at(L, (mesh.strand, mesh.strand), ls)
    return 0.5 * (L + L.T)


def via_resistance(geo: CoilGeometry, rho: float, freq: float) -> Tuple[np.ndarray, float]:
    """AC resistance added to each strand by its vias, and the resistance of one via (ohm)."""
    p = geo.params
    z = p.stackup.layer_z()
    delta = math.sqrt(rho / (math.pi * freq * MU0))
    d_o = p.via.drill * 1e-3
    tp = min(p.via.plating * 1e-3, d_o / 2 * 0.999)
    area = math.pi / 4 * (d_o ** 2 - (d_o - 2 * tp) ** 2)
    Dl = tp / delta
    F = Dl * (math.sinh(2 * Dl) + math.sin(2 * Dl)) / (math.cosh(2 * Dl) - math.cos(2 * Dl))
    per = np.zeros(geo.n_strands)
    r_one = []
    for v in geo.vias:
        h = abs(z[v.la] - z[v.lb]) * 1e-3
        r = rho * h / area * F
        per[v.strand] += r
        r_one.append(r)
    return per, (float(np.mean(r_one)) if r_one else 0.0)
