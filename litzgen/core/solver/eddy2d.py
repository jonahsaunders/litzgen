"""Axisymmetric (ring-filament) eddy-current solver for coil cross-sections.

A cross-section of the coil at one azimuth contains, for every turn, the
rectangular strand conductors (and lane positions in transposition zones).
Treating each conductor locally as a full ring ("local axisymmetric"
approximation) every conductor is subdivided into n_w x n_t ring filaments
with cosine grading toward the edges (skin depth at 6.78 MHz is 25 um).

Per radian of azimuth the filament impedance matrix is

    Z_f = diag(rho r_i / A_i) + j w M_ij / (2 pi),

    M_ij = mu0 sqrt(r_i r_j) [ (2/k - k) K(k) - (2/k) E(k) ],
    k^2  = 4 r_i r_j / ((r_i + r_j)^2 + (z_i - z_j)^2),
    M_ii = mu0 r_i [ ln(8 r_i / GMD_i) - 2 ],  GMD = 0.2235 (dw + dt).

Filaments of one conductor share one voltage, so with the incidence matrix B
the conductor impedance matrix is Z_c = (B^T Z_f^-1 B)^-1.  Re(Z_c) contains
DC resistance, skin effect, proximity effect and the mutual screening of
all strands of all turns; Im(Z_c)/w is the conductor inductance matrix.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional, Tuple

import numpy as np

from ..params import MU0


def ellipke_agm(k2: np.ndarray, kp2: np.ndarray, iters: int = 12) -> Tuple[np.ndarray, np.ndarray]:
    """Complete elliptic integrals K(k), E(k) (parameter m = k^2) by the AGM.

    kp2 = 1 - k^2 is passed separately to avoid cancellation for close filaments.
    """
    a = np.ones_like(k2)
    b = np.sqrt(kp2)
    c2sum = 0.5 * k2              # 2^{-1} c0^2
    pw = 0.5
    for _ in range(iters):
        an = 0.5 * (a + b)
        bn = np.sqrt(a * b)
        c = 0.5 * (a - b)
        pw *= 2.0
        c2sum = c2sum + pw * c * c
        a, b = an, bn
    K = math.pi / (2.0 * a)
    E = K * (1.0 - c2sum)
    return K, E


def ring_mutual(r1, z1, r2, z2) -> np.ndarray:
    """Mutual inductance of coaxial circular filaments (H), vectorised."""
    num = (r1 + r2) ** 2 + (z1 - z2) ** 2
    k2 = 4.0 * r1 * r2 / num
    kp2 = ((r1 - r2) ** 2 + (z1 - z2) ** 2) / num
    K, E = ellipke_agm(k2, kp2)
    k = np.sqrt(k2)
    return MU0 * np.sqrt(r1 * r2) * ((2.0 / k - k) * K - (2.0 / k) * E)


def graded_edges(n: int) -> np.ndarray:
    """n cells on [-0.5, 0.5], cosine-graded (fine at both edges)."""
    return -0.5 * np.cos(np.pi * np.arange(n + 1) / n)


@dataclass
class CrossSection:
    """Conductors: arrays of centre radius r, centre z, width w (radial), thickness t (all metres)."""
    r: np.ndarray
    z: np.ndarray
    w: np.ndarray
    t: np.ndarray

    @property
    def n(self) -> int:
        return len(self.r)


@dataclass
class EddyResult:
    Zc: np.ndarray               # conductor impedance per radian (ohm/rad), N_c x N_c
    fil_r: np.ndarray
    fil_z: np.ndarray
    fil_dw: np.ndarray
    fil_dt: np.ndarray
    fil_cond: np.ndarray
    X: Optional[np.ndarray]      # Z_f^-1 B (filament current per unit conductor voltage)

    def filament_currents(self, I_c: np.ndarray) -> np.ndarray:
        V = self.Zc @ I_c
        return self.X @ V

    def current_density(self, I_c: np.ndarray) -> np.ndarray:
        return self.filament_currents(I_c) / (self.fil_dw * self.fil_dt)


def mesh_filaments(cs: CrossSection, nw: int, nt: int):
    ew = graded_edges(nw)
    et = graded_edges(nt)
    cw = 0.5 * (ew[1:] + ew[:-1]); dw = np.diff(ew)
    ct = 0.5 * (et[1:] + et[:-1]); dt = np.diff(et)
    CW, CT = np.meshgrid(cw, ct, indexing="ij")
    DW, DT = np.meshgrid(dw, dt, indexing="ij")
    CW, CT, DW, DT = CW.ravel(), CT.ravel(), DW.ravel(), DT.ravel()
    m = len(CW)
    fr = (cs.r[:, None] + cs.w[:, None] * CW[None, :]).ravel()
    fz = (cs.z[:, None] + cs.t[:, None] * CT[None, :]).ravel()
    fdw = (cs.w[:, None] * DW[None, :]).ravel()
    fdt = (cs.t[:, None] * DT[None, :]).ravel()
    fc = np.repeat(np.arange(cs.n), m)
    return fr, fz, fdw, fdt, fc


def _F4(x, y):
    """Fourth antiderivative: d^4 F / dx^2 dy^2 = ln(x^2 + y^2)."""
    x2, y2 = x * x, y * y
    r2 = x2 + y2
    with np.errstate(divide="ignore", invalid="ignore"):
        L = np.where(r2 > 0, np.log(np.where(r2 > 0, r2, 1.0)), 0.0)
        t1 = np.where(x != 0, np.arctan(y / np.where(x != 0, x, 1.0)), 0.0)
        t2 = np.where(y != 0, np.arctan(x / np.where(y != 0, y, 1.0)), 0.0)
    return (0.25 * x2 * y2 * L - (x2 * x2 + y2 * y2) * L / 24.0
            + (x2 * x * y * t1 + x * y2 * y * t2) / 3.0 - 25.0 / 24.0 * x2 * y2)


def log_gmd_rect(xa, wa, ya, ta, xb, wb, yb, tb):
    """Exact mean of ln(distance) between two axis-aligned rectangles (centre, width, height).

    Vectorised; coordinates are rescaled per pair to limit cancellation.
    """
    sc = np.maximum.reduce([wa, ta, wb, tb])
    xa, wa, ya, ta = xa / sc, wa / sc, ya / sc, ta / sc
    xb, wb, yb, tb = xb / sc, wb / sc, yb / sc, tb / sc
    a, b = xa - wa / 2, xa + wa / 2
    c, d = xb - wb / 2, xb + wb / 2
    e, f = ya - ta / 2, ya + ta / 2
    g, h = yb - tb / 2, yb + tb / 2
    U = [(b - c, 1.0), (a - c, -1.0), (b - d, -1.0), (a - d, 1.0)]
    V = [(f - g, 1.0), (e - g, -1.0), (f - h, -1.0), (e - h, 1.0)]
    tot = 0.0
    for u, su in U:
        for v, sv in V:
            tot = tot + su * sv * _F4(u, v)
    mean_ln_r2 = tot / (wa * ta * wb * tb)
    return 0.5 * mean_ln_r2 + np.log(sc)


def filament_inductance(fr, fz, fdw, fdt, near_cells: float = 4.0) -> np.ndarray:
    """Ring-filament inductance matrix (H per full ring).

    Far pairs use the filament formula at the centres; pairs closer than
    near_cells cell sizes get the exact rectangle-to-rectangle geometric mean
    distance in the local logarithmic term, which keeps elongated cells accurate.
    """
    n = len(fr)
    M = np.empty((n, n))
    blk = 1024
    for a in range(0, n, blk):
        b = min(n, a + blk)
        M[a:b] = _mut_block(fr[a:b], fz[a:b], fr, fz)
    # self terms with exact self-GMD
    lg = log_gmd_rect(fr, fdw, fz, fdt, fr, fdw, fz, fdt)
    M[np.arange(n), np.arange(n)] = MU0 * fr * (np.log(8.0 * fr) - lg - 2.0)
    # near-pair correction
    from ..neighbors import query_pairs
    size = np.maximum(fdw, fdt)
    pairs = query_pairs(np.stack([fr, fz], axis=1), near_cells * size.max())
    if len(pairs):
        i, j = pairs[:, 0], pairs[:, 1]
        d = np.hypot(fr[i] - fr[j], fz[i] - fz[j])
        sel = d < near_cells * np.maximum(size[i], size[j])
        i, j, d = i[sel], j[sel], d[sel]
        lg = log_gmd_rect(fr[i], fdw[i], fz[i], fdt[i], fr[j], fdw[j], fz[j], fdt[j])
        corr = MU0 * np.sqrt(fr[i] * fr[j]) * (np.log(d) - lg)
        M[i, j] += corr
        M[j, i] += corr
    return M


def _mut_block(r1, z1, r2, z2):
    R1 = r1[:, None]; Z1 = z1[:, None]
    R2 = r2[None, :]; Z2 = z2[None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        return ring_mutual(R1, Z1, R2, Z2)


def solve_cross_section(cs: CrossSection, freq: float, rho: float, nw: int = 10, nt: int = 3,
                        keep_currents: bool = False, M: Optional[np.ndarray] = None) -> EddyResult:
    fr, fz, fdw, fdt, fc = mesh_filaments(cs, nw, nt)
    if M is None:
        M = filament_inductance(fr, fz, fdw, fdt)
    w = 2.0 * math.pi * freq
    Zf = (1j * w / (2.0 * math.pi)) * M
    Zf[np.diag_indices_from(Zf)] += rho * fr / (fdw * fdt)
    B = np.zeros((len(fr), cs.n))
    B[np.arange(len(fr)), fc] = 1.0
    X = np.linalg.solve(Zf, B)
    Y = B.T @ X
    Zc = np.linalg.inv(Y)
    Zc = 0.5 * (Zc + Zc.T)
    return EddyResult(Zc, fr, fz, fdw, fdt, fc, X if keep_currents else None)
