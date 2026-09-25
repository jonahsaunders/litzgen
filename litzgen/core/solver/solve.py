"""Coil simulation: combines the 3D PEEC inductance model with the
axisymmetric eddy-current resistance model and solves for strand current
sharing.

Strand-level circuit (all strands in parallel between the two terminals):

    (R_eddy + R_via + j w L_3D) I = V 1,        Z_eq = 1 / (1^T Z^-1 1)

R_eddy (S x S) is assembled by sweeping the azimuth: at each sampled angle
the position of every strand in every turn is read from its path, the
cross-section is solved (cached per occupancy pattern and interpolated
between a few representative angles), and Re(Z_c) is scattered into the
strand matrix.  Off-diagonal terms of R_eddy are the proximity coupling
between strands; they are what makes transposition matter.
"""
from __future__ import annotations

import math
import time
from dataclasses import dataclass, field
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np

from ..geometry import CoilGeometry, lateral_offset, strand_state_at
from ..params import CoilParams, SimParams
from .eddy2d import CrossSection, EddyResult, filament_inductance, mesh_filaments, solve_cross_section
from .peec3d import build_mesh, strand_inductance, via_resistance

Progress = Optional[Callable[[str, float], None]]


@dataclass
class SimResult:
    frequency: float
    L: float                         # H
    R: float                         # ohm
    Q: float
    R_dc: float
    strand_currents: np.ndarray      # complex, normalised to total current 1 A
    ring_of_strand: List[int]
    L_3d_strand: np.ndarray
    R_strand: np.ndarray
    R_via_strand: np.ndarray
    R_via_one: float
    R_ideal_sharing: float           # R if every strand carried exactly 1/S
    L_axisym: float                  # coil inductance from the ring model (cross-check)
    loss_split: Dict[str, float]     # fraction of loss: dc, eddy (skin+proximity), via, circulating
    timings: Dict[str, float] = field(default_factory=dict)
    patterns: int = 0
    sweep: List[Tuple[float, float, float, float]] = field(default_factory=list)  # (f, L, R, Q)
    cross_section: Optional[dict] = None

    def ring_share(self) -> Dict[int, complex]:
        out: Dict[int, complex] = {}
        for i, r in zip(self.strand_currents, self.ring_of_strand):
            out[r] = out.get(r, 0) + i
        return out

    def to_dict(self) -> dict:
        I = self.strand_currents
        S = len(I)
        return {
            "frequency_hz": self.frequency,
            "L_uH": self.L * 1e6,
            "ESR_ohm": self.R,
            "Q": self.Q,
            "R_dc_ohm": self.R_dc,
            "Rac_over_Rdc": self.R / self.R_dc,
            "R_ideal_sharing_ohm": self.R_ideal_sharing,
            "L_axisymmetric_uH": self.L_axisym * 1e6,
            "via_R_mohm_each": self.R_via_one * 1e3,
            "loss_split": self.loss_split,
            "strand_current_rel_mag": [float(abs(x) * S) for x in I],
            "strand_current_phase_deg": [float(np.degrees(np.angle(x))) for x in I],
            "ring_current_share": {str(k): float(abs(v)) for k, v in sorted(self.ring_share().items())},
            "patterns_solved": self.patterns,
            "timings_s": self.timings,
            "sweep": [{"f_hz": f, "L_uH": l * 1e6, "ESR_ohm": r, "Q": q} for f, l, r, q in self.sweep],
        }


class EddyAssembler:
    """Builds the strand resistance matrix from cached cross-section solves."""

    def __init__(self, geo: CoilGeometry, sim: SimParams):
        self.geo = geo
        self.p = geo.params
        self.sim = sim
        self.n_rep = max(1, sim.eddy_angles)
        self.rep = [2 * math.pi * j / self.n_rep for j in range(self.n_rep)]
        self.cache: Dict[tuple, List[Optional[EddyResult]]] = {}
        self.M_cache: Dict[tuple, np.ndarray] = {}
        self.n_turn_slots = int(math.ceil(self.p.n_turns - 1e-9))
        z = self.p.stackup.layer_z()
        self.z = z

    def _cross_section(self, key: tuple, theta: float) -> CrossSection:
        p = self.p
        r, zz = [], []
        for m, turn in enumerate(key):
            if turn is None:
                continue
            rc = p.radius(theta + 2 * math.pi * m)
            for (u, l) in turn:
                r.append(rc + lateral_offset(p, u))
                zz.append(self.z[l])
        n = len(r)
        mm = 1e-3
        return CrossSection(np.array(r) * mm, np.array(zz) * mm,
                            np.full(n, p.strand_width * mm), np.full(n, p.stackup.copper_thickness * mm))

    def _union(self, j: int):
        """Filament inductance matrix for every possible conductor position at rep angle j."""
        if j in self.M_cache:
            return self.M_cache[j]
        p = self.p
        C, R = p.n_cols, p.stackup.n_layers
        allpos = tuple(tuple((u, l) for l in range(R) for u in range(-1, C + 1)) for _ in range(self.n_turn_slots))
        n_fil = sum(len(t) for t in allpos) * self.sim.eddy_nw * self.sim.eddy_nt
        if n_fil * n_fil * 8 * self.n_rep > self.sim.max_cache_mb * 1e6:
            raise KeyError("union matrix too large to cache")
        cs = self._cross_section(allpos, self.rep[j])
        fr, fz, fdw, fdt, fc = mesh_filaments(cs, self.sim.eddy_nw, self.sim.eddy_nt)
        M = filament_inductance(fr, fz, fdw, fdt)
        per = self.sim.eddy_nw * self.sim.eddy_nt
        index = {}
        c = 0
        for m, turn in enumerate(allpos):
            for pos in turn:
                index[(m, pos)] = c
                c += 1
        self.M_cache[j] = (M, index, per)
        return self.M_cache[j]

    def solve(self, key: tuple, j: int, freq: float) -> EddyResult:
        ent = self.cache.setdefault((key, freq), [None] * self.n_rep)
        if ent[j] is None:
            cs = self._cross_section(key, self.rep[j])
            M = None
            try:
                Mu, index, per = self._union(j)
                cidx = [index[(m, pos)] for m, turn in enumerate(key) if turn is not None for pos in turn]
                fi = (np.array(cidx)[:, None] * per + np.arange(per)[None, :]).ravel()
                M = Mu[np.ix_(fi, fi)]
            except KeyError:
                M = None      # position outside the union (custom schemes): compute directly
            ent[j] = solve_cross_section(cs, freq, self.sim.rho, self.sim.eddy_nw, self.sim.eddy_nt, M=M)
        return ent[j]

    def assemble(self, freq: float, progress: Progress = None):
        p, geo = self.p, self.geo
        S = geo.n_strands
        d = math.radians(self.sim.slice_deg)
        n_s = int(round(2 * math.pi / d))
        d = 2 * math.pi / n_s
        R = np.zeros((S, S))
        Lax = np.zeros((S, S))
        w = 2 * math.pi * freq
        for k in range(n_s):
            th = (k + 0.5) * d
            key_turns = []
            owners = []       # per conductor: (strand, length factor)
            for m in range(self.n_turn_slots):
                tm = th + 2 * math.pi * m
                if tm < 0 or tm > p.theta_end:
                    key_turns.append(None)
                    continue
                pos = []
                for s in range(S):
                    st = strand_state_at(p, geo.pieces[s], tm)
                    if st is None:
                        continue
                    u, l, f = st
                    pos.append(((int(round(u)), l), s, f))
                pos.sort()
                key_turns.append(tuple(x[0] for x in pos))
                owners += [(s, f) for (_, s, f) in pos]
            key = tuple(key_turns)
            # periodic linear interpolation between representative angles
            x = th / (2 * math.pi) * self.n_rep
            j0 = int(math.floor(x)) % self.n_rep
            j1 = (j0 + 1) % self.n_rep
            a = x - math.floor(x)
            Z0 = self.solve(key, j0, freq).Zc
            Z1 = Z0 if self.n_rep == 1 else self.solve(key, j1, freq).Zc
            Zc = (1 - a) * Z0 + a * Z1
            idx = np.array([o[0] for o in owners])
            fac = np.sqrt(np.array([o[1] for o in owners]))
            Rc = Zc.real * fac[:, None] * fac[None, :]
            Lc = Zc.imag / w * fac[:, None] * fac[None, :]
            P = np.zeros((len(idx), S))
            P[np.arange(len(idx)), idx] = 1.0
            R += d * (P.T @ Rc @ P)
            Lax += d * (P.T @ Lc @ P)
            if progress and k % 40 == 0:
                progress("eddy-current slices", k / n_s)
        return 0.5 * (R + R.T), 0.5 * (Lax + Lax.T)

    @property
    def n_patterns(self) -> int:
        return len({k for (k, f) in self.cache})


def _dc_resistance(geo: CoilGeometry, rho: float) -> np.ndarray:
    p = geo.params
    A = p.strand_width * p.stackup.copper_thickness * 1e-6
    Ls = geo.strand_lengths() * 1e-3
    return rho * Ls / A


def _solve_network(Z: np.ndarray):
    S = Z.shape[0]
    one = np.ones(S)
    y = np.linalg.solve(Z, one)
    Zeq = 1.0 / y.sum()
    I = y * Zeq                   # currents for 1 A total
    return Zeq, I


def simulate(geo: CoilGeometry, sim: Optional[SimParams] = None, sweep_freqs: Optional[List[float]] = None,
             progress: Progress = None, keep_cross_section: bool = True) -> SimResult:
    sim = sim or SimParams()
    p = geo.params
    tim = {}
    t0 = time.time()
    if progress:
        progress("3D inductance mesh", 0.0)
    mesh = build_mesh(geo, sim.max_seg_len)
    L = strand_inductance(mesh, geo.n_strands, sim.near_factor, sim.mid_factor)
    tim["inductance_3d"] = time.time() - t0

    def one_freq(freq: float, asm: EddyAssembler):
        t1 = time.time()
        R_via, r_one = via_resistance(geo, sim.rho, freq)
        if sim.include_eddy:
            R, Lax = asm.assemble(freq, progress)
        else:
            R = np.diag(_dc_resistance(geo, sim.rho))
            Lax = np.full_like(L, np.nan)
        w = 2 * math.pi * freq
        Rt = R + np.diag(R_via)
        Z = Rt + 1j * w * L
        Zeq, I = _solve_network(Z)
        return Zeq, I, R, Rt, R_via, r_one, Lax, time.time() - t1

    asm = EddyAssembler(geo, sim)
    Zeq, I, R, Rt, R_via, r_one, Lax, dt = one_freq(sim.frequency, asm)
    tim["eddy_and_network"] = dt
    w = 2 * math.pi * sim.frequency
    S = geo.n_strands
    Rdc_s = _dc_resistance(geo, sim.rho)
    R_dc = 1.0 / np.sum(1.0 / Rdc_s)
    one = np.ones(S) / S
    R_ideal = float(one @ Rt @ one)
    # loss split for 1 A total (peak): P = 1/2 I^H R I
    P_tot = 0.5 * float(np.real(np.conj(I) @ Rt @ I))
    P_dc = 0.5 * float(np.sum(np.abs(I) ** 2 * Rdc_s))
    P_via = 0.5 * float(np.sum(np.abs(I) ** 2 * R_via))
    P_cu = 0.5 * float(np.real(np.conj(I) @ R @ I))
    split = {
        "dc": P_dc / P_tot,
        "skin_and_proximity": (P_cu - P_dc) / P_tot,
        "vias": P_via / P_tot,
        "sharing_penalty_vs_equal_currents": float(Zeq.real) / R_ideal - 1.0,
    }
    L_ax = float(1.0 / np.sum(np.linalg.inv(Lax))) if sim.include_eddy else float("nan")
    res = SimResult(sim.frequency, float(Zeq.imag / w), float(Zeq.real), float(Zeq.imag / Zeq.real),
                    float(R_dc), I, list(geo.strand_ring), L, R, R_via, r_one, R_ideal, L_ax, split, tim,
                    asm.n_patterns)
    if sweep_freqs:
        t2 = time.time()
        for f in sweep_freqs:
            if progress:
                progress("frequency sweep %.3g MHz" % (f / 1e6), 0.0)
            Zf, _, _, _, _, _, _, _ = one_freq(f, asm)
            wf = 2 * math.pi * f
            res.sweep.append((f, float(Zf.imag / wf), float(Zf.real), float(Zf.imag / Zf.real)))
        tim["sweep"] = time.time() - t2
    if keep_cross_section and sim.include_eddy:
        res.cross_section = cross_section_snapshot(geo, sim, asm, I)
    return res


def cross_section_snapshot(geo: CoilGeometry, sim: SimParams, asm: EddyAssembler, I_strand: np.ndarray,
                           theta: float = None) -> dict:
    """Current density in the cross-section at one angle (for the Fig. 6-style plot)."""
    p = geo.params
    S = geo.n_strands
    if theta is None:
        # middle of the first dwell after half a turn
        theta = math.pi
        if geo.zones:
            zs = [z for z in geo.zones if z.t_start > math.pi]
            if len(zs) >= 2:
                theta = 0.5 * (zs[0].t_end + zs[1].t_start)
    key_turns, owners = [], []
    for m in range(asm.n_turn_slots):
        tm = theta + 2 * math.pi * m
        if tm > p.theta_end:
            key_turns.append(None)
            continue
        pos = []
        for s in range(S):
            st = strand_state_at(p, geo.pieces[s], tm)
            if st is None:
                continue
            u, l, f = st
            pos.append(((int(round(u)), l), s))
        pos.sort()
        key_turns.append(tuple(x[0] for x in pos))
        owners += [s for (_, s) in pos]
    key = tuple(key_turns)
    cs = asm._cross_section(key, theta)
    res = solve_cross_section(cs, sim.frequency, sim.rho, sim.eddy_nw, sim.eddy_nt, keep_currents=True)
    Ic = np.array([I_strand[s] for s in owners])
    J = res.current_density(Ic)
    return {
        "theta": theta,
        "r": res.fil_r.tolist(), "z": res.fil_z.tolist(),
        "dw": res.fil_dw.tolist(), "dt": res.fil_dt.tolist(),
        "J_abs": np.abs(J).tolist(),
        "cond_strand": [owners[c] for c in res.fil_cond],
    }
