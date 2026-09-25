"""Strand path generation and mapping onto an Archimedean spiral.

Every strand is described in *path space* ``(theta, u, layer)``: ``theta`` is
the unwrapped spiral angle (0 at the outer terminal, 2*pi*n at the inner
one), ``u`` the lateral slot index across the bundle and ``layer`` the copper
layer.  The path is a list of pieces:

* ``Dwell``  constant (u, layer) over [theta_a, theta_b]
* ``Jog``    linear change of u on one layer over [theta_a, theta_b]
* ``ViaHop`` layer change at one theta (an adjacent-layer via)

Mapping to board coordinates only happens at the very end, so the same
path description feeds the KiCad writer, the DRC checker and the solver.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

import numpy as np

from .params import CoilParams
from .transposition import Slot, build_rings, initial_slots, step_plan


# ----------------------------------------------------------------------------- pieces
@dataclass
class Dwell:
    u: float
    layer: int
    ta: float
    tb: float


@dataclass
class Jog:
    ua: float
    ub: float
    layer: int
    ta: float
    tb: float


@dataclass
class ViaHop:
    u: float
    la: int
    lb: int
    t: float


Piece = object


@dataclass
class Zone:
    """One transposition zone (angles in rad, unwrapped)."""
    t_start: float
    t_x_end: float
    t_vias: List[float]
    t_r_start: float
    t_end: float


# ----------------------------------------------------------------------------- output primitives
@dataclass
class Track:
    layer: int
    start: Tuple[float, float]
    end: Tuple[float, float]
    width: float
    strand: int


@dataclass
class ArcTrack:
    layer: int
    start: Tuple[float, float]
    mid: Tuple[float, float]
    end: Tuple[float, float]
    width: float
    strand: int


@dataclass
class Via:
    pos: Tuple[float, float]
    la: int
    lb: int
    drill: float
    pad: float
    strand: int


@dataclass
class Terminal:
    name: str
    polygon: List[Tuple[float, float]]
    holes: List[Tuple[float, float]]
    anchor: Tuple[float, float]
    drill: float
    hole_pad: float


@dataclass
class CoilGeometry:
    params: CoilParams
    strand_slots0: List[Slot]
    strand_ring: List[int]
    pieces: List[List[Piece]]
    zones: List[Zone]
    tracks: List[Track] = field(default_factory=list)
    arcs: List[ArcTrack] = field(default_factory=list)
    vias: List[Via] = field(default_factory=list)
    terminals: List[Terminal] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)

    @property
    def n_strands(self) -> int:
        return len(self.pieces)

    # --------------------------------------------------------------- statistics
    def via_counts(self) -> Dict[str, int]:
        n = self.params.stackup.n_layers
        out: Dict[str, int] = {}
        for v in self.vias:
            lo, hi = sorted((v.la, v.lb))
            kind = "blind" if (lo == 0 or hi == n - 1) else "buried"
            key = "%s %s-%s" % (kind, self.params.stackup.kicad_layers[lo], self.params.stackup.kicad_layers[hi])
            out[key] = out.get(key, 0) + 1
        return out

    def vias_per_ring(self) -> Dict[int, int]:
        out: Dict[int, int] = {}
        for v in self.vias:
            r = self.strand_ring[v.strand]
            out[r] = out.get(r, 0) + 1
        return out

    def strand_lengths(self) -> np.ndarray:
        """Physical length of every strand (mm), including via barrels."""
        p = self.params
        z = p.stackup.layer_z()
        L = np.zeros(self.n_strands)
        for s, pcs in enumerate(self.pieces):
            for pc in pcs:
                if isinstance(pc, Dwell):
                    L[s] += _dwell_length(p, pc)
                elif isinstance(pc, Jog):
                    a = pos_xy(p, pc.ta, pc.ua)
                    b = pos_xy(p, pc.tb, pc.ub)
                    L[s] += math.hypot(b[0] - a[0], b[1] - a[1])
                else:
                    L[s] += abs(z[pc.lb] - z[pc.la])
        return L

    def exposure(self) -> np.ndarray:
        """Fraction of each strand's dwell length spent in each slot: array [strand, u, layer]."""
        p = self.params
        C, R = p.n_cols, p.stackup.n_layers
        E = np.zeros((self.n_strands, C + 2, R))
        for s, pcs in enumerate(self.pieces):
            for pc in pcs:
                if isinstance(pc, Dwell):
                    E[s, int(round(pc.u)) + 1, pc.layer] += _dwell_length(p, pc)
        E /= E.sum(axis=(1, 2), keepdims=True)
        return E[:, 1:C + 1, :]

    def summary(self) -> dict:
        p = self.params
        L = self.strand_lengths()
        return {
            "strands": self.n_strands,
            "rings": sorted(set(self.strand_ring)),
            "zones": len(self.zones),
            "vias_total": len(self.vias),
            "vias_by_type": self.via_counts(),
            "vias_by_ring": {str(k): v for k, v in sorted(self.vias_per_ring().items())},
            "strand_length_mm": {"min": float(L.min()), "max": float(L.max()), "mean": float(L.mean())},
            "derived_d_in_mm": p.derived_d_in,
            "turn_pitch_mm": p.turn_pitch,
            "turn_gap_mm": p.effective_turn_gap,
            "tracks": len(self.tracks),
            "arcs": len(self.arcs),
        }


# ----------------------------------------------------------------------------- mapping helpers
def lateral_offset(p: CoilParams, u: float) -> float:
    return (u - (p.n_cols - 1) / 2.0) * p.pitch_lat


def phi_of(p: CoilParams, theta: float) -> float:
    d = -1.0 if p.clockwise else 1.0
    return math.radians(p.start_angle_deg) + d * theta


def pos_xy(p: CoilParams, theta: float, u: float) -> Tuple[float, float]:
    """Board coordinates (KiCad convention: y grows downward)."""
    r = p.radius(theta) + lateral_offset(p, u)
    ph = phi_of(p, theta)
    return (p.center_x + r * math.cos(ph), p.center_y - r * math.sin(ph))


def pos_xy_r(p: CoilParams, theta: float, r: float) -> Tuple[float, float]:
    ph = phi_of(p, theta)
    return (p.center_x + r * math.cos(ph), p.center_y - r * math.sin(ph))


def _dwell_length(p: CoilParams, d: Dwell) -> float:
    # arc length of an Archimedean spiral offset: integrate sqrt(r^2 + (dr/dth)^2)
    n = max(2, int(abs(d.tb - d.ta) / 0.05) + 2)
    th = np.linspace(d.ta, d.tb, n)
    r = p.r_start - p.turn_pitch * th / (2 * math.pi) + lateral_offset(p, d.u)
    drdt = p.turn_pitch / (2 * math.pi)
    f = np.sqrt(r * r + drdt * drdt)
    return float(np.sum((f[1:] + f[:-1]) * 0.5 * np.diff(th)))


# ----------------------------------------------------------------------------- zone layout
@dataclass
class ZoneLengths:
    jog: float
    via_gap: float
    via_pitch: float
    n_via: int

    @property
    def total(self) -> float:
        return 2 * self.jog + 2 * self.via_gap + max(0, self.n_via - 1) * self.via_pitch


def zone_lengths(p: CoilParams) -> ZoneLengths:
    pl = p.pitch_lat
    # Parallel jogs of neighbouring strands are p*cos(alpha) apart; a jog ending next to a
    # neighbour's straight run sees (w + w_jog)/2.  Jogs are emitted as three-point arcs that
    # follow the spiral, so chord sag does not eat into the clearance.
    need = (p.strand_width + p.w_jog) / 2.0 + p.clearance + 0.003
    cos_a = need / pl
    if cos_a >= 1.0:
        raise ValueError("no jog angle satisfies the clearance: increase strand_gap or reduce jog_width")
    # the spiral itself drifts inward at angle psi to the tangent; an inward jog adds to it
    r_min = max(1.0, p.r_end - p.bundle_width / 2.0 - pl)
    psi = math.atan(p.turn_pitch / (2.0 * math.pi * r_min))
    alpha = math.acos(cos_a) - psi
    if alpha <= 0.05:
        raise ValueError("no jog angle satisfies the clearance: increase strand_gap or reduce jog_width")
    jog = 1.03 * pl / math.tan(alpha)
    # consecutive vias in one lane share a layer: the trace arriving at via i+1 and the trace
    # leaving via i end in round caps at the via centres, so the pitch is set by the trace width
    via_pitch = 1.02 * max(p.strand_width + p.clearance, p.via.pad + p.clearance,
                           p.via.pad / 2 + p.strand_width / 2 + p.clearance, p.via.drill + p.hole_clearance)
    via_gap = p.via.pad / 2.0 + 0.05
    _, n_via = step_plan(p.n_cols, p.stackup.n_layers)
    return ZoneLengths(jog, via_gap, via_pitch, n_via)


def max_steps_per_turn(p: CoilParams) -> int:
    """Largest steps_per_turn whose zones still fit at the innermost radius."""
    zl = zone_lengths(p)
    r_min = p.r_end + lateral_offset(p, -1.0) - p.pitch_lat * 0.5
    if r_min <= 0:
        return 0
    return max(0, int(2 * math.pi / ((zl.total + 0.5) / r_min) * 0.999))


def layout_zones(p: CoilParams) -> Tuple[List[Zone], List[str]]:
    zl = zone_lengths(p)
    warn: List[str] = []
    N = int(p.steps_per_turn)
    if N < 0:
        N = max_steps_per_turn(p)
        warn.append("auto transposition density: %d steps per turn" % N)
    if N <= 0:
        return [], warn
    dth = 2 * math.pi / N
    margin = math.radians(p.terminal_margin_deg)
    t_end = p.theta_end
    zones: List[Zone] = []
    i = 0
    min_dwell_mm = 0.5
    while True:
        base = (i + 0.5) * dth
        turn = int(base // (2 * math.pi))
        if p.stagger_alternate_turns and turn % 2 == 1:
            base += dth / 2.0
        i += 1
        if base > t_end + dth:
            break
        r_eff = p.radius(base) + lateral_offset(p, -1.0) - p.pitch_lat * 0.5
        if r_eff <= 0:
            break
        span = zl.total / r_eff
        ta = base - span / 2.0
        tb = base + span / 2.0
        if ta < margin or tb > t_end - margin:
            continue
        if span > dth - min_dwell_mm / r_eff:
            raise ValueError(
                "transposition zones do not fit at radius %.1f mm: zone needs %.1f deg, spacing is %.1f deg. "
                "Reduce steps_per_turn (max ~%d here) or the jog length." %
                (r_eff, math.degrees(span), math.degrees(dth), int(2 * math.pi / (span + min_dwell_mm / r_eff))))
        tx = ta + zl.jog / r_eff
        tv = [tx + (zl.via_gap + k * zl.via_pitch) / r_eff for k in range(zl.n_via)]
        tr = tx + (2 * zl.via_gap + max(0, zl.n_via - 1) * zl.via_pitch) / r_eff
        zones.append(Zone(ta, tx, tv, tr, tr + zl.jog / r_eff))
    zones.sort(key=lambda z: z.t_start)
    # staggering can bring the last zone of one turn close to the first of the next: drop overlaps
    kept: List[Zone] = []
    for z in zones:
        if kept:
            r_eff = p.radius(z.t_start) + lateral_offset(p, -1.0) - p.pitch_lat * 0.5
            if z.t_start < kept[-1].t_end + min_dwell_mm / r_eff:
                continue
        kept.append(z)
    return kept, warn


# ----------------------------------------------------------------------------- path construction
def _append_dwell(pcs: List[Piece], u: float, layer: int, ta: float, tb: float) -> None:
    if tb <= ta + 1e-12:
        return
    if pcs and isinstance(pcs[-1], Dwell) and pcs[-1].u == u and pcs[-1].layer == layer and abs(pcs[-1].tb - ta) < 1e-12:
        pcs[-1].tb = tb
    else:
        pcs.append(Dwell(u, layer, ta, tb))


def build_paths(p: CoilParams) -> Tuple[List[List[Piece]], List[Slot], List[int], List[Zone], List[str]]:
    C, R = p.n_cols, p.stackup.n_layers
    plan, n_via = step_plan(C, R)
    slots0, ring_of = initial_slots(C, R)
    zones, warn = layout_zones(p)
    pieces: List[List[Piece]] = [[] for _ in slots0]
    cur = list(slots0)
    t_prev = 0.0
    for z in zones:
        for s, (u, l) in enumerate(cur):
            m = plan[(u, l)]
            pcs = pieces[s]
            _append_dwell(pcs, u, l, t_prev, z.t_start)
            u1 = u + m.x_du
            if m.x_du:
                pcs.append(Jog(u, u1, l, z.t_start, z.t_x_end))
            else:
                _append_dwell(pcs, u, l, z.t_start, z.t_x_end)
            l1 = l
            if m.via_index is not None:
                tv = z.t_vias[m.via_index]
                _append_dwell(pcs, u1, l, z.t_x_end, tv)
                l1 = l + m.dl
                pcs.append(ViaHop(u1, l, l1, tv))
                _append_dwell(pcs, u1, l1, tv, z.t_r_start)
            else:
                _append_dwell(pcs, u1, l, z.t_x_end, z.t_r_start)
            u2 = u1 + m.r_du
            if m.r_du:
                pcs.append(Jog(u1, u2, l1, z.t_r_start, z.t_end))
            else:
                _append_dwell(pcs, u1, l1, z.t_r_start, z.t_end)
            cur[s] = (u2, l1)
        t_prev = z.t_end
    for s, (u, l) in enumerate(cur):
        _append_dwell(pieces[s], u, l, t_prev, p.theta_end)
    return pieces, slots0, ring_of, zones, warn


# ----------------------------------------------------------------------------- primitive emission
def _emit_dwell(p: CoilParams, d: Dwell, s: int, out_arcs: List[ArcTrack]) -> None:
    span = d.tb - d.ta
    n = max(1, int(math.ceil(abs(span) / math.radians(p.max_arc_deg))))
    for k in range(n):
        a = d.ta + span * k / n
        b = d.ta + span * (k + 1) / n
        out_arcs.append(ArcTrack(d.layer, pos_xy(p, a, d.u), pos_xy(p, (a + b) / 2, d.u), pos_xy(p, b, d.u),
                                 p.strand_width, s))


def _terminal(p: CoilParams, name: str, ta: float, tb: float) -> Terminal:
    half = p.bundle_width / 2.0
    n = max(3, int(abs(tb - ta) / math.radians(1.0)) + 2)
    ts = np.linspace(ta, tb, n)
    outer = [pos_xy_r(p, t, p.radius(t) + half) for t in ts]
    inner = [pos_xy_r(p, t, p.radius(t) - half) for t in ts[::-1]]
    poly = outer + inner
    tm = 0.5 * (ta + tb)
    anchor = pos_xy_r(p, tm, p.radius(tm))
    # stitching holes along the bar centre-line
    hole_pad = p.terminal_drill + 0.5
    arc_len = abs(tb - ta) * p.radius(tm)
    nh = max(1, int((arc_len - hole_pad) // (hole_pad + p.clearance)) + 1)
    nh = min(nh, 6)
    holes = []
    rows = 2 if p.bundle_width >= 2 * hole_pad + p.clearance else 1
    for j in range(rows):
        rr = p.radius(tm) + (0 if rows == 1 else (-1 if j == 0 else 1) * (p.bundle_width / 4.0))
        for i in range(nh):
            f = (i + 0.5) / nh
            t = ta + (tb - ta) * (0.15 + 0.7 * f)
            holes.append(pos_xy_r(p, t, rr))
    return Terminal(name, poly, holes, anchor, p.terminal_drill, hole_pad)


def generate(p: CoilParams) -> CoilGeometry:
    """Build the complete coil: paths, arcs/tracks/vias, terminals."""
    warn = p.validate()
    pieces, slots0, ring_of, zones, w2 = build_paths(p)
    if p.omit_rings:
        keep = [i for i, r in enumerate(ring_of) if r not in set(p.omit_rings)]
        if not keep:
            raise ValueError("omit_rings removes every strand")
        pieces = [pieces[i] for i in keep]
        slots0 = [slots0[i] for i in keep]
        ring_of = [ring_of[i] for i in keep]
    geo = CoilGeometry(p, slots0, ring_of, pieces, zones, warnings=warn + w2)
    if len(set(ring_of)) and -1 in ring_of:
        geo.warnings.append("%d centre strand(s) cannot be transposed with this matrix size and run straight"
                            % ring_of.count(-1))
    z = p.stackup.layer_z()
    for s, pcs in enumerate(pieces):
        for pc in pcs:
            if isinstance(pc, Dwell):
                _emit_dwell(p, pc, s, geo.arcs)
            elif isinstance(pc, Jog):
                tm = 0.5 * (pc.ta + pc.tb)
                geo.arcs.append(ArcTrack(pc.layer, pos_xy(p, pc.ta, pc.ua), pos_xy(p, tm, 0.5 * (pc.ua + pc.ub)),
                                         pos_xy(p, pc.tb, pc.ub), p.w_jog, s))
            elif isinstance(pc, ViaHop):
                geo.vias.append(Via(pos_xy(p, pc.t, pc.u), pc.la, pc.lb, p.via.drill, p.via.pad, s))
    if not p.per_strand_nets:
        eps_o = 0.6 / p.r_start
        T_o = p.terminal_length / p.r_start
        eps_i = 0.6 / p.r_end
        T_i = p.terminal_length / p.r_end
        geo.terminals.append(_terminal(p, "A", -T_o, eps_o))
        geo.terminals.append(_terminal(p, "B", p.theta_end - eps_i, p.theta_end + T_i))
    if len(zones) == 0:
        geo.warnings.append("no transposition zones: strands run untransposed")
    return geo


# ----------------------------------------------------------------------------- path queries used by the solver
def strand_state_at(p: CoilParams, pcs: List[Piece], theta: float) -> Optional[Tuple[float, int, float]]:
    """(u, layer, length factor) of a strand at unwrapped angle theta; None outside the coil.

    The length factor is the ratio of physical length to the length of a
    dwell at the same angle (> 1 inside jogs).
    """
    lo, hi = 0, len(pcs) - 1
    # pieces are ordered by angle; binary search on start angle
    while lo < hi:
        mid = (lo + hi + 1) // 2
        pm = pcs[mid]
        ta = pm.t if isinstance(pm, ViaHop) else pm.ta
        if ta <= theta:
            lo = mid
        else:
            hi = mid - 1
    pc = pcs[lo]
    if isinstance(pc, ViaHop):
        pc = pcs[min(lo + 1, len(pcs) - 1)]
    if isinstance(pc, Dwell):
        if theta < pc.ta - 1e-12 or theta > pc.tb + 1e-12:
            return None
        return pc.u, pc.layer, 1.0
    if isinstance(pc, Jog):
        f = (theta - pc.ta) / (pc.tb - pc.ta)
        u = pc.ua + f * (pc.ub - pc.ua)
        ds = (p.radius(theta) + lateral_offset(p, u)) * (pc.tb - pc.ta)
        du = (pc.ub - pc.ua) * p.pitch_lat
        return u, pc.layer, math.sqrt(1.0 + (du / ds) ** 2)
    return None
