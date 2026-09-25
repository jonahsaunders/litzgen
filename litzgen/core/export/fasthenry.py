"""FastHenry (.inp) export for independent quasi-static cross-checks.

Every strand becomes a chain of rectangular segments (sub-divided into
nwinc x nhinc filaments by FastHenry for skin / proximity effects), vias
become vertical bars of equal perimeter, the start nodes of all strands are
shorted with ``.equiv`` (terminal A) and likewise the end nodes (terminal B),
and one port ``.external`` is defined between the two terminals.  Works with
FastHenry 3.0wr and FastHenry2 (FastFieldSolvers).
"""
from __future__ import annotations

import math
from typing import List, Optional

from ..geometry import CoilGeometry, Dwell, Jog, ViaHop, lateral_offset, phi_of
from ..params import SimParams


def fasthenry_text(geo: CoilGeometry, sim: Optional[SimParams] = None, max_seg_len: float = 3.0,
                   nwinc: int = 7, nhinc: int = 3, fmin: float = None, fmax: float = None,
                   ndec: int = 1) -> str:
    sim = sim or SimParams()
    p = geo.params
    z = p.stackup.layer_z()
    sigma_mm = 1.0 / sim.rho * 1e-3          # S/mm
    fmin = fmin or sim.frequency
    fmax = fmax or sim.frequency
    out: List[str] = []
    out.append("* LitzGen export: %d strands, %d vias, d_out=%.2f mm, %.3g turns"
               % (geo.n_strands, len(geo.vias), p.d_out, p.n_turns))
    out.append("* strand %.3f x %.3f mm, filaments per segment %d x %d" % (p.strand_width, p.stackup.copper_thickness, nwinc, nhinc))
    out.append(".units mm")
    out.append(".default sigma=%.6g" % sigma_mm)
    via_side = math.pi * p.via.drill / 4.0   # perimeter-equivalent square bar
    starts, ends = [], []

    def xyz(theta, u, layer):
        r = p.radius(theta) + lateral_offset(p, u)
        ph = phi_of(p, theta)
        return (p.center_x + r * math.cos(ph), p.center_y - r * math.sin(ph), z[layer])

    for s, pcs in enumerate(geo.pieces):
        nodes = []
        segs = []   # (i0, i1, kind)
        cur = None

        def add(pt):
            nodes.append(pt)
            return len(nodes) - 1

        for pc in pcs:
            if isinstance(pc, (Dwell, Jog)):
                ua, ub = (pc.u, pc.u) if isinstance(pc, Dwell) else (pc.ua, pc.ub)
                rr = p.radius(0.5 * (pc.ta + pc.tb)) + lateral_offset(p, 0.5 * (ua + ub))
                nseg = max(1, int(math.ceil(abs(pc.tb - pc.ta) * rr / max_seg_len)))
                for k in range(nseg + 1):
                    f = k / nseg
                    pt = xyz(pc.ta + f * (pc.tb - pc.ta), ua + f * (ub - ua), pc.layer)
                    if cur is not None and math.dist(nodes[cur], pt) < 1e-6:
                        continue
                    i = add(pt)
                    if cur is not None:
                        segs.append((cur, i, "t"))
                    cur = i
            elif isinstance(pc, ViaHop):
                pt = xyz(pc.t, pc.u, pc.lb)
                i = add(pt)
                segs.append((cur, i, "v"))
                cur = i
        for k, (x, y, zz) in enumerate(nodes):
            out.append("N%d_%d x=%.6f y=%.6f z=%.6f" % (s, k, x, y, zz))
        for k, (a, b, kind) in enumerate(segs):
            if kind == "t":
                out.append("E%d_%d N%d_%d N%d_%d w=%.6f h=%.6f nwinc=%d nhinc=%d rw=2 rh=2"
                           % (s, k, s, a, s, b, p.strand_width, p.stackup.copper_thickness, nwinc, nhinc))
            else:
                out.append("E%d_%d N%d_%d N%d_%d w=%.6f h=%.6f wx=1 wy=0 wz=0 nwinc=3 nhinc=3"
                           % (s, k, s, a, s, b, via_side, via_side))
        starts.append("N%d_0" % s)
        ends.append("N%d_%d" % (s, len(nodes) - 1))
    out.append(".equiv " + " ".join(starts))
    out.append(".equiv " + " ".join(ends))
    out.append(".external %s %s coil" % (starts[0], ends[0]))
    out.append(".freq fmin=%.6g fmax=%.6g ndec=%d" % (fmin, fmax, ndec))
    out.append(".end")
    return "\n".join(out) + "\n"
