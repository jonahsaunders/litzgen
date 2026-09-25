"""Self-contained HTML report (inline SVG, no plotting library needed)."""
from __future__ import annotations

import html
import math
from typing import List, Optional

import numpy as np

from .geometry import CoilGeometry
from .solver.solve import SimResult

SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
SERIES_DK = ["#3987e5", "#d95926", "#199e70", "#c98500", "#d55181", "#008300", "#9085e9", "#e66767"]
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

CSS = """
.viz-root{color-scheme:light;--surface-1:#fcfcfb;--page:#f9f9f7;--text-primary:#0b0b0b;--text-secondary:#52514e;
--muted:#898781;--grid:#e1e0d9;--axis:#c3c2b7;--border:rgba(11,11,11,.10);%s}
@media (prefers-color-scheme:dark){:root:where(:not([data-theme="light"])) .viz-root{color-scheme:dark;--surface-1:#1a1a19;
--page:#0d0d0d;--text-primary:#fff;--text-secondary:#c3c2b7;--grid:#2c2c2a;--axis:#383835;--border:rgba(255,255,255,.10);%s}}
:root[data-theme="dark"] .viz-root{color-scheme:dark;--surface-1:#1a1a19;--page:#0d0d0d;--text-primary:#fff;--text-secondary:#c3c2b7;
--grid:#2c2c2a;--axis:#383835;--border:rgba(255,255,255,.10);%s}
body{margin:0;background:var(--page);}
.viz-root{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;color:var(--text-primary);background:var(--page);
max-width:1100px;margin:0 auto;padding:24px 16px 48px}
h1{font-size:22px;margin:0 0 4px}h2{font-size:16px;margin:28px 0 8px}
.sub{color:var(--text-secondary);font-size:13px;margin-bottom:16px}
.tiles{display:grid;grid-template-columns:repeat(auto-fit,minmax(140px,1fr));gap:12px}
.tile{background:var(--surface-1);border:1px solid var(--border);border-radius:8px;padding:12px}
.tile .k{font-size:12px;color:var(--text-secondary)}.tile .v{font-size:24px;margin-top:4px}
.card{background:var(--surface-1);border:1px solid var(--border);border-radius:8px;padding:12px;overflow-x:auto}
.grid2{display:grid;grid-template-columns:repeat(auto-fit,minmax(320px,1fr));gap:12px}
table{border-collapse:collapse;font-size:13px;width:100%%}td,th{padding:4px 8px;border-bottom:1px solid var(--grid);text-align:left}
td.n{text-align:right;font-variant-numeric:tabular-nums}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:12px;color:var(--text-secondary);margin:4px 0 8px}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}
svg text{fill:var(--muted);font-size:11px}
.note{font-size:12px;color:var(--text-secondary)}
"""


SEQ_DK = ["#1f2a3a", "#104281", "#1c5cab", "#2a78d6", "#5598e7", "#86b6ef", "#cde2fb"]


def _css() -> str:
    lt = "".join("--series-%d:%s;" % (i + 1, c) for i, c in enumerate(SERIES))
    lt += "".join("--seq-%d:%s;" % (i, c) for i, c in enumerate(SEQ))
    dk = "".join("--series-%d:%s;" % (i + 1, c) for i, c in enumerate(SERIES_DK))
    dk += "".join("--seq-%d:%s;" % (i, c) for i, c in enumerate(SEQ_DK))
    return CSS % (lt, dk, dk)


def _top_view(geo: CoilGeometry, size: int = 520) -> str:
    p = geo.params
    R = p.d_out / 2 + 3
    sc = size / (2 * R)

    def X(x):
        return (x - p.center_x + R) * sc

    def Y(y):
        return (y - p.center_y + R) * sc
    paths = {l: [] for l in range(p.stackup.n_layers)}
    for a in geo.arcs:
        paths[a.layer].append("M%.1f %.1fQ%.1f %.1f %.1f %.1f" % (
            X(a.start[0]), Y(a.start[1]), 2 * X(a.mid[0]) - 0.5 * (X(a.start[0]) + X(a.end[0])),
            2 * Y(a.mid[1]) - 0.5 * (Y(a.start[1]) + Y(a.end[1])), X(a.end[0]), Y(a.end[1])))
    for t in geo.tracks:
        paths[t.layer].append("M%.1f %.1fL%.1f %.1f" % (X(t.start[0]), Y(t.start[1]), X(t.end[0]), Y(t.end[1])))
    out = ['<svg viewBox="0 0 %d %d" width="100%%" style="max-width:%dpx" role="img" aria-label="Coil top view">' % (size, size, size)]
    for l in reversed(range(p.stackup.n_layers)):
        out.append('<path d="%s" fill="none" stroke="var(--series-%d)" stroke-width="%.2f" stroke-opacity="0.9"><title>%s</title></path>'
                   % ("".join(paths[l]), l % 8 + 1, max(0.4, p.strand_width * sc), html.escape(p.stackup.kicad_layers[l])))
    vx = "".join('<circle cx="%.1f" cy="%.1f" r="%.2f"/>' % (X(v.pos[0]), Y(v.pos[1]), max(0.6, v.pad * sc)) for v in geo.vias)
    out.append('<g fill="var(--text-primary)" fill-opacity="0.55">%s</g>' % vx)
    for t in geo.terminals:
        pts = " ".join("%.1f,%.1f" % (X(x), Y(y)) for x, y in t.polygon)
        out.append('<polygon points="%s" fill="var(--text-secondary)"><title>Terminal %s</title></polygon>' % (pts, t.name))
    out.append("</svg>")
    return "".join(out)


def _legend(names: List[str]) -> str:
    return '<div class="legend">%s</div>' % "".join(
        '<span><span class="sw" style="background:var(--series-%d)"></span>%s</span>' % (i % 8 + 1, html.escape(n))
        for i, n in enumerate(names))


def _bars(values: List[float], cats: List[int], labels: List[str], unit: str, ref: Optional[float] = None,
          w: int = 1000, h: int = 260) -> str:
    n = len(values)
    vmax = max(max(values), ref or 0) * 1.1 or 1
    ml, mb, mt = 40, 22, 8
    bw = (w - ml - 8) / n
    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (w, h)]
    for k in range(5):
        v = vmax * k / 4
        y = h - mb - (h - mb - mt) * k / 4
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="var(--grid)"/>' % (ml, w - 4, y, y))
        out.append('<text x="%d" y="%.1f" text-anchor="end">%.2f</text>' % (ml - 4, y + 4, v))
    for i, v in enumerate(values):
        x = ml + i * bw + 1
        bh = (h - mb - mt) * v / vmax
        y = h - mb - bh
        out.append('<path d="M%.1f %.1fv%.1fa2 2 0 0 1 2 -2h%.1fa2 2 0 0 1 2 2v%.1fz" fill="var(--series-%d)"><title>%s: %.3f %s</title></path>'
                   % (x, h - mb, -(bh - 2) if bh > 2 else -bh, max(0.5, bw - 6), (bh - 2) if bh > 2 else bh,
                      cats[i] % 8 + 1, html.escape(labels[i]), v, unit))
        if n <= 24:
            out.append('<text x="%.1f" y="%d" text-anchor="middle">%s</text>' % (x + bw / 2 - 1, h - 6, html.escape(labels[i])))
    if ref:
        y = h - mb - (h - mb - mt) * ref / vmax
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="var(--text-secondary)" stroke-dasharray="4 3"/>' % (ml, w - 4, y, y))
        out.append('<text x="%d" y="%.1f" text-anchor="end">equal share</text>' % (w - 6, y - 4))
    out.append('<line x1="%d" x2="%d" y1="%d" y2="%d" stroke="var(--axis)"/></svg>' % (ml, w - 4, h - mb, h - mb))
    return "".join(out)


def _line(xs, ys, xlabel, ylabel, w=520, h=220, mark_x=None) -> str:
    ml, mb, mt, mr = 56, 30, 10, 22
    x0, x1 = min(xs), max(xs)
    y0, y1 = min(ys), max(ys)
    pad = (y1 - y0) * 0.15 or abs(y1) * 0.05 or 1
    y0 -= pad; y1 += pad

    def X(x):
        return ml + (w - ml - mr) * (x - x0) / ((x1 - x0) or 1)

    def Y(y):
        return h - mb - (h - mb - mt) * (y - y0) / ((y1 - y0) or 1)
    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img">' % (w, h)]
    for k in range(5):
        v = y0 + (y1 - y0) * k / 4
        dig = max(0, int(math.ceil(-math.log10(max((y1 - y0) / 4, 1e-12)))) + 1)
        out.append('<line x1="%d" x2="%d" y1="%.1f" y2="%.1f" stroke="var(--grid)"/><text x="%d" y="%.1f" text-anchor="end">%.*f</text>'
                   % (ml, w - mr, Y(v), Y(v), ml - 4, Y(v) + 4, dig, v))
    for x in xs:
        out.append('<text x="%.1f" y="%d" text-anchor="middle">%.2f</text>' % (X(x), h - 12, x))
    out.append('<text x="%d" y="%d" text-anchor="middle">%s</text>' % ((w + ml) / 2, h, html.escape(xlabel)))
    if mark_x is not None:
        out.append('<line x1="%.1f" x2="%.1f" y1="%d" y2="%d" stroke="var(--axis)" stroke-dasharray="3 3"/>' % (X(mark_x), X(mark_x), mt, h - mb))
    d = "M" + "L".join("%.1f %.1f" % (X(x), Y(y)) for x, y in zip(xs, ys))
    out.append('<path d="%s" fill="none" stroke="var(--series-1)" stroke-width="2"/>' % d)
    for x, y in zip(xs, ys):
        out.append('<circle cx="%.1f" cy="%.1f" r="4" fill="var(--series-1)" stroke="var(--surface-1)" stroke-width="2"><title>%.3f MHz: %.4g %s</title></circle>'
                   % (X(x), Y(y), x, y, html.escape(ylabel)))
    out.append("</svg>")
    return "".join(out)


def _xsec(cs: dict, geo: CoilGeometry, w=1000) -> str:
    """One panel per turn, true aspect ratio, log colour scale (2 decades)."""
    p = geo.params
    r = np.array(cs["r"]) * 1e3; z = np.array(cs["z"]) * 1e3
    dw = np.array(cs["dw"]) * 1e3; dt = np.array(cs["dt"]) * 1e3
    J = np.array(cs["J_abs"])
    th = cs["theta"]
    turns = []
    for m in range(int(math.ceil(p.n_turns))):
        if th + 2 * math.pi * m <= p.theta_end:
            turns.append(p.radius(th + 2 * math.pi * m))
    turns = sorted(turns)
    half_w = p.bundle_width / 2 + 1.5 * p.pitch_lat
    ztop = z.max() + dt.max(); zbot = z.min() - dt.max()
    gap = 14
    ncol = min(3, max(1, len(turns)))
    nrow = int(math.ceil(len(turns) / ncol))
    pw = (w - gap * (ncol + 1)) / ncol
    sc = pw / (2 * half_w)
    ph = (ztop - zbot) * sc
    rowh = ph + 30
    h = int(nrow * rowh + 6)
    jmax = J.max()
    lvl = np.clip(np.log10(np.maximum(J, 1e-30) / jmax) / 2.0 + 1.0, 0, 1)
    out = ['<svg viewBox="0 0 %d %d" width="100%%" role="img" aria-label="Current density per turn">' % (w, h)]
    for k, rc in enumerate(turns):
        x0 = gap + (k % ncol) * (pw + gap)
        y0 = 4 + (k // ncol) * rowh
        out.append('<rect x="%.1f" y="%.1f" width="%.1f" height="%.1f" fill="none" stroke="var(--grid)"/>' % (x0, y0, pw, ph))
        sel = np.nonzero(np.abs(r - rc) < half_w)[0]
        for i in sel:
            c = "var(--seq-%d)" % int(round(lvl[i] * (len(SEQ) - 1)))
            x = x0 + (r[i] - dw[i] / 2 - (rc - half_w)) * sc
            y = y0 + (ztop - (z[i] + dt[i] / 2)) * sc
            out.append('<rect x="%.2f" y="%.2f" width="%.2f" height="%.2f" fill="%s"><title>strand %d</title></rect>'
                       % (x, y, dw[i] * sc + 0.25, dt[i] * sc + 0.25, c, cs["cond_strand"][i] + 1))
        out.append('<text x="%.1f" y="%.1f" text-anchor="middle">turn at r = %.1f mm</text>' % (x0 + pw / 2, y0 + ph + 16, rc))
    out.append("</svg>")
    leg = '<div class="legend">|J| (log scale): ' + "".join('<span class="sw" style="background:var(--seq-%d)"></span>' % k for k in range(len(SEQ))) + \
          ' 1 %% &rarr; 100 %% of max (%.3g A/mm&sup2; per 1 A coil current); inner radius on the left</div>' % (jmax * 1e-6)
    return leg + "".join(out)


def html_report(geo: CoilGeometry, res: Optional[SimResult], drc_violations: int = 0, title: str = "LitzGen coil report") -> str:
    p = geo.params
    s = geo.summary()
    tiles = []
    if res is not None:
        tiles += [("Inductance", "%.3f &micro;H" % (res.L * 1e6)), ("ESR @ %.2f MHz" % (res.frequency / 1e6), "%.0f m&Omega;" % (res.R * 1e3)),
                  ("Q", "%.0f" % res.Q), ("R<sub>ac</sub>/R<sub>dc</sub>", "%.1f" % (res.R / res.R_dc))]
    tiles += [("Strands", str(geo.n_strands)), ("Vias", str(s["vias_total"])), ("DRC violations", str(drc_violations))]
    H = ['<!doctype html><html><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
         '<title>%s</title><style>%s</style></head><body><div class="viz-root">' % (html.escape(title), _css()),
         '<h1>%s</h1><div class="sub">%s &middot; d<sub>out</sub> %.1f mm &middot; %g turns &middot; %d&times;%d strand matrix</div>'
         % (html.escape(title), html.escape(p.group_name), p.d_out, p.n_turns, p.n_cols, p.stackup.n_layers),
         '<div class="tiles">%s</div>' % "".join('<div class="tile"><div class="k">%s</div><div class="v">%s</div></div>' % t for t in tiles)]
    H.append('<h2>Layout</h2><div class="card">%s%s</div>' % (_legend(list(p.stackup.kicad_layers)), _top_view(geo)))
    if res is not None:
        S = geo.n_strands
        mags = [abs(x) * S for x in res.strand_currents]
        rings = list(res.ring_of_strand)
        rnames = sorted(set(rings))
        H.append('<h2>Strand current sharing</h2><div class="card">%s%s<div class="note">Current magnitude per strand relative to an equal 1/%d share. '
                 'Ring currents: %s.</div></div>' % (
                     _legend(["ring %d (%s)" % (r, "outer bundle" if r == 0 else "inner bundle" if r == 1 else "ring") for r in rnames]),
                     _bars(mags, [rnames.index(r) for r in rings], [str(i + 1) for i in range(S)], "x equal share", ref=1.0), S,
                     ", ".join("ring %s %.1f %%" % (k, abs(v) * 100) for k, v in sorted(res.ring_share().items()))))
        if res.cross_section:
            H.append('<h2>Cross-section current density</h2><div class="card">%s<div class="note">Axisymmetric eddy-current solution at %.0f&deg; '
                     'with the solved strand currents (compare Fig. 6 of Kale &amp; Wicht).</div></div>'
                     % (_xsec(res.cross_section, geo), math.degrees(res.cross_section["theta"]) % 360))
        if res.sweep:
            pts = sorted(res.sweep + [(res.frequency, res.L, res.R, res.Q)])
            f = [x[0] / 1e6 for x in pts]
            H.append('<h2>Frequency response</h2><div class="grid2"><div class="card"><div class="k">ESR (&Omega;)</div>%s</div>'
                     '<div class="card"><div class="k">Q</div>%s</div></div><div class="note">Quasi-static model: inductance is '
                     'frequency-independent here and self-resonance is not included (see the openEMS export).</div>'
                     % (_line(f, [x[2] for x in pts], "MHz", "ohm", mark_x=res.frequency / 1e6),
                        _line(f, [x[3] for x in pts], "MHz", "Q", mark_x=res.frequency / 1e6)))
        ls = res.loss_split
        H.append('<h2>Loss breakdown</h2><div class="card"><table><tr><th>Part</th><th>Share of coil loss</th></tr>'
                 '<tr><td>DC resistance of copper</td><td class="n">%.1f %%</td></tr>'
                 '<tr><td>Skin + proximity (eddy) excess</td><td class="n">%.1f %%</td></tr>'
                 '<tr><td>Vias</td><td class="n">%.1f %%</td></tr></table>'
                 '<div class="note">ESR with the solved sharing is %+.1f %% versus forcing equal strand currents '
                 '(%.0f m&Omega;). Ring-model inductance cross-check: %.3f &micro;H.</div></div>'
                 % (ls["dc"] * 100, ls["skin_and_proximity"] * 100, ls["vias"] * 100,
                    ls["sharing_penalty_vs_equal_currents"] * 100, res.R_ideal_sharing * 1e3, res.L_axisym * 1e6))
    rows = "".join("<tr><td>%s</td><td class='n'>%d</td></tr>" % (html.escape(k), v) for k, v in s["vias_by_type"].items())
    H.append('<h2>Geometry</h2><div class="grid2"><div class="card"><table><tr><th>Via type</th><th>Count</th></tr>%s</table></div>'
             '<div class="card"><table>'
             '<tr><td>Turn pitch / gap</td><td class="n">%.2f / %.2f mm</td></tr>'
             '<tr><td>Inner diameter</td><td class="n">%.1f mm</td></tr>'
             '<tr><td>Bundle width</td><td class="n">%.2f mm</td></tr>'
             '<tr><td>Transposition zones</td><td class="n">%d</td></tr>'
             '<tr><td>Strand length (min / max)</td><td class="n">%.1f / %.1f mm</td></tr></table></div></div>'
             % (rows, p.turn_pitch, p.effective_turn_gap, p.derived_d_in, p.bundle_width, s["zones"],
                s["strand_length_mm"]["min"], s["strand_length_mm"]["max"]))
    if geo.warnings:
        H.append('<h2>Warnings</h2><div class="card"><ul>%s</ul></div>' % "".join("<li>%s</li>" % html.escape(w) for w in geo.warnings))
    H.append('<h2>Parameters</h2><div class="card"><pre style="font-size:11px;white-space:pre-wrap">%s</pre></div>' % html.escape(p.to_json()))
    H.append("</div></body></html>")
    return "".join(H)
