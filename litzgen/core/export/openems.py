"""openEMS script export (bundle-level model for self-resonance / capacitance).

A full-wave FDTD model cannot resolve 0.2 mm strand gaps on a curved 160 mm
coil at a few MHz in reasonable time (the Cartesian staircase alone would
short neighbouring strands unless the mesh is ~0.05 mm).  The generated
script therefore models each copper layer's bundle as one conducting sheet
following the spiral (the strand detail does not matter for inter-turn
capacitance), joins the layers at the two terminal bars, closes the loop
with a measurement lead above the board and excites it with a lumped port.
It reports Z(f), the low-frequency inductance and the first self-resonance.

The script requires the openEMS Python bindings (CSXCAD + openEMS) and is
emitted as plain text; LitzGen itself does not import openEMS.
"""
from __future__ import annotations

import json
import math

from ..geometry import CoilGeometry, pos_xy_r

TEMPLATE = r'''#!/usr/bin/env python3
"""openEMS bundle-level model of a LitzGen coil (generated file).

Run:  python3 {name}.py            (needs CSXCAD + openEMS python bindings)
Result: {name}_Z.csv (f, Re Z, Im Z) and the first self-resonance printed.
Runtime scales with mesh_res^-3; start with 1.0 mm.
"""
import os, sys, json
import numpy as np
from CSXCAD import ContinuousStructure
from openEMS import openEMS

geo = json.loads({geo_json!r})
mesh_res = float(os.environ.get("LITZ_MESH", geo["mesh_res"]))
f_max = geo["f_max"]
sim_path = os.path.abspath("{name}_sim")

FDTD = openEMS(NrTS=int(geo["nrts"]), EndCriteria=1e-5)
FDTD.SetGaussExcite(f_max / 2, f_max / 2)
FDTD.SetBoundaryCond(["MUR"] * 6)
CSX = ContinuousStructure()
FDTD.SetCSX(CSX)
mesh = CSX.GetGrid()
mesh.SetDeltaUnit(1e-3)

# substrate
fr4 = CSX.AddMaterial("FR4", epsilon=geo["eps_r"])
R = geo["r_board"]
fr4.AddBox([geo["cx"] - R, geo["cy"] - R, geo["z_bot"]], [geo["cx"] + R, geo["cy"] + R, geo["z_top"]], priority=0)

# one conducting sheet per copper layer, following the spiral band
cu = CSX.AddConductingSheet("copper", conductivity=geo["sigma"], thickness=geo["t_cu"])
for z, poly in zip(geo["z_layers"], geo["band"]):
    pts = np.array(poly).T
    cu.AddPolygon(points=pts, norm_dir=2, elevation=z, priority=10)

# terminal bars join all layers; the lead closes the loop above the board
metal = CSX.AddMetal("pec")
for t in geo["terminals"]:
    metal.AddBox([t[0] - 1, t[1] - 1, geo["z_bot"]], [t[0] + 1, t[1] + 1, geo["z_top"]], priority=20)
(ax, ay), (bx, by) = geo["terminals"]
zl = geo["z_top"] + geo["lead_height"]
metal.AddBox([bx - 0.5, by - 0.5, geo["z_top"]], [bx + 0.5, by + 0.5, zl], priority=20)
# horizontal part of the lead: a chain of small boxes along the straight line B -> A
nlead = int(np.hypot(ax - bx, ay - by) / 0.5) + 1
for k in range(nlead + 1):
    px, py = bx + (ax - bx) * k / nlead, by + (ay - by) * k / nlead
    metal.AddBox([px - 0.5, py - 0.5, zl - 0.5], [px + 0.5, py + 0.5, zl + 0.5], priority=20)
gap = 1.0
metal.AddBox([ax - 0.5, ay - 0.5, geo["z_top"] + gap], [ax + 0.5, ay + 0.5, zl], priority=20)
port = FDTD.AddLumpedPort(1, 50, [ax - 0.5, ay - 0.5, geo["z_top"]], [ax + 0.5, ay + 0.5, geo["z_top"] + gap],
                          "z", 1.0, priority=30)

# mesh
m = R + geo["air"]
x = np.concatenate([np.arange(geo["cx"] - m, geo["cx"] + m + 1e-9, mesh_res), [ax, bx]])
y = np.concatenate([np.arange(geo["cy"] - m, geo["cy"] + m + 1e-9, mesh_res), [ay, by]])
z = np.concatenate([geo["z_layers"], [geo["z_bot"], geo["z_top"], geo["z_top"] + gap, zl,
                    geo["z_bot"] - geo["air"], zl + geo["air"]]])
mesh.AddLine("x", np.unique(np.round(x, 4)))
mesh.AddLine("y", np.unique(np.round(y, 4)))
mesh.AddLine("z", np.unique(np.round(z, 4)))
mesh.SmoothMeshLines("z", mesh_res)

if "--dry" not in sys.argv:
    FDTD.Run(sim_path, cleanup=True)
    f = np.linspace(geo["f_min"], f_max, 2001)
    port.CalcPort(sim_path, f)
    Z = port.uf_tot / port.if_tot
    np.savetxt("{name}_Z.csv", np.c_[f, Z.real, Z.imag], delimiter=",", header="f_Hz,ReZ,ImZ")
    k = np.nonzero((Z.imag[:-1] > 0) & (Z.imag[1:] <= 0))[0]
    L0 = Z.imag[1] / (2 * np.pi * f[1])
    print("L(%.2f MHz) = %.3f uH (includes the lead)" % (f[1] / 1e6, L0 * 1e6))
    print("first self-resonance: %s" % (("%.2f MHz" % (f[k[0]] / 1e6)) if len(k) else "above f_max"))
'''


def openems_script(geo: CoilGeometry, name: str = "litz_coil", mesh_res: float = 1.0, f_max: float = 60e6,
                   f_min: float = 1e6, eps_r: float = 4.5, lead_height: float = 10.0, rho: float = 1.72e-8) -> str:
    p = geo.params
    st = p.stackup
    half = p.bundle_width / 2.0
    n = max(8, int(math.degrees(p.theta_end) / 2.0))
    ts = [p.theta_end * k / n for k in range(n + 1)]
    band_xy = [pos_xy_r(p, t, p.radius(t) + half) for t in ts] + [pos_xy_r(p, t, p.radius(t) - half) for t in reversed(ts)]
    # openEMS uses a right-handed frame; flip KiCad's downward y
    band = [[x, -y] for x, y in band_xy]
    z_layers = st.layer_z()
    a = pos_xy_r(p, 0.0, p.radius(0.0))
    b = pos_xy_r(p, p.theta_end, p.radius(p.theta_end))
    geo_d = {
        "cx": p.center_x, "cy": -p.center_y,
        "r_board": p.d_out / 2.0 + 5.0, "air": 40.0,
        "z_layers": z_layers, "z_top": z_layers[0] + st.copper_thickness / 2, "z_bot": z_layers[-1] - st.copper_thickness / 2,
        "band": [band] * st.n_layers, "t_cu": st.copper_thickness * 1e-3, "sigma": 1.0 / rho,
        "terminals": [[a[0], -a[1]], [b[0], -b[1]]], "lead_height": lead_height,
        "eps_r": eps_r, "mesh_res": mesh_res, "f_max": f_max, "f_min": f_min, "nrts": 3e6,
    }
    return TEMPLATE.replace("{geo_json!r}", repr(json.dumps(geo_d))).replace("{name}", name)
