"""Command-line interface: generate, check, simulate, sweep and export without KiCad.

Examples
  python -m litzgen params > coil.json                     # paper defaults
  python -m litzgen generate -p coil.json -o coil.kicad_pcb --drc
  python -m litzgen simulate -p coil.json --report coil.html --json coil.json.out --sweep
  python -m litzgen sweep --set n_cols=4 --vary steps_per_turn=12,18,24 --quality fast
  python -m litzgen export fasthenry -o coil.inp
"""
from __future__ import annotations

import argparse
import json
import sys
import time

from .core.drc import check
from .core.export.fasthenry import fasthenry_text
from .core.export.kicad_sexpr import board_text
from .core.export.openems import openems_script
from .core.geometry import generate
from .core.params import CoilParams, SimParams


def _coerce(v: str):
    try:
        return json.loads(v)
    except Exception:
        return v


def apply_sets(d: dict, sets):
    for s in sets or []:
        k, v = s.split("=", 1)
        cur = d
        parts = k.split(".")
        for part in parts[:-1]:
            cur = cur[part]
        if parts[-1] not in cur:
            raise SystemExit("unknown parameter: %s" % k)
        cur[parts[-1]] = _coerce(v)
    return d


def load_params(args) -> CoilParams:
    d = CoilParams().to_dict()
    if getattr(args, "params", None):
        with open(args.params) as f:
            d.update(json.load(f))
    return CoilParams.from_dict(apply_sets(d, getattr(args, "set", None)))


def _progress(stage, frac):
    sys.stderr.write("\r%-40s %3.0f%%" % (stage, frac * 100))
    sys.stderr.flush()


def cmd_params(args):
    print(load_params(args).to_json())


def cmd_generate(args):
    p = load_params(args)
    g = generate(p)
    print(json.dumps(g.summary(), indent=2))
    for w in g.warnings:
        print("warning:", w)
    if args.drc:
        v = check(g)
        print("DRC: %d violation(s)" % len(v))
        for x in v[:20]:
            print("  ", x)
        if v and args.strict:
            raise SystemExit(2)
    if args.out:
        with open(args.out, "w") as f:
            f.write(board_text(g))
        print("wrote", args.out)


def cmd_simulate(args):
    from .core.report import html_report
    from .core.solver.solve import simulate
    p = load_params(args)
    g = generate(p)
    sim = SimParams.preset(args.quality, frequency=args.freq, temperature_c=args.temp)
    sweep = None
    if args.sweep:
        f0 = args.freq
        sweep = [f0 * k for k in (0.94, 0.97, 1.03, 1.06)]
    t = time.time()
    res = simulate(g, sim, sweep_freqs=sweep, progress=None if args.quiet else _progress)
    sys.stderr.write("\n")
    d = res.to_dict()
    d["geometry"] = g.summary()
    d["sim_params"] = sim.to_dict()
    print("L = %.4f uH   ESR = %.1f mOhm   Q = %.1f   Rac/Rdc = %.2f   (%.0f s)" %
          (res.L * 1e6, res.R * 1e3, res.Q, res.R / res.R_dc, time.time() - t))
    print("ring current share:", ", ".join("ring %s: %.1f %%" % (k, abs(v) * 100) for k, v in sorted(res.ring_share().items())))
    if args.json:
        with open(args.json, "w") as f:
            json.dump(d, f, indent=2)
    if args.report:
        v = check(g)
        with open(args.report, "w") as f:
            f.write(html_report(g, res, len(v)))
        print("wrote", args.report)


def cmd_sweep(args):
    from .core.solver.solve import simulate
    key, vals = args.vary.split("=", 1)
    vals = [_coerce(x) for x in vals.split(",")]
    print("%-20s %8s %9s %7s %7s %6s" % (key, "L_uH", "ESR_mOhm", "Q", "vias", "DRC"))
    for v in vals:
        a = argparse.Namespace(params=args.params, set=(args.set or []) + ["%s=%s" % (key, json.dumps(v))])
        p = load_params(a)
        try:
            g = generate(p)
        except ValueError as e:
            print("%-20s rejected: %s" % (v, e))
            continue
        nv = len(check(g))
        if args.no_sim:
            print("%-20s %8s %9s %7s %7d %6d" % (v, "-", "-", "-", len(g.vias), nv))
            continue
        r = simulate(g, SimParams.preset(args.quality, frequency=args.freq), keep_cross_section=False)
        print("%-20s %8.4f %9.1f %7.1f %7d %6d" % (v, r.L * 1e6, r.R * 1e3, r.Q, len(g.vias), nv))


def cmd_export(args):
    p = load_params(args)
    g = generate(p)
    if args.kind == "fasthenry":
        txt = fasthenry_text(g, SimParams(frequency=args.freq), nwinc=args.nwinc, nhinc=args.nhinc)
    elif args.kind == "openems":
        txt = openems_script(g, name=args.name)
    else:
        txt = board_text(g)
    with open(args.out, "w") as f:
        f.write(txt)
    print("wrote", args.out)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="litzgen", description="Multi-bundle PCB Litz coil generator and solver")
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("-p", "--params", help="JSON parameter file (missing keys use paper defaults)")
        sp.add_argument("--set", action="append", metavar="KEY=VALUE", help="override, e.g. --set stackup.copper_thickness=0.035")

    sp = sub.add_parser("params", help="print parameters as JSON"); common(sp); sp.set_defaults(fn=cmd_params)
    sp = sub.add_parser("generate", help="generate geometry, optional DRC and .kicad_pcb"); common(sp)
    sp.add_argument("-o", "--out"); sp.add_argument("--drc", action="store_true"); sp.add_argument("--strict", action="store_true")
    sp.set_defaults(fn=cmd_generate)
    sp = sub.add_parser("simulate", help="simulate L, ESR, Q and current sharing"); common(sp)
    sp.add_argument("--freq", type=float, default=6.78e6); sp.add_argument("--temp", type=float, default=25.0)
    sp.add_argument("--quality", choices=["fast", "standard", "fine"], default="standard")
    sp.add_argument("--sweep", action="store_true", help="also solve at +-3 %% and +-6 %% of --freq")
    sp.add_argument("--report"); sp.add_argument("--json"); sp.add_argument("-q", "--quiet", action="store_true")
    sp.set_defaults(fn=cmd_simulate)
    sp = sub.add_parser("sweep", help="vary one parameter"); common(sp)
    sp.add_argument("--vary", required=True, metavar="KEY=V1,V2,...")
    sp.add_argument("--quality", choices=["fast", "standard", "fine"], default="fast")
    sp.add_argument("--freq", type=float, default=6.78e6); sp.add_argument("--no-sim", action="store_true")
    sp.set_defaults(fn=cmd_sweep)
    sp = sub.add_parser("export", help="export to other tools"); common(sp)
    sp.add_argument("kind", choices=["fasthenry", "openems", "kicad"]); sp.add_argument("-o", "--out", required=True)
    sp.add_argument("--freq", type=float, default=6.78e6); sp.add_argument("--nwinc", type=int, default=7)
    sp.add_argument("--nhinc", type=int, default=3); sp.add_argument("--name", default="litz_coil")
    sp.set_defaults(fn=cmd_export)
    args = ap.parse_args(argv)
    try:
        args.fn(args)
    except ValueError as e:
        raise SystemExit("error: %s" % e)


if __name__ == "__main__":
    main()
