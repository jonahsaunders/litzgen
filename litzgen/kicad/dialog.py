"""wxPython dialog for the LitzGen action plugin."""
from __future__ import annotations

import json
import os
import threading
import traceback
import webbrowser

import wx

import pcbnew

from ..core.drc import check
from ..core.export.fasthenry import fasthenry_text
from ..core.export.openems import openems_script
from ..core.geometry import generate
from ..core.params import CoilParams, SimParams
from . import pcbnew_writer as W

# (key path, label, kind, tooltip)
FIELDS = [
    ("Coil", None, None, None),
    ("d_out", "Outer diameter (mm)", float, "Outer edge of the outer turn"),
    ("n_turns", "Turns", float, ""),
    ("d_in", "Inner diameter (mm), blank = use gap", "optfloat", "Inner edge of the inner turn"),
    ("turn_gap", "Turn gap (mm), blank = use d_in", "optfloat", "Edge-to-edge gap between turns"),
    ("clockwise", "Clockwise (top view)", bool, ""),
    ("center_x", "Centre X (mm)", float, ""),
    ("center_y", "Centre Y (mm)", float, ""),
    ("start_angle_deg", "Outer terminal angle (deg)", float, ""),
    ("Strands", None, None, None),
    ("n_cols", "Strands per layer", int, "Columns of the strand matrix"),
    ("strand_width", "Strand width (mm)", float, ""),
    ("strand_gap", "Strand gap (mm)", float, ""),
    ("stackup.n_layers", "Copper layers used", int, ""),
    ("stackup.kicad_layers", "KiCad layers (comma list)", "strlist", "Top to bottom"),
    ("stackup.copper_thickness", "Copper thickness (mm)", float, ""),
    ("stackup.dielectric", "Dielectric thicknesses (mm, comma list)", "floatlist", "Between consecutive coil layers"),
    ("omit_rings", "Leave rings empty (comma list)", "intlist", "e.g. 1 = hollow bundle"),
    ("Transposition", None, None, None),
    ("steps_per_turn", "Transposition steps per turn", int, "0 = untransposed"),
    ("stagger_alternate_turns", "Stagger zones on odd turns", bool, ""),
    ("jog_width", "Jog width (mm), blank = strand width", "optfloat", "Neck traces inside jogs"),
    ("clearance", "Clearance (mm)", float, ""),
    ("hole_clearance", "Hole-to-hole clearance (mm)", float, ""),
    ("via.drill", "Via drill (mm)", float, ""),
    ("via.pad", "Via pad (mm)", float, ""),
    ("via.plating", "Via plating (mm)", float, ""),
    ("terminal_margin_deg", "Zone-free margin at terminals (deg)", float, ""),
    ("terminal_length", "Terminal bar length (mm)", float, ""),
    ("terminal_drill", "Terminal hole drill (mm)", float, ""),
    ("KiCad", None, None, None),
    ("net_prefix", "Net name", str, ""),
    ("per_strand_nets", "One net per strand (DRC board, no terminals)", bool, ""),
    ("group_name", "Group name", str, "Used to find and replace the coil"),
]


def _get(d, path):
    for k in path.split("."):
        d = d[k]
    return d


def _set(d, path, v):
    ks = path.split(".")
    for k in ks[:-1]:
        d = d[k]
    d[ks[-1]] = v


class LitzDialog(wx.Dialog):
    def __init__(self, parent, board):
        super().__init__(parent, title="LitzGen - multi-bundle PCB Litz coil", size=(760, 820),
                         style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.board = board
        self.result = None
        d = None
        for name in W.list_coils(board):
            d = W.read_params_from_board(board, name)
            break
        self.params = CoilParams.from_dict(d) if d else CoilParams()
        self._build()
        self._load(self.params)

    # ----------------------------------------------------------------- UI
    def _build(self):
        outer = wx.BoxSizer(wx.VERTICAL)
        split = wx.BoxSizer(wx.HORIZONTAL)
        sw = wx.ScrolledWindow(self, style=wx.VSCROLL)
        sw.SetScrollRate(0, 12)
        grid = wx.FlexGridSizer(0, 2, 4, 8)
        grid.AddGrowableCol(1)
        self.ctrls = {}
        for key, label, kind, tip in FIELDS:
            if label is None:
                t = wx.StaticText(sw, label=key)
                f = t.GetFont(); f.SetWeight(wx.FONTWEIGHT_BOLD); t.SetFont(f)
                grid.Add(t, 0, wx.TOP, 8); grid.Add((0, 0))
                continue
            grid.Add(wx.StaticText(sw, label=label), 0, wx.ALIGN_CENTER_VERTICAL)
            c = wx.CheckBox(sw) if kind is bool else wx.TextCtrl(sw)
            if tip:
                c.SetToolTip(tip)
            self.ctrls[key] = (c, kind)
            grid.Add(c, 1, wx.EXPAND)
        sim_box = wx.StaticText(sw, label="Simulation")
        f = sim_box.GetFont(); f.SetWeight(wx.FONTWEIGHT_BOLD); sim_box.SetFont(f)
        grid.Add(sim_box, 0, wx.TOP, 8); grid.Add((0, 0))
        grid.Add(wx.StaticText(sw, label="Frequency (MHz)"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.freq = wx.TextCtrl(sw, value="6.78"); grid.Add(self.freq, 1, wx.EXPAND)
        grid.Add(wx.StaticText(sw, label="Quality"), 0, wx.ALIGN_CENTER_VERTICAL)
        self.quality = wx.Choice(sw, choices=["fast", "standard", "fine"]); self.quality.SetSelection(1)
        grid.Add(self.quality, 1, wx.EXPAND)
        self.do_sweep = wx.CheckBox(sw, label="Frequency sweep (+-6 %)")
        grid.Add((0, 0)); grid.Add(self.do_sweep)
        sw.SetSizer(grid)
        split.Add(sw, 1, wx.EXPAND | wx.ALL, 6)
        outer.Add(split, 1, wx.EXPAND)
        self.log = wx.TextCtrl(self, style=wx.TE_MULTILINE | wx.TE_READONLY, size=(-1, 170))
        outer.Add(self.log, 0, wx.EXPAND | wx.LEFT | wx.RIGHT, 6)
        self.gauge = wx.Gauge(self, range=100)
        outer.Add(self.gauge, 0, wx.EXPAND | wx.ALL, 6)
        btns = wx.WrapSizer()
        for label, fn in [("Paper defaults", self.on_defaults), ("Load JSON...", self.on_load), ("Save JSON...", self.on_save),
                          ("Check DRC", self.on_check), ("Generate", self.on_generate), ("Simulate", self.on_simulate),
                          ("Export FastHenry...", self.on_fasthenry), ("Export openEMS...", self.on_openems),
                          ("Remove coil", self.on_remove), ("Close", lambda e: self.EndModal(wx.ID_CLOSE))]:
            b = wx.Button(self, label=label)
            b.Bind(wx.EVT_BUTTON, fn)
            btns.Add(b, 0, wx.ALL, 3)
        outer.Add(btns, 0, wx.ALL, 3)
        self.SetSizer(outer)

    def _load(self, p: CoilParams):
        d = p.to_dict()
        for key, (c, kind) in self.ctrls.items():
            v = _get(d, key)
            if kind is bool:
                c.SetValue(bool(v))
            elif kind in ("strlist", "floatlist", "intlist"):
                c.SetValue(", ".join(str(x) for x in v))
            else:
                c.SetValue("" if v is None else str(v))

    def _read(self) -> CoilParams:
        d = CoilParams().to_dict()
        for key, (c, kind) in self.ctrls.items():
            if kind is bool:
                v = c.GetValue()
            else:
                s = c.GetValue().strip()
                if kind == "optfloat":
                    v = float(s) if s else None
                elif kind == "strlist":
                    v = [x.strip() for x in s.split(",") if x.strip()]
                elif kind == "floatlist":
                    v = [float(x) for x in s.split(",") if x.strip()]
                elif kind == "intlist":
                    v = [int(x) for x in s.split(",") if x.strip()]
                else:
                    v = kind(s)
            _set(d, key, v)
        return CoilParams.from_dict(d)

    def say(self, msg):
        self.log.AppendText(msg + "\n")

    def _safe(self, fn):
        try:
            return fn()
        except Exception as e:  # show every failure in the log instead of crashing KiCad
            self.say("ERROR: %s" % e)
            self.say(traceback.format_exc(limit=2))
            return None

    def _geo(self):
        p = self._read()
        g = generate(p)
        for w in g.warnings:
            self.say("warning: " + w)
        return g

    # ----------------------------------------------------------------- actions
    def on_defaults(self, e):
        self._load(CoilParams())
        self.say("Loaded Kale & Wicht (WPTCE 2026) Table I defaults.")

    def on_load(self, e):
        with wx.FileDialog(self, "Load parameters", wildcard="JSON (*.json)|*.json", style=wx.FD_OPEN) as fd:
            if fd.ShowModal() == wx.ID_OK:
                self._safe(lambda: self._load(CoilParams.from_json(open(fd.GetPath()).read())))

    def on_save(self, e):
        with wx.FileDialog(self, "Save parameters", wildcard="JSON (*.json)|*.json",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fd:
            if fd.ShowModal() == wx.ID_OK:
                self._safe(lambda: open(fd.GetPath(), "w").write(self._read().to_json()))

    def on_check(self, e):
        def run():
            g = self._geo()
            s = g.summary()
            self.say("%d strands, %d zones, %d vias %s" % (s["strands"], s["zones"], s["vias_total"], json.dumps(s["vias_by_type"])))
            v = check(g)
            self.say("DRC (strand-to-strand, hole-to-hole): %d violation(s)" % len(v))
            for x in v[:10]:
                self.say("  " + str(x))
            return g
        self._safe(run)

    def on_generate(self, e):
        def run():
            g = self._geo()
            v = check(g)
            if v:
                if wx.MessageBox("%d clearance violations. Place anyway?" % len(v), "LitzGen",
                                 wx.YES_NO | wx.ICON_WARNING) != wx.YES:
                    return
            n = W.place_coil(self.board, g, replace=True)
            self.say("Placed %s: %s" % (g.params.group_name, n))
        self._safe(run)

    def on_remove(self, e):
        self._safe(lambda: self.say("Removed %d items." % W.remove_coil(self.board, self._read().group_name)))

    def _board_dir(self):
        fn = self.board.GetFileName()
        return os.path.dirname(fn) if fn else os.path.expanduser("~")

    def on_simulate(self, e):
        g = self._safe(self._geo)
        if g is None:
            return
        sim = SimParams.preset(self.quality.GetStringSelection(), frequency=float(self.freq.GetValue()) * 1e6)
        sweep = [sim.frequency * k for k in (0.94, 0.97, 1.03, 1.06)] if self.do_sweep.GetValue() else None
        self.say("Simulating (%s quality) ..." % self.quality.GetStringSelection())

        def prog(stage, frac):
            wx.CallAfter(self.gauge.SetValue, int(frac * 100))

        def work():
            try:
                from ..core.report import html_report
                from ..core.solver.solve import simulate
                res = simulate(g, sim, sweep_freqs=sweep, progress=prog)
                path = os.path.join(self._board_dir(), "%s_report.html" % g.params.group_name)
                with open(path, "w") as f:
                    f.write(html_report(g, res, len(check(g))))
                wx.CallAfter(self._sim_done, res, path)
            except Exception as ex:
                wx.CallAfter(self.say, "ERROR: %s\n%s" % (ex, traceback.format_exc(limit=3)))
        threading.Thread(target=work, daemon=True).start()

    def _sim_done(self, res, path):
        self.gauge.SetValue(100)
        self.say("L = %.4f uH, ESR = %.1f mOhm, Q = %.1f, Rac/Rdc = %.2f" % (res.L * 1e6, res.R * 1e3, res.Q, res.R / res.R_dc))
        self.say("Ring current share: " + ", ".join("ring %s %.1f %%" % (k, abs(v) * 100) for k, v in sorted(res.ring_share().items())))
        self.say("Report: " + path)
        webbrowser.open("file://" + path)

    def _export(self, title, wildcard, fn):
        with wx.FileDialog(self, title, defaultDir=self._board_dir(), wildcard=wildcard,
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fd:
            if fd.ShowModal() == wx.ID_OK:
                def run():
                    with open(fd.GetPath(), "w") as f:
                        f.write(fn(self._geo()))
                    self.say("Wrote " + fd.GetPath())
                self._safe(run)

    def on_fasthenry(self, e):
        self._export("Export FastHenry", "FastHenry (*.inp)|*.inp",
                     lambda g: fasthenry_text(g, SimParams(frequency=float(self.freq.GetValue()) * 1e6)))

    def on_openems(self, e):
        self._export("Export openEMS script", "Python (*.py)|*.py", lambda g: openems_script(g))
