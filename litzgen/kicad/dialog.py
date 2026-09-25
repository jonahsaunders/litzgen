"""wxPython dialog for the LitzGen action plugin.

Layout follows Apple's Human Interface Guidelines for settings windows: parameters
are grouped into tabs, labels are short and end in a colon, units sit after the
field, the primary action is the trailing default button, destructive actions are
separated and confirmed, and results are reported in plain language in a status
line (the full activity log is one click away).
"""
from __future__ import annotations

import json
import os
import re
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

# Tabs of (key path, label, kind, unit, tooltip, placeholder); a 1-tuple is a section heading.
PAGES = [
    ("Coil", [
        ("Size",),
        ("d_out", "Outer diameter", float, "mm", "Outer edge of the outer turn.", None),
        ("n_turns", "Turns", float, "", "Use whole or half turns so both terminals land at defined angles.", None),
        ("d_in", "Inner diameter", "optfloat", "mm",
         "Inner edge of the inner turn. Leave blank to derive it from the turn gap.", "Auto"),
        ("turn_gap", "Turn gap", "optfloat", "mm",
         "Edge-to-edge gap between turns. Leave blank to derive it from the inner diameter.", "Auto"),
        ("clockwise", "Wind clockwise (seen from the top)", bool, "", "Spiral direction as seen on F.Cu.", None),
        ("Placement",),
        ("center_x", "Centre X", float, "mm", "Board position of the coil centre.", None),
        ("center_y", "Centre Y", float, "mm", "Board position of the coil centre.", None),
        ("start_angle_deg", "Outer terminal angle", float, "°", "Angle of the outer terminal around the centre.", None),
    ]),
    ("Strands & Layers", [
        ("Strand matrix",),
        ("n_cols", "Strands per layer", int, "", "Columns of the strand matrix (at least 2).", None),
        ("strand_width", "Strand width", float, "mm", "Track width of each strand.", None),
        ("strand_gap", "Strand gap", float, "mm", "Edge-to-edge gap between neighbouring strands.", None),
        ("omit_rings", "Empty rings", "intlist", "",
         "Comma-separated ring numbers to leave empty, e.g. 1 for a hollow bundle.", "None"),
        ("Stack-up",),
        ("stackup.n_layers", "Copper layers", int, "", "Copper layers used by the coil (at least 2).", None),
        ("stackup.kicad_layers", "KiCad layers", "strlist", "",
         "One KiCad layer per copper layer, top to bottom, comma-separated. Missing inner layers are added.", None),
        ("stackup.copper_thickness", "Copper thickness", float, "mm", "", None),
        ("stackup.dielectric", "Dielectric thicknesses", "floatlist", "mm",
         "Between consecutive coil layers, comma-separated: one fewer than the copper layers.", None),
    ]),
    ("Transposition & Vias", [
        ("Transposition",),
        ("steps_per_turn", "Zones per turn", int, "",
         "Transposition zones per turn. 0 = untransposed, -1 = the densest that fits.", None),
        ("stagger_alternate_turns", "Stagger zones on odd turns", bool, "", "", None),
        ("jog_width", "Jog width", "optfloat", "mm", "Trace width inside the transposition jogs.", "Strand width"),
        ("clearance", "Clearance", float, "mm", "Copper-to-copper clearance used by the generator and the check.", None),
        ("hole_clearance", "Hole clearance", float, "mm", "Hole edge to hole edge.", None),
        ("Vias",),
        ("via.drill", "Drill", float, "mm", "", None),
        ("via.pad", "Pad diameter", float, "mm", "", None),
        ("via.plating", "Barrel plating", float, "mm", "Plating thickness, used for via resistance.", None),
    ]),
    ("Terminals & Board", [
        ("Terminals",),
        ("terminal_margin_deg", "Zone-free margin", float, "°", "No transposition zones this close to a terminal.", None),
        ("terminal_length", "Bar length", float, "mm", "Arc length of each terminal bar.", None),
        ("terminal_drill", "Hole drill", float, "mm", "Stitching and wire holes in the terminal bars.", None),
        ("Board",),
        ("net_prefix", "Net name", str, "", "", None),
        ("group_name", "Group name", str, "", "LitzGen finds, updates and removes the coil by this group name.", None),
        ("per_strand_nets", "One net per strand (no terminals)", bool, "",
         "For a throwaway DRC board: lets KiCad's own DRC confirm there are no shorts between strands.", None),
    ]),
]
FIELDS = [f for _, rows in PAGES for f in rows if len(f) > 1]
LABELS = {f[0].split(".")[-1]: f[1] for f in FIELDS}
QUALITIES = [("fast", "Fast"), ("standard", "Standard"), ("fine", "Fine")]


def _get(d, path):
    for k in path.split("."):
        d = d[k]
    return d


def _set(d, path, v):
    ks = path.split(".")
    for k in ks[:-1]:
        d = d[k]
    d[ks[-1]] = v


def humanize(msg: str) -> str:
    """Replace parameter identifiers in core error messages with the labels shown in the dialog."""
    for key in sorted(LABELS, key=len, reverse=True):
        msg = re.sub(r"\b%s\b" % re.escape(key), LABELS[key], msg)
    msg = msg.replace(">=", "≥").replace("<=", "≤").strip()
    if msg and msg[-1] not in ".!?":
        msg += "."
    return msg[:1].upper() + msg[1:]


class _FieldError(ValueError):
    def __init__(self, key, message):
        super().__init__(message)
        self.key = key


class _Cancelled(Exception):
    pass


class LitzDialog(wx.Dialog):
    def __init__(self, parent, board):
        super().__init__(parent, title="LitzGen", style=wx.DEFAULT_DIALOG_STYLE | wx.RESIZE_BORDER)
        self.board = board
        self.result = None
        self._busy = False
        self._cancel = None
        self._report = None
        d = None
        for name in W.list_coils(board):
            d = W.read_params_from_board(board, name)
            break
        self.params = CoilParams.from_dict(d) if d else CoilParams()
        self._build()
        self._load(self.params)
        self._refresh_state()
        self.Bind(wx.EVT_CLOSE, self.on_close)

    # ----------------------------------------------------------------- UI
    def _heading(self, parent, text):
        t = wx.StaticText(parent, label=text)
        t.SetFont(t.GetFont().Bold())
        return t

    def _build(self):
        dip = self.FromDIP
        outer = wx.BoxSizer(wx.VERTICAL)

        self.title = wx.StaticText(self, label="Multi-bundle PCB Litz coil")
        self.title.SetFont(self.title.GetFont().Bold().Scaled(1.25))
        self.subtitle = wx.StaticText(self)
        self.subtitle.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT))
        outer.Add(self.title, 0, wx.LEFT | wx.RIGHT | wx.TOP, dip(16))
        outer.AddSpacer(dip(2))
        outer.Add(self.subtitle, 0, wx.LEFT | wx.RIGHT, dip(16))
        outer.AddSpacer(dip(12))

        self.book = wx.Notebook(self)
        self.ctrls = {}
        self.labels = {}
        self.tips = {}
        self.page_of = {}
        for title, rows in PAGES:
            page = wx.Panel(self.book)
            self._form(page, rows)
            self.book.AddPage(page, title.replace("&", "&&"))
        self.book.AddPage(self._sim_page(), "Simulation")
        outer.Add(self.book, 1, wx.EXPAND | wx.LEFT | wx.RIGHT, dip(16))

        # status line: icon + plain-language result, progress while busy
        status = wx.BoxSizer(wx.HORIZONTAL)
        self.status_icon = wx.StaticBitmap(self, bitmap=self._icon(wx.ART_INFORMATION))
        self.status = wx.StaticText(self, label="")
        self.open_report = wx.Button(self, label="Open Report")
        self.open_report.Bind(wx.EVT_BUTTON, lambda e: self._report and webbrowser.open("file://" + self._report))
        self.open_report.Hide()
        status.Add(self.status_icon, 0, wx.ALIGN_TOP | wx.TOP, dip(1))
        status.Add(self.status, 1, wx.LEFT, dip(8))
        status.Add(self.open_report, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, dip(8))
        outer.Add(status, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, dip(16))

        self.progress_row = wx.BoxSizer(wx.HORIZONTAL)
        self.gauge = wx.Gauge(self, range=100, size=(-1, dip(8)))
        self.pulse = wx.Timer(self)
        self.Bind(wx.EVT_TIMER, lambda e: self.gauge.Pulse(), self.pulse)
        self.cancel_btn = wx.Button(self, label="Cancel")
        self.cancel_btn.Bind(wx.EVT_BUTTON, self.on_cancel)
        self.progress_row.Add(self.gauge, 1, wx.ALIGN_CENTER_VERTICAL)
        self.progress_row.Add(self.cancel_btn, 0, wx.LEFT, dip(8))
        outer.Add(self.progress_row, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, dip(16))

        self.details = wx.CollapsiblePane(self, label="Activity log")
        self.log = wx.TextCtrl(self.details.GetPane(), style=wx.TE_MULTILINE | wx.TE_READONLY, size=(-1, dip(140)))
        self.log.SetName("Activity log")
        ps = wx.BoxSizer(wx.VERTICAL)
        ps.Add(self.log, 1, wx.EXPAND | wx.TOP, dip(4))
        self.details.GetPane().SetSizer(ps)
        self.details.Bind(wx.EVT_COLLAPSIBLEPANE_CHANGED, self._on_details)
        outer.Add(self.details, 0, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, dip(12))

        # button bar: destructive + menu on the leading edge, actions and the default on the trailing edge
        bar = wx.BoxSizer(wx.HORIZONTAL)
        self.remove_btn = wx.Button(self, label="Remove Coil")
        self.remove_btn.Bind(wx.EVT_BUTTON, self.on_remove)
        self.more_btn = wx.Button(self, label="More ▾")
        self.more_btn.Bind(wx.EVT_BUTTON, self.on_more)
        self.check_btn = wx.Button(self, label="Check Clearances")
        self.check_btn.Bind(wx.EVT_BUTTON, self.on_check)
        self.sim_btn = wx.Button(self, label="Simulate")
        self.sim_btn.Bind(wx.EVT_BUTTON, self.on_simulate)
        self.close_btn = wx.Button(self, wx.ID_CLOSE, "Close")
        self.close_btn.Bind(wx.EVT_BUTTON, lambda e: self.Close())
        self.place_btn = wx.Button(self, wx.ID_OK, "Place Coil")
        self.place_btn.Bind(wx.EVT_BUTTON, self.on_generate)
        self.place_btn.SetDefault()
        self.SetEscapeId(wx.ID_CLOSE)
        bar.Add(self.remove_btn)
        bar.Add(self.more_btn, 0, wx.LEFT, dip(8))
        bar.AddStretchSpacer()
        bar.Add(self.check_btn)
        bar.Add(self.sim_btn, 0, wx.LEFT, dip(8))
        bar.AddSpacer(dip(24))
        bar.Add(self.close_btn)
        bar.Add(self.place_btn, 0, wx.LEFT, dip(8))
        outer.Add(bar, 0, wx.EXPAND | wx.ALL, dip(16))

        self.SetSizer(outer)
        self._set_status("info", "Set the coil parameters, then place the coil on the board.")
        self.Fit()  # sized with the progress row showing, so it never squeezes the tabs
        self.SetMinSize(self.GetSize())
        self._show_progress(False)
        self.Bind(wx.EVT_SIZE, lambda e: (e.Skip(), wx.CallAfter(self._ui, self._wrap_status)))
        self.ctrls[FIELDS[0][0]][0].SetFocus()

    def _form(self, page, rows):
        dip = self.FromDIP
        grid = wx.FlexGridSizer(0, 2, dip(8), dip(8))
        grid.AddGrowableCol(0, 2)
        grid.AddGrowableCol(1, 3)
        first = True
        for row in rows:
            if len(row) == 1:
                if not first:
                    grid.Add((0, dip(6))); grid.Add((0, 0))
                grid.Add(self._heading(page, row[0]), 0, wx.ALIGN_RIGHT)
                grid.Add((0, 0))
                first = False
                continue
            key, label, kind, unit, tip, hint = row
            if kind is bool:
                c = wx.CheckBox(page, label=label)
                grid.Add((0, 0)); grid.Add(c, 0, wx.ALIGN_CENTER_VERTICAL)
            else:
                lab = wx.StaticText(page, label=label + ":")
                wide = kind in ("strlist", "floatlist", "intlist", str)
                c = wx.TextCtrl(page, size=(dip(220 if wide else 110), -1))
                if hint:
                    c.SetHint(hint)
                c.Bind(wx.EVT_TEXT, lambda e, k=key: self._validate_field(k))
                self.labels[key] = lab
                grid.Add(lab, 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
                grid.Add(self._with_unit(page, c, unit), 0, wx.ALIGN_CENTER_VERTICAL)
            c.SetName(label)
            self.tips[key] = tip or None
            c.SetToolTip(self.tips[key])
            self.ctrls[key] = (c, kind)
            self.page_of[key] = self.book.GetPageCount()
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(grid, 1, wx.EXPAND | wx.ALL, dip(16))
        page.SetSizer(s)

    def _with_unit(self, page, ctrl, unit):
        row = wx.BoxSizer(wx.HORIZONTAL)
        row.Add(ctrl, 0, wx.ALIGN_CENTER_VERTICAL)
        if unit:
            row.Add(wx.StaticText(page, label=unit), 0, wx.ALIGN_CENTER_VERTICAL | wx.LEFT, self.FromDIP(6))
        return row

    def _sim_page(self):
        dip = self.FromDIP
        page = wx.Panel(self.book)
        grid = wx.FlexGridSizer(0, 2, dip(8), dip(8))
        grid.AddGrowableCol(0, 2)
        grid.AddGrowableCol(1, 3)
        grid.Add(self._heading(page, "Solver"), 0, wx.ALIGN_RIGHT); grid.Add((0, 0))
        self.labels["freq"] = wx.StaticText(page, label="Frequency:")
        grid.Add(self.labels["freq"], 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.freq = wx.TextCtrl(page, value="6.78", size=(dip(110), -1), name="Frequency")
        self.freq.Bind(wx.EVT_TEXT, lambda e: self._validate_field("freq"))
        grid.Add(self._with_unit(page, self.freq, "MHz"), 0, wx.ALIGN_CENTER_VERTICAL)
        grid.Add(wx.StaticText(page, label="Quality:"), 0, wx.ALIGN_RIGHT | wx.ALIGN_CENTER_VERTICAL)
        self.quality = wx.Choice(page, choices=[t for _, t in QUALITIES], name="Quality")
        self.quality.SetSelection(1)
        self.quality.SetToolTip("Fast for exploring, Fine for final numbers (slower).")
        grid.Add(self.quality, 0)
        grid.Add((0, 0))
        self.do_sweep = wx.CheckBox(page, label="Also sweep ±6 % around the frequency")
        grid.Add(self.do_sweep)
        grid.Add((0, 0))
        note = wx.StaticText(page, label="Simulate writes an HTML report next to the board file.")
        note.SetForegroundColour(wx.SystemSettings.GetColour(wx.SYS_COLOUR_GRAYTEXT))
        grid.Add(note)
        s = wx.BoxSizer(wx.VERTICAL)
        s.Add(grid, 1, wx.EXPAND | wx.ALL, dip(16))
        page.SetSizer(s)
        self.tips["freq"] = None
        self.page_of["freq"] = self.book.GetPageCount()
        return page

    def _icon(self, art):
        return wx.ArtProvider.GetBitmapBundle(art, wx.ART_MESSAGE_BOX, wx.Size(16, 16))

    def _on_details(self, e):
        self.Layout()
        self.SetSize((self.GetSize().width, max(self.GetBestSize().height, self.GetMinSize().height)))

    def _set_status(self, level, text):
        art = {"info": wx.ART_INFORMATION, "ok": wx.ART_INFORMATION, "warn": wx.ART_WARNING, "error": wx.ART_ERROR}[level]
        self.status_icon.SetBitmap(self._icon(art))
        self.status.SetLabel(text)
        self._wrap_status()

    def _wrap_status(self):
        room = self.book.GetSize().width - self.FromDIP(24)
        if self.open_report.IsShown():
            room -= self.open_report.GetSize().width + self.FromDIP(8)
        self.status.Wrap(max(room, self.FromDIP(200)))
        self.Layout()

    def _show_progress(self, on):
        if not on:
            self.pulse.Stop()
        self.progress_row.ShowItems(on)
        self.Layout()

    def _refresh_state(self):
        name = self.ctrls["group_name"][0].GetValue().strip() or "LitzCoil"
        on_board = W.find_group(self.board, name) is not None
        self.subtitle.SetLabel(("Editing “%s” on this board." % name) if on_board
                               else "New coil “%s”. Not yet on this board." % name)
        self.place_btn.SetLabel("Update Coil" if on_board else "Place Coil")
        self.remove_btn.Enable(on_board and not self._busy)
        self.Layout()

    def _set_busy(self, busy):
        self._busy = busy
        for b in (self.check_btn, self.sim_btn, self.place_btn, self.more_btn):
            b.Enable(not busy)
        self._refresh_state()
        self._show_progress(busy)

    # ----------------------------------------------------------------- values
    def _load(self, p: CoilParams):
        d = p.to_dict()
        for key, (c, kind) in self.ctrls.items():
            v = _get(d, key)
            if kind is bool:
                c.SetValue(bool(v))
            elif kind in ("strlist", "floatlist", "intlist"):
                c.ChangeValue(", ".join(str(x) for x in v))
            else:
                c.ChangeValue("" if v is None else str(v))
        for key in list(self.ctrls) + ["freq"]:
            self._validate_field(key)
        self._refresh_state()

    def _parse(self, key):
        if key == "freq":
            s = self.freq.GetValue().strip()
            try:
                v = float(s)
            except ValueError:
                raise _FieldError(key, "Frequency must be a number in MHz, e.g. 6.78.")
            if v <= 0:
                raise _FieldError(key, "Frequency must be greater than zero.")
            return v
        c, kind = self.ctrls[key]
        label = next(f[1] for f in FIELDS if f[0] == key)
        if kind is bool:
            return c.GetValue()
        s = c.GetValue().strip()
        try:
            if kind == "optfloat":
                return float(s) if s else None
            if kind == "strlist":
                return [x.strip() for x in s.split(",") if x.strip()]
            if kind == "floatlist":
                return [float(x) for x in s.split(",") if x.strip()]
            if kind == "intlist":
                return [int(x) for x in s.split(",") if x.strip()]
            if kind is str:
                if not s:
                    raise ValueError
                return s
            return kind(s)
        except ValueError:
            what = {int: "a whole number", float: "a number", "optfloat": "a number or blank",
                    "floatlist": "numbers separated by commas", "intlist": "whole numbers separated by commas",
                    str: "filled in"}.get(kind, "valid")
            raise _FieldError(key, "%s must be %s." % (label, what))

    def _validate_field(self, key):
        c = self.freq if key == "freq" else self.ctrls[key][0]
        lab = self.labels.get(key)
        try:
            self._parse(key)
            ok, tip = True, None
        except _FieldError as e:
            ok, tip = False, str(e)
        if lab is not None:
            lab.SetForegroundColour(wx.NullColour if ok else wx.Colour(0xD7, 0x3A, 0x2F))
            lab.Refresh()
        if key == "group_name":
            self._refresh_state()
        c.SetToolTip(tip or self.tips.get(key))
        return ok

    def _read(self) -> CoilParams:
        d = CoilParams().to_dict()
        for key in self.ctrls:
            _set(d, key, self._parse(key))
        return CoilParams.from_dict(d)

    def _focus(self, key):
        self.book.SetSelection(self.page_of[key])
        c = self.freq if key == "freq" else self.ctrls[key][0]
        c.SetFocus()
        if isinstance(c, wx.TextCtrl):
            c.SelectAll()

    # ----------------------------------------------------------------- feedback
    def say(self, msg):
        self.log.AppendText(msg + "\n")

    def _safe(self, fn):
        try:
            return fn()
        except _FieldError as e:
            self._focus(e.key)
            self._set_status("error", str(e))
        except Exception as e:  # never crash KiCad: explain in the status line, keep details in the log
            self._set_status("error", humanize(str(e)) or type(e).__name__)
            self.say("Error: %s" % e)
            self.say(traceback.format_exc(limit=3))
        return None

    def _geo(self):
        g = generate(self._read())
        for w in g.warnings:
            self.say("Warning: " + humanize(w))
        return g

    def _confirm(self, message, detail, action):
        dlg = wx.MessageDialog(self, message, "LitzGen", wx.YES_NO | wx.NO_DEFAULT | wx.ICON_WARNING)
        dlg.SetExtendedMessage(detail)
        dlg.SetYesNoLabels(action, "Cancel")
        ok = dlg.ShowModal() == wx.ID_YES
        dlg.Destroy()
        return ok

    def _board_dir(self):
        fn = self.board.GetFileName()
        return os.path.dirname(fn) if fn else os.path.expanduser("~")

    # ----------------------------------------------------------------- actions
    def on_more(self, e):
        menu = wx.Menu()
        items = [("Load Parameters…", self.on_load), ("Save Parameters…", self.on_save), None,
                 ("Export FastHenry…", self.on_fasthenry), ("Export openEMS Script…", self.on_openems), None,
                 ("Reset to Paper Defaults", self.on_defaults)]
        for it in items:
            if it is None:
                menu.AppendSeparator()
            else:
                mi = menu.Append(wx.ID_ANY, it[0])
                self.Bind(wx.EVT_MENU, it[1], mi)
        self.PopupMenu(menu, self.more_btn.GetPosition() + (0, self.more_btn.GetSize().height))
        menu.Destroy()

    def on_defaults(self, e):
        if self._confirm("Reset all parameters to the paper's defaults?",
                         "Every field is replaced with the Kale & Wicht (WPTCE 2026) Table I values. "
                         "The coil on the board is not changed until you place it.", "Reset"):
            self._load(CoilParams())
            self._set_status("info", "Parameters reset to the paper's defaults.")

    def on_load(self, e):
        with wx.FileDialog(self, "Load Parameters", wildcard="LitzGen parameters (*.json)|*.json",
                           style=wx.FD_OPEN | wx.FD_FILE_MUST_EXIST) as fd:
            if fd.ShowModal() == wx.ID_OK:
                def run():
                    with open(fd.GetPath(), encoding="utf-8") as f:
                        self._load(CoilParams.from_json(f.read()))
                    self._set_status("info", "Loaded %s." % os.path.basename(fd.GetPath()))
                self._safe(run)

    def on_save(self, e):
        with wx.FileDialog(self, "Save Parameters", defaultDir=self._board_dir(),
                           defaultFile="%s.json" % self.ctrls["group_name"][0].GetValue().strip(),
                           wildcard="LitzGen parameters (*.json)|*.json",
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fd:
            if fd.ShowModal() == wx.ID_OK:
                def run():
                    text = self._read().to_json()
                    with open(fd.GetPath(), "w", encoding="utf-8") as f:
                        f.write(text)
                    self._set_status("ok", "Saved %s." % os.path.basename(fd.GetPath()))
                self._safe(run)

    def on_check(self, e):
        def run():
            g = self._geo()
            s = g.summary()
            v = check(g)
            vias = ", ".join("%d %s" % (n, k) for k, n in s["vias_by_type"].items())
            self.say("%d strands, %d zones, %d vias (%s)" % (s["strands"], s["zones"], s["vias_total"], vias))
            self.say("Clearance check (strand-to-strand, hole-to-hole): %d violation(s)" % len(v))
            for x in v[:10]:
                self.say("  " + str(x))
            if v:
                self._set_status("warn", "%d clearance violation%s. The first ones are listed in the activity log."
                                 % (len(v), "" if len(v) == 1 else "s"))
                self.details.Expand(); self._on_details(None)
            else:
                self._set_status("ok", "No clearance violations. %d strands, %d transposition zones, %d vias."
                                 % (s["strands"], s["zones"], s["vias_total"]))
        self._safe(run)

    def on_generate(self, e):
        def run():
            g = self._geo()
            v = check(g)
            if v and not self._confirm(
                    "The coil has %d clearance violation%s." % (len(v), "" if len(v) == 1 else "s"),
                    "Placing it anyway adds copper that will fail DRC. "
                    "Choose Check Clearances to list them.", "Place Anyway"):
                return
            replacing = W.find_group(self.board, g.params.group_name) is not None
            n = W.place_coil(self.board, g, replace=True)
            counts = ", ".join("%d %s" % (c, k) for k, c in n.items() if c)
            self.say("Placed %s: %s" % (g.params.group_name, counts))
            self._set_status("ok", "%s “%s” on the board (%s)." % ("Updated" if replacing else "Placed",
                                                                 g.params.group_name, counts))
            self._refresh_state()
        self._safe(run)

    def on_remove(self, e):
        name = self.ctrls["group_name"][0].GetValue().strip()
        if not self._confirm("Remove “%s” from the board?" % name,
                             "Every track, via and terminal in the group is deleted.", "Remove"):
            return

        def run():
            n = W.remove_coil(self.board, name)
            self.say("Removed %d items." % n)
            self._set_status("info", "Removed “%s” (%d items)." % (name, n))
            self._refresh_state()
        self._safe(run)

    def on_simulate(self, e):
        def prepare():
            g = self._geo()
            return g, self._parse("freq")
        r = self._safe(prepare)
        if r is None:
            return
        g, mhz = r
        quality = QUALITIES[self.quality.GetSelection()][0]
        sim = SimParams.preset(quality, frequency=mhz * 1e6)
        sweep = [sim.frequency * k for k in (0.94, 0.97, 1.03, 1.06)] if self.do_sweep.GetValue() else None
        path = os.path.join(self._board_dir(), "%s_report.html" % g.params.group_name)
        cancel = self._cancel = threading.Event()
        self._report = None
        self.open_report.Hide()
        self.gauge.SetValue(0)
        self._set_busy(True)
        self.pulse.Start(80)
        self._set_status("info", "Simulating at %g MHz (%s quality)…" % (mhz, QUALITIES[self.quality.GetSelection()][1]))
        self.say("Simulating at %g MHz, %s quality" % (mhz, quality))

        def prog(stage, frac):
            if cancel.is_set():
                raise _Cancelled()
            wx.CallAfter(self._ui, self._progress, stage, frac)

        def work():
            try:
                from ..core.report import html_report
                from ..core.solver.solve import simulate
                res = simulate(g, sim, sweep_freqs=sweep, progress=prog)
                if cancel.is_set():
                    raise _Cancelled()
                with open(path, "w", encoding="utf-8") as f:
                    f.write(html_report(g, res, len(check(g))))
                wx.CallAfter(self._ui, self._sim_done, cancel, res, path)
            except _Cancelled:
                pass
            except Exception as ex:
                wx.CallAfter(self._ui, self._sim_failed, cancel, ex, traceback.format_exc(limit=3))
        threading.Thread(target=work, daemon=True).start()

    def _ui(self, fn, *args):
        if self:  # the dialog may have been closed while the worker was running
            fn(*args)

    def _progress(self, stage, frac):
        if self._busy:
            if frac > 0:  # determinate while a stage reports fractions, indeterminate otherwise
                self.pulse.Stop()
                self.gauge.SetValue(int(frac * 100))
            elif not self.pulse.IsRunning():
                self.pulse.Start(80)
            self.status.SetLabel("Simulating: %s…" % stage)

    def on_cancel(self, e):
        if self._cancel:
            self._cancel.set()
        self._set_busy(False)
        self._set_status("info", "Simulation cancelled.")
        self.say("Simulation cancelled")

    def _sim_failed(self, cancel, ex, tb):
        if cancel is not self._cancel or cancel.is_set():
            return
        self._set_busy(False)
        self._set_status("error", "Simulation failed: %s" % humanize(str(ex)))
        self.say("Error: %s\n%s" % (ex, tb))

    def _sim_done(self, cancel, res, path):
        if cancel is not self._cancel or cancel.is_set():
            return
        self._set_busy(False)
        self._report = path
        summary = "L = %.4f µH · ESR = %.1f mΩ · Q = %.1f · Rac/Rdc = %.2f" % (
            res.L * 1e6, res.R * 1e3, res.Q, res.R / res.R_dc)
        self.say(summary)
        self.say("Ring current share: " + ", ".join("ring %s %.1f %%" % (k, abs(v) * 100)
                                                    for k, v in sorted(res.ring_share().items())))
        self.say("Report: " + path)
        self.open_report.Show()
        self._set_status("ok", summary)

    def _export(self, title, wildcard, default, fn):
        with wx.FileDialog(self, title, defaultDir=self._board_dir(), defaultFile=default, wildcard=wildcard,
                           style=wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT) as fd:
            if fd.ShowModal() == wx.ID_OK:
                def run():
                    text = fn(self._geo())
                    with open(fd.GetPath(), "w", encoding="utf-8") as f:
                        f.write(text)
                    self.say("Wrote " + fd.GetPath())
                    self._set_status("ok", "Exported %s." % os.path.basename(fd.GetPath()))
                self._safe(run)

    def on_fasthenry(self, e):
        mhz = self._safe(lambda: self._parse("freq"))
        if mhz is None:
            return
        name = self.ctrls["group_name"][0].GetValue().strip()
        self._export("Export FastHenry", "FastHenry input (*.inp)|*.inp", name + ".inp",
                     lambda g: fasthenry_text(g, SimParams(frequency=mhz * 1e6)))

    def on_openems(self, e):
        name = self.ctrls["group_name"][0].GetValue().strip()
        self._export("Export openEMS Script", "Python script (*.py)|*.py", name + "_openems.py",
                     lambda g: openems_script(g))

    def on_close(self, e):
        self.pulse.Stop()
        if self._cancel:
            self._cancel.set()
        if self.IsModal():
            self.EndModal(wx.ID_CLOSE)
        else:
            self.Destroy()
