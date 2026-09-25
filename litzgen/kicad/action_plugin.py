"""pcbnew ActionPlugin entry point."""
from __future__ import annotations

import os

import pcbnew


class LitzGenPlugin(pcbnew.ActionPlugin):
    def defaults(self):
        self.name = "LitzGen: multi-bundle PCB Litz coil"
        self.category = "Generate"
        self.description = "Generate, check and simulate transposed multi-layer PCB Litz coils"
        self.show_toolbar_button = True
        here = os.path.dirname(os.path.dirname(__file__))
        icon = os.path.join(here, "icon.png")
        if os.path.exists(icon):
            self.icon_file_name = icon

    def Run(self):
        import wx
        try:
            import numpy  # noqa: F401
        except ImportError:
            wx.MessageBox("LitzGen needs numpy in KiCad's Python.\n"
                          "Open the KiCad Command Prompt (Windows) or a terminal and run:\n"
                          "    python -m pip install numpy", "LitzGen", wx.ICON_ERROR)
            return
        from .dialog import LitzDialog
        board = pcbnew.GetBoard()
        parent = wx.FindWindowByName("PcbFrame")
        dlg = LitzDialog(parent, board)
        dlg.ShowModal()
        dlg.Destroy()
        pcbnew.Refresh()
