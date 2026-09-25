"""LitzGen: generator and solver for transposed multi-bundle PCB Litz coils.

Inside KiCad's pcbnew this package registers an action plugin; elsewhere it
is a plain library / CLI (``python -m litzgen --help``).
"""
__version__ = "1.0.0"

try:  # register only inside pcbnew with a GUI
    import pcbnew  # noqa: F401
    import wx  # noqa: F401
    if hasattr(pcbnew, "ActionPlugin"):
        from .kicad.action_plugin import LitzGenPlugin
        LitzGenPlugin().register()
except Exception:  # pragma: no cover - not running inside KiCad
    pass
