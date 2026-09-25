"""Place a generated coil on the open board through the pcbnew (SWIG) API.

Written against KiCad 8 and 9; every call that changed between versions
goes through a small adapter below.  Terminal bars are built by writing the
footprint as a .kicad_mod file and loading it with FootprintLoad, which
avoids the version-specific pad / padstack API entirely.
"""
from __future__ import annotations

import json
import os
import tempfile
from typing import Dict, List, Optional

import pcbnew  # noqa: F401  (only importable inside KiCad)

from ..core.export.kicad_sexpr import strand_net_names, terminal_footprint_text
from ..core.geometry import CoilGeometry

PARAM_TAG = "LITZGEN "


# ----------------------------------------------------------------------------- version adapters
def _v(x_mm: float, y_mm: float):
    return pcbnew.VECTOR2I(pcbnew.FromMM(x_mm), pcbnew.FromMM(y_mm))


def _via_type_blind_buried(lo_is_outer: bool):
    # KiCad 8/9: one enum for blind and buried; newer versions split them
    if hasattr(pcbnew, "VIATYPE_BLIND_BURIED"):
        return pcbnew.VIATYPE_BLIND_BURIED
    if lo_is_outer and hasattr(pcbnew, "VIATYPE_BLIND"):
        return pcbnew.VIATYPE_BLIND
    if hasattr(pcbnew, "VIATYPE_BURIED"):
        return pcbnew.VIATYPE_BURIED
    return pcbnew.VIATYPE_BLIND


def _set_via_size(via, pad_mm: float, drill_mm: float) -> None:
    via.SetDrill(pcbnew.FromMM(drill_mm))
    w = pcbnew.FromMM(pad_mm)
    try:
        via.SetWidth(w)
    except TypeError:                     # padstack-era signature SetWidth(width, layer)
        via.SetWidth(w, pcbnew.F_Cu)


def _ensure_copper_layers(board, names: List[str]) -> None:
    need = 2
    for nm in names:
        if nm.startswith("In"):
            need = max(need, int(nm[2:-3]) + 2)
    if board.GetCopperLayerCount() < need:
        board.SetCopperLayerCount(need)
        try:
            board.SetEnabledLayers(board.GetEnabledLayers() | pcbnew.LSET.AllCuMask(need))
        except Exception:
            pass


def _get_net(board, name: str):
    net = board.FindNet(name)
    if net is None:
        net = pcbnew.NETINFO_ITEM(board, name)
        board.Add(net)
    return net


# ----------------------------------------------------------------------------- remove / find
def find_group(board, name: str):
    for g in board.Groups():
        if g.GetName() == name:
            return g
    return None


def _cast(item):
    try:
        return item.Cast()
    except Exception:
        return item


def _uuid(obj) -> str:
    """UUID of a group or item. KiCad 10's GetParentGroup() returns an EDA_GROUP, which has
    no m_Uuid of its own; the board item behind it is AsEdaItem()."""
    uid = getattr(obj, "m_Uuid", None)
    if uid is None:
        uid = obj.AsEdaItem().m_Uuid
    return uid.AsString()


def _group_members(board, group) -> list:
    """Items of a group (via GetParentGroup, which is reliable across SWIG versions)."""
    out = []
    gid = _uuid(group)
    pools = [board.GetTracks(), board.GetFootprints(), board.GetDrawings()]
    for pool in pools:
        for it in pool:
            try:
                pg = it.GetParentGroup()
            except Exception:
                pg = None
            if pg is not None and _uuid(pg) == gid:
                out.append(it)
    return out


def read_params_from_board(board, group_name: str) -> Optional[dict]:
    g = find_group(board, group_name)
    items = _group_members(board, g) if g is not None else list(board.GetDrawings())
    for it in items:
        it = _cast(it)
        get = getattr(it, "GetText", None)
        if get is None:
            continue
        txt = get()
        if isinstance(txt, str) and txt.startswith(PARAM_TAG):
            return json.loads(txt[len(PARAM_TAG):])
    return None


def list_coils(board) -> List[str]:
    return [g.GetName() for g in board.Groups() if read_params_from_board(board, g.GetName()) is not None]


def remove_coil(board, group_name: str) -> int:
    g = find_group(board, group_name)
    if g is None:
        return 0
    items = _group_members(board, g)
    try:
        g.RemoveAll()
    except AttributeError:
        for it in items:
            try:
                g.RemoveItem(it)
            except Exception:
                pass
    # Delete, not Remove: after Remove() Python owns each item, and in KiCad 10 freeing those
    # wrappers corrupts pcbnew's SWIG state (later calls return bare SwigPyObjects).
    for it in items:
        board.Delete(it)
    board.Delete(g)
    return len(items)


# ----------------------------------------------------------------------------- place
def place_coil(board, geo: CoilGeometry, replace: bool = True) -> Dict[str, int]:
    p = geo.params
    lay_names = list(p.stackup.kicad_layers)
    _ensure_copper_layers(board, lay_names)
    lid = [board.GetLayerID(nm) for nm in lay_names]
    if replace:
        remove_coil(board, p.group_name)
    group = pcbnew.PCB_GROUP(board)
    group.SetName(p.group_name)
    board.Add(group)

    nets = strand_net_names(geo)
    net_obj = {nm: _get_net(board, nm) for nm in set(nets)}
    count = {"tracks": 0, "arcs": 0, "vias": 0, "terminals": 0}

    for t in geo.tracks:
        tr = pcbnew.PCB_TRACK(board)
        tr.SetStart(_v(*t.start)); tr.SetEnd(_v(*t.end))
        tr.SetWidth(pcbnew.FromMM(t.width)); tr.SetLayer(lid[t.layer]); tr.SetNet(net_obj[nets[t.strand]])
        board.Add(tr); group.AddItem(tr); count["tracks"] += 1
    for a in geo.arcs:
        arc = pcbnew.PCB_ARC(board)
        arc.SetStart(_v(*a.start)); arc.SetMid(_v(*a.mid)); arc.SetEnd(_v(*a.end))
        arc.SetWidth(pcbnew.FromMM(a.width)); arc.SetLayer(lid[a.layer]); arc.SetNet(net_obj[nets[a.strand]])
        board.Add(arc); group.AddItem(arc); count["arcs"] += 1
    n = len(lay_names)
    for v in geo.vias:
        lo, hi = sorted((v.la, v.lb))
        via = pcbnew.PCB_VIA(board)
        via.SetPosition(_v(*v.pos))
        if lo == 0 and hi == n - 1 and lay_names[0] == "F.Cu" and lay_names[-1] == "B.Cu":
            via.SetViaType(pcbnew.VIATYPE_THROUGH)
        else:
            via.SetViaType(_via_type_blind_buried(lo == 0 or hi == n - 1))
        via.SetLayerPair(lid[lo], lid[hi])
        _set_via_size(via, v.pad, v.drill)
        via.SetNet(net_obj[nets[v.strand]])
        board.Add(via); group.AddItem(via); count["vias"] += 1

    if geo.terminals:
        tmp = tempfile.mkdtemp(prefix="litzgen_")
        lib = os.path.join(tmp, "LitzGen.pretty")
        os.makedirs(lib, exist_ok=True)
        for k, t in enumerate(geo.terminals):
            name = "LitzTerminal_%s" % t.name
            with open(os.path.join(lib, name + ".kicad_mod"), "w") as f:
                f.write(terminal_footprint_text(t, name))
            fp = pcbnew.FootprintLoad(lib, name)
            fp.SetPosition(_v(*t.anchor))
            fp.SetReference("LT%d" % (k + 1))
            fp.SetValue("%s_%s" % (p.net_prefix, t.name))
            for pad in fp.Pads():
                pad.SetNet(net_obj[nets[0]])
            board.Add(fp); group.AddItem(fp); count["terminals"] += 1

    # parameter record so the coil can be re-opened and regenerated
    txt = pcbnew.PCB_TEXT(board)
    txt.SetText(PARAM_TAG + json.dumps(p.to_dict(), separators=(",", ":")))
    txt.SetLayer(pcbnew.Cmts_User)
    txt.SetPosition(_v(p.center_x - p.d_out / 2, p.center_y + p.d_out / 2 + 4))
    try:
        txt.SetTextSize(pcbnew.VECTOR2I(pcbnew.FromMM(0.5), pcbnew.FromMM(0.5)))
    except Exception:
        pass
    board.Add(txt); group.AddItem(txt)
    try:
        board.BuildConnectivity()
    except Exception:
        pass
    pcbnew.Refresh()
    return count
