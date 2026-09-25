"""Minimal stand-in for the pcbnew SWIG module, used to exercise the board writer
without KiCad. It records calls and enforces the method names LitzGen relies on."""
import itertools
import uuid as _uuid

from litzgen.core.export.kicad_sexpr import parse

_ids = itertools.count()
F_Cu, B_Cu, Cmts_User = 0, 31, 41
VIATYPE_THROUGH, VIATYPE_BLIND_BURIED, VIATYPE_MICROVIA = 3, 2, 1


def FromMM(x):
    return int(round(x * 1e6))


class VECTOR2I:
    def __init__(self, x, y):
        self.x, self.y = x, y


class KIID:
    def __init__(self):
        self._s = str(_uuid.uuid4())

    def AsString(self):
        return self._s


class LSET:
    @staticmethod
    def AllCuMask(n):
        return set(range(n))


class _Item:
    def __init__(self, board=None):
        self.m_Uuid = KIID()
        self._group = None
        self.net = None
        self.layer = None

    def GetParentGroup(self):
        return self._group

    def SetLayer(self, l):
        self.layer = l

    def SetNet(self, n):
        self.net = n

    def Cast(self):
        return self


class PCB_TRACK(_Item):
    def SetStart(self, v): self.start = v
    def SetEnd(self, v): self.end = v
    def SetWidth(self, w): self.width = w


class PCB_ARC(PCB_TRACK):
    def SetMid(self, v): self.mid = v


class PCB_VIA(_Item):
    def SetPosition(self, v): self.pos = v
    def SetViaType(self, t): self.vtype = t
    def SetLayerPair(self, a, b): self.pair = (a, b)
    def SetDrill(self, d): self.drill = d
    def SetWidth(self, w): self.width = w


class PCB_TEXT(_Item):
    def SetText(self, t): self.text = t
    def GetText(self): return self.text
    def SetPosition(self, v): self.pos = v
    def SetTextSize(self, v): self.size = v


class PCB_GROUP(_Item):
    def __init__(self, board):
        super().__init__(board)
        self.items = []
        self.name = ""

    def SetName(self, n): self.name = n
    def GetName(self): return self.name

    def AddItem(self, it):
        it._group = self
        self.items.append(it)

    def RemoveItem(self, it):
        self.items.remove(it)
        it._group = None


class NETINFO_ITEM:
    def __init__(self, board, name):
        self.name = name


class PAD(_Item):
    pass


class FOOTPRINT(_Item):
    def __init__(self, n_pads):
        super().__init__()
        self._pads = [PAD() for _ in range(n_pads)]

    def Pads(self): return self._pads
    def SetPosition(self, v): self.pos = v
    def SetReference(self, r): self.ref = r
    def SetValue(self, v): self.value = v


def FootprintLoad(lib, name):
    import os
    tree = parse(open(os.path.join(lib, name + ".kicad_mod")).read())[0]
    assert tree[0] == "footprint"
    return FOOTPRINT(sum(1 for x in tree if isinstance(x, list) and x[0] == "pad"))


def Refresh():
    pass


class BOARD:
    def __init__(self):
        self.items = []
        self.nets = {}
        self.cu = 2
        self.names = {"F.Cu": 0, "B.Cu": 31}

    def GetCopperLayerCount(self): return self.cu

    def SetCopperLayerCount(self, n):
        self.cu = n
        for i in range(1, n - 1):
            self.names["In%d.Cu" % i] = i

    def GetEnabledLayers(self): return set()
    def SetEnabledLayers(self, s): pass
    def GetLayerID(self, name): return self.names[name]
    def FindNet(self, name): return self.nets.get(name)

    def Add(self, it):
        if isinstance(it, NETINFO_ITEM):
            self.nets[it.name] = it
        else:
            self.items.append(it)

    def Remove(self, it): self.items.remove(it)
    def Groups(self): return [i for i in self.items if isinstance(i, PCB_GROUP)]
    def GetTracks(self): return [i for i in self.items if isinstance(i, (PCB_TRACK, PCB_VIA))]
    def GetFootprints(self): return [i for i in self.items if isinstance(i, FOOTPRINT)]
    def GetDrawings(self): return [i for i in self.items if isinstance(i, PCB_TEXT)]
    def BuildConnectivity(self): pass
    def GetFileName(self): return ""
