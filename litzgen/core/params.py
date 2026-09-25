"""Parameter model for multi-bundle PCB Litz coils.

All lengths are in millimetres, frequencies in hertz, conductivity in S/m.
Defaults reproduce the coil of Kale & Wicht, "A Dual-Bundle PCB Litz Coil
Achieving 1 kW, 6.78 MHz WPT with 97.8 % AC-AC Efficiency", WPTCE 2026
(Table I): d_out = 160 mm, 5 turns, 4 x 4 strand matrix of 0.8 mm traces
with 0.2 mm spacing (bundle width 3.8 mm), 70 um copper, 4-layer stackup
with 0.4 mm prepreg / 0.5 mm core, blind + buried adjacent-layer vias.
"""
from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, fields
from typing import List, Optional

MU0 = 4e-7 * math.pi
RHO_CU_20C = 1.72e-8          # ohm*m
ALPHA_CU = 0.00393            # 1/K


@dataclass
class Stackup:
    """Copper layers used by the coil, top to bottom.

    ``dielectric`` has one entry fewer than the number of copper layers and
    gives the dielectric thickness between consecutive copper layers.
    ``kicad_layers`` maps each coil layer to a KiCad copper layer name.
    """
    n_layers: int = 4
    copper_thickness: float = 0.070
    dielectric: List[float] = field(default_factory=lambda: [0.4, 0.5, 0.4])
    kicad_layers: List[str] = field(default_factory=lambda: ["F.Cu", "In1.Cu", "In2.Cu", "B.Cu"])

    def layer_z(self) -> List[float]:
        """z of each copper layer's mid-plane (mm), top layer at 0, downward negative."""
        z = [0.0]
        for i in range(1, self.n_layers):
            z.append(z[-1] - (self.copper_thickness + self.dielectric[i - 1]))
        return z

    @property
    def total_thickness(self) -> float:
        return self.n_layers * self.copper_thickness + sum(self.dielectric)

    def validate(self) -> None:
        if self.n_layers < 2:
            raise ValueError("at least two copper layers are needed for transposition")
        if len(self.dielectric) != self.n_layers - 1:
            raise ValueError("dielectric must have n_layers-1 entries")
        if len(self.kicad_layers) != self.n_layers:
            raise ValueError("kicad_layers must have n_layers entries")


@dataclass
class ViaSpec:
    """Adjacent-layer via (blind for outer pairs, buried for inner pairs)."""
    drill: float = 0.20
    pad: float = 0.45
    plating: float = 0.020            # barrel plating thickness (mm)
    filled: bool = True               # IPC-4761 type VI (epoxy filled)


@dataclass
class CoilParams:
    # --- coil envelope --------------------------------------------------
    d_out: float = 160.0              # outer diameter, outer edge of outer turn
    n_turns: float = 5.0
    turn_gap: Optional[float] = None  # edge-to-edge gap between turns; None -> derive from d_in
    d_in: Optional[float] = 69.0      # inner diameter, inner edge of inner turn (used if turn_gap is None)
    clockwise: bool = False           # spiral direction seen from the top (F.Cu)
    center_x: float = 100.0           # board position of coil centre (mm)
    center_y: float = 100.0
    start_angle_deg: float = 0.0      # angle of the outer terminal

    # --- strand matrix ----------------------------------------------------
    n_cols: int = 4                   # strands per layer across the bundle
    strand_width: float = 0.8
    strand_gap: float = 0.2
    stackup: Stackup = field(default_factory=Stackup)

    # --- transposition ----------------------------------------------------
    steps_per_turn: int = 24          # transposition zones per turn (0 = none, -1 = densest that fits)
    stagger_alternate_turns: bool = False
    jog_width: Optional[float] = None # trace width inside jogs (None = strand_width)
    clearance: float = 0.127          # copper-to-copper clearance rule used by the generator + DRC
    hole_clearance: float = 0.20      # hole-to-hole (edge) clearance
    via: ViaSpec = field(default_factory=ViaSpec)
    omit_rings: List[int] = field(default_factory=list)  # leave these concentric rings empty (e.g. [1] = hollow bundle)
    terminal_margin_deg: float = 6.0  # no transposition zones this close to a terminal
    terminal_length: float = 6.0      # arc length of the terminal bar (mm)
    terminal_drill: float = 1.2       # stitching / wire holes in terminal bars
    max_arc_deg: float = 10.0         # longest emitted arc track

    # --- KiCad ------------------------------------------------------------
    net_prefix: str = "LITZ1"
    per_strand_nets: bool = False     # True: one net per strand, no terminal bars (DRC-only boards)
    group_name: str = "LitzCoil"

    # --- derived helpers ----------------------------------------------------
    @property
    def pitch_lat(self) -> float:
        """Centre-to-centre distance of neighbouring strands in one layer."""
        return self.strand_width + self.strand_gap

    @property
    def bundle_width(self) -> float:
        return self.n_cols * self.strand_width + (self.n_cols - 1) * self.strand_gap

    @property
    def r_start(self) -> float:
        """Centre-line radius of the bundle at the outer terminal."""
        return self.d_out / 2.0 - self.bundle_width / 2.0

    @property
    def turn_pitch(self) -> float:
        """Radial advance of the spiral per turn (centre-line)."""
        if self.turn_gap is not None:
            return self.bundle_width + self.turn_gap
        if self.d_in is None:
            raise ValueError("give either turn_gap or d_in")
        r_end = self.d_in / 2.0 + self.bundle_width / 2.0
        return (self.r_start - r_end) / self.n_turns

    @property
    def effective_turn_gap(self) -> float:
        return self.turn_pitch - self.bundle_width

    @property
    def r_end(self) -> float:
        return self.r_start - self.turn_pitch * self.n_turns

    @property
    def derived_d_in(self) -> float:
        return 2.0 * (self.r_end - self.bundle_width / 2.0)

    @property
    def theta_end(self) -> float:
        return 2.0 * math.pi * self.n_turns

    @property
    def w_jog(self) -> float:
        return self.jog_width if self.jog_width is not None else self.strand_width

    def radius(self, theta: float) -> float:
        """Centre-line radius at unwrapped spiral angle theta (rad, 0 = outer terminal)."""
        return self.r_start - self.turn_pitch * theta / (2.0 * math.pi)

    def validate(self) -> List[str]:
        """Raise on impossible input; return a list of warnings."""
        self.stackup.validate()
        warn: List[str] = []
        if self.n_cols < 2:
            raise ValueError("n_cols must be >= 2")
        if self.turn_pitch <= self.bundle_width:
            raise ValueError("turns overlap: turn pitch %.3f <= bundle width %.3f" % (self.turn_pitch, self.bundle_width))
        if self.r_end - self.bundle_width / 2.0 - self.pitch_lat <= 1.0:
            raise ValueError("inner radius too small for the requested turns")
        p = self.pitch_lat
        if self.strand_gap < self.clearance:
            raise ValueError("strand_gap %.3f < clearance %.3f" % (self.strand_gap, self.clearance))
        if p < self.via.pad / 2 + self.strand_width / 2 + self.clearance:
            raise ValueError("lateral pitch too small to fit a via pad beside a strand")
        if p < self.via.pad + self.clearance:
            raise ValueError("lateral pitch too small for side-by-side via pads")
        if p < self.via.drill + self.hole_clearance:
            raise ValueError("lateral pitch too small for hole-to-hole clearance")
        if self.w_jog + self.clearance >= p:
            raise ValueError("jog width + clearance must be smaller than the lateral pitch")
        if self.effective_turn_gap < 2 * p + self.clearance:
            raise ValueError("turn gap %.2f mm leaves no room for the transposition lanes (need > %.2f mm)"
                             % (self.effective_turn_gap, 2 * p + self.clearance))
        if abs(self.n_turns * 2 - round(self.n_turns * 2)) > 1e-9:
            warn.append("non half-integer turn count: terminal angles are arbitrary")
        return warn

    # --- serialisation --------------------------------------------------------
    def to_dict(self) -> dict:
        return asdict(self)

    def to_json(self) -> str:
        return json.dumps(self.to_dict(), indent=2)

    @classmethod
    def from_dict(cls, d: dict) -> "CoilParams":
        d = dict(d)
        if "stackup" in d and isinstance(d["stackup"], dict):
            d["stackup"] = Stackup(**d["stackup"])
        if "via" in d and isinstance(d["via"], dict):
            d["via"] = ViaSpec(**d["via"])
        known = {f.name for f in fields(cls)}
        unknown = set(d) - known
        if unknown:
            raise ValueError("unknown parameters: %s" % ", ".join(sorted(unknown)))
        return cls(**d)

    @classmethod
    def from_json(cls, s: str) -> "CoilParams":
        return cls.from_dict(json.loads(s))


@dataclass
class SimParams:
    frequency: float = 6.78e6
    temperature_c: float = 25.0
    # 3D PEEC mesh
    max_seg_len: float = 3.0          # mm, longest straight filament segment
    near_factor: float = 1.5          # pairs closer than this x length get high-order quadrature
    mid_factor: float = 5.0           # pairs closer than this x length get 4-point quadrature
    # axisymmetric eddy-current mesh (per strand cross-section)
    eddy_nw: int = 10
    eddy_nt: int = 3
    eddy_angles: int = 3              # representative angles for the ring solves
    slice_deg: float = 0.5            # angular sampling of strand positions
    include_eddy: bool = True
    max_cache_mb: float = 1200.0      # memory allowed for cached filament matrices

    @classmethod
    def preset(cls, quality: str = "standard", **kw) -> "SimParams":
        q = {
            "fast": dict(eddy_nw=8, eddy_nt=3, eddy_angles=2, slice_deg=1.0, max_seg_len=4.0),
            "standard": dict(),
            "fine": dict(eddy_nw=14, eddy_nt=4, eddy_angles=4, slice_deg=0.25, max_seg_len=2.0),
        }[quality]
        q.update(kw)
        return cls(**q)

    def to_dict(self) -> dict:
        return asdict(self)

    @property
    def rho(self) -> float:
        return RHO_CU_20C * (1.0 + ALPHA_CU * (self.temperature_c - 20.0))

    def skin_depth(self) -> float:
        """Skin depth in mm."""
        return math.sqrt(self.rho / (math.pi * self.frequency * MU0)) * 1e3
