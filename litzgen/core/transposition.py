"""Concentric-ring transposition of an R x C strand matrix.

Cross-section model
-------------------
The bundle cross-section is a matrix of *slots* ``(u, layer)``: ``u`` is the
lateral index across the bundle (0 = innermost radius, C-1 = outermost) and
``layer`` the copper layer (0 = top).  The matrix is peeled into concentric
rings: ring 0 is the perimeter (12 slots for 4 x 4), ring 1 the next
perimeter inside it (the 2 x 2 "inner bundle" for 4 x 4), and so on.  This is
the "concentric multi-bundle" architecture of Kale & Wicht generalised to any
layer / column count.

A transposition *step* rotates every ring by one slot along its cycle
(top row -> right column -> bottom row -> left column).  A full ring rotation
is not realisable as a sequence of single moves (every slot is occupied, so
each move waits for another: a deadlock).  The step therefore uses *lanes*:

* ring 0's lanes are one pitch outside the bundle (in the turn gap);
* ring k's lanes are ring k-1's column slots, which ring k-1 vacates first.

Each step has three phases along the path:

X  (exit)     every strand that must change layer jogs one pitch sideways
              into its lane; at the same time the rows shift sideways.
              On every layer the moving groups diverge, so the jogs are
              clash-free parallel shifts.
V  (vias)     strands in a lane change layer with adjacent-layer vias
              (blind L1-L2 / L3-L4, buried inner pairs).  Vias in one lane
              are sequenced bottom-first (right lane) / top-first (left
              lane) so each via lands in a lane slot that is already empty.
R  (re-entry) lane strands jog back into the now-vacant column slots.

Per step, ring k places 2*(R-2k-1) vias; a 4 x 4 matrix therefore uses
6 (outer bundle) + 2 (inner bundle) = 8 vias per step.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

Slot = Tuple[int, int]   # (u, layer)


@dataclass(frozen=True)
class SlotMove:
    """What a strand sitting in a given slot does during one step."""
    x_du: int                      # lateral shift in phase X (-1, 0, +1)
    via_index: Optional[int]       # position of its via in the V sequence (None = no via)
    dl: int                        # layer change in phase V (-1, 0, +1)
    r_du: int                      # lateral shift in phase R


@dataclass
class Ring:
    index: int
    cycle: List[Slot]              # slots in rotation order

    @property
    def size(self) -> int:
        return len(self.cycle)

    def successor(self) -> Dict[Slot, Slot]:
        n = len(self.cycle)
        return {self.cycle[i]: self.cycle[(i + 1) % n] for i in range(n)}


def build_rings(n_cols: int, n_layers: int) -> Tuple[List[Ring], List[Slot]]:
    """Peel the matrix into rings. Returns (rings, fixed_slots).

    fixed_slots are left-over single rows/columns in the centre that cannot be
    transposed without extra layers; strands there run straight.
    """
    C, R = n_cols, n_layers
    rings: List[Ring] = []
    covered = set()
    k = 0
    while C - 2 * k >= 2 and R - 2 * k >= 2:
        top = [(c, k) for c in range(k, C - k)]
        right = [(C - 1 - k, l) for l in range(k + 1, R - k)]
        bottom = [(c, R - 1 - k) for c in range(C - 2 - k, k - 1, -1)]
        left = [(k, l) for l in range(R - 2 - k, k, -1)]
        cyc = top + right + bottom + left
        rings.append(Ring(k, cyc))
        covered.update(cyc)
        k += 1
    fixed = [(c, l) for l in range(R) for c in range(C) if (c, l) not in covered]
    return rings, fixed


def step_plan(n_cols: int, n_layers: int) -> Tuple[Dict[Slot, SlotMove], int]:
    """Slot-level choreography of one step. Returns (plan, n_via_slots)."""
    C, R = n_cols, n_layers
    rings, fixed = build_rings(C, R)
    plan: Dict[Slot, SlotMove] = {}
    n_via = 0
    for ring in rings:
        k = ring.index
        rows = R - 2 * k
        n_via = max(n_via, rows - 1)
        # top row: shift +1 in X (top-right corner goes to the lane)
        for c in range(k, C - k):
            plan[(c, k)] = SlotMove(+1, None, 0, 0)
        # bottom row: shift -1 in X (bottom-left corner goes to the lane)
        for c in range(k, C - k):
            plan[(c, R - 1 - k)] = SlotMove(-1, None, 0, 0)
        # right movers: (C-1-k, l) for l = k .. R-2-k descend; order bottom-first
        for l in range(k, R - 1 - k):
            idx = (R - 2 - k) - l
            plan[(C - 1 - k, l)] = SlotMove(+1, idx, +1, -1)
        # left movers: (k, l) for l = R-1-k .. k+1 ascend; order top-first
        for l in range(k + 1, R - k):
            idx = l - (k + 1)
            plan[(k, l)] = SlotMove(-1, idx, -1, +1)
    for s in fixed:
        plan[s] = SlotMove(0, None, 0, 0)
    return plan, n_via


def apply_step(slots: List[Slot], plan: Dict[Slot, SlotMove]) -> List[Slot]:
    out = []
    for (u, l) in slots:
        m = plan[(u, l)]
        out.append((u + m.x_du + m.r_du, l + m.dl))
    return out


def phase_occupancy(slots: List[Slot], plan: Dict[Slot, SlotMove], n_via: int):
    """Yield (label, positions) snapshots through one step; used for self-checks."""
    pos = list(slots)
    yield "start", list(pos)
    pos = [(u + plan[s].x_du, l) for s, (u, l) in zip(slots, pos)]
    yield "after X", list(pos)
    for i in range(n_via):
        pos = [(u, l + (plan[s].dl if plan[s].via_index == i else 0)) for s, (u, l) in zip(slots, pos)]
        yield "after via %d" % i, list(pos)
    pos = [(u + plan[s].r_du, l) for s, (u, l) in zip(slots, pos)]
    yield "after R", list(pos)


def initial_slots(n_cols: int, n_layers: int) -> Tuple[List[Slot], List[int]]:
    """Strand numbering: ring 0 in cycle order, then ring 1, ..., then fixed slots.

    Returns (slot per strand, ring index per strand (-1 = fixed)).
    """
    rings, fixed = build_rings(n_cols, n_layers)
    slots: List[Slot] = []
    ring_of: List[int] = []
    for r in rings:
        slots += r.cycle
        ring_of += [r.index] * r.size
    slots += fixed
    ring_of += [-1] * len(fixed)
    return slots, ring_of


def verify_plan(n_cols: int, n_layers: int) -> None:
    """Raise AssertionError if the choreography is not a clash-free ring rotation."""
    rings, _ = build_rings(n_cols, n_layers)
    plan, n_via = step_plan(n_cols, n_layers)
    slots, _ = initial_slots(n_cols, n_layers)
    for label, pos in phase_occupancy(slots, plan, n_via):
        assert len(set(pos)) == len(pos), "slot clash %s: %s" % (label, pos)
    nxt = apply_step(slots, plan)
    succ = {}
    for r in rings:
        succ.update(r.successor())
    for s, e in zip(slots, nxt):
        assert succ.get(s, s) == e, "slot %s goes to %s, expected %s" % (s, e, succ.get(s, s))
