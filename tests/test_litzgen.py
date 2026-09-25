"""LitzGen test-suite (stdlib unittest; run: python -m unittest discover -s tests -v)."""
import json
import math
import os
import random
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from litzgen.core.drc import check, seg_seg_distance                      # noqa: E402
from litzgen.core.export.fasthenry import fasthenry_text                  # noqa: E402
from litzgen.core.export.kicad_sexpr import board_text, parse, terminal_footprint_text  # noqa: E402
from litzgen.core.export.openems import openems_script                    # noqa: E402
from litzgen.core.geometry import Dwell, Jog, ViaHop, generate, pos_xy   # noqa: E402
from litzgen.core.neighbors import query_pairs                            # noqa: E402
from litzgen.core.params import MU0, CoilParams, SimParams, Stackup       # noqa: E402
from litzgen.core.solver.eddy2d import (CrossSection, ellipke_agm, log_gmd_rect,  # noqa: E402
                                        ring_mutual, solve_cross_section)
from litzgen.core.solver.peec3d import SegmentMesh, self_partial, strand_inductance  # noqa: E402
from litzgen.core.transposition import build_rings, step_plan, verify_plan  # noqa: E402

RHO = 1.72e-8


def ring_ref(r1, r2, dz):
    """Coaxial ring mutual via numerical Neumann integral (independent of the AGM code)."""
    phi = np.linspace(0, 2 * np.pi, 20001)[:-1]
    d = np.sqrt(r1 * r1 + r2 * r2 - 2 * r1 * r2 * np.cos(phi) + dz * dz)
    return MU0 / (4 * np.pi) * 2 * np.pi * r1 * r2 * np.mean(np.cos(phi) / d) * 2 * np.pi


class TestTransposition(unittest.TestCase):
    def test_plan_is_clash_free_rotation(self):
        for C in range(2, 9):
            for R in range(2, 9):
                verify_plan(C, R)

    def test_paper_matrix(self):
        rings, fixed = build_rings(4, 4)
        self.assertEqual([r.size for r in rings], [12, 4])
        self.assertEqual(fixed, [])
        plan, nv = step_plan(4, 4)
        self.assertEqual(nv, 3)
        self.assertEqual(sum(m.via_index is not None for m in plan.values()), 8)   # 6 outer + 2 inner

    def test_full_rotation_returns_home(self):
        from litzgen.core.transposition import apply_step, initial_slots
        slots, _ = initial_slots(4, 4)
        plan, _ = step_plan(4, 4)
        cur = list(slots)
        for _ in range(12):
            cur = apply_step(cur, plan)
        self.assertEqual(cur, slots)


class TestGeometry(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.p = CoilParams()
        cls.g = generate(cls.p)

    def test_paper_dimensions(self):
        p = self.p
        self.assertAlmostEqual(p.bundle_width, 3.8)
        self.assertAlmostEqual(p.derived_d_in, 69.0)
        self.assertEqual(self.g.n_strands, 16)

    def test_paths_are_continuous(self):
        for pcs in self.g.pieces:
            last = None
            for pc in pcs:
                if isinstance(pc, ViaHop):
                    self.assertEqual(last[1], pc.la)
                    last = (pc.u, pc.lb, pc.t)
                    continue
                ua = pc.u if isinstance(pc, Dwell) else pc.ua
                ub = pc.u if isinstance(pc, Dwell) else pc.ub
                if last is not None:
                    self.assertAlmostEqual(last[0], ua)
                    self.assertEqual(last[1], pc.layer)
                    self.assertAlmostEqual(last[2], pc.ta, places=12)
                last = (ub, pc.layer, pc.tb)
            self.assertAlmostEqual(last[2], self.p.theta_end, places=9)

    def test_vias_only_adjacent_layers(self):
        for v in self.g.vias:
            self.assertEqual(abs(v.la - v.lb), 1)

    def test_strand_lengths_matched(self):
        L = self.g.strand_lengths()
        self.assertLess((L.max() - L.min()) / L.mean(), 0.005)

    def test_drc_clean_default(self):
        self.assertEqual(check(self.g), [])

    def test_drc_detects_short(self):
        p = CoilParams(strand_gap=0.2, clearance=0.127)
        g = generate(p)
        g.arcs[0].width = 1.2          # widen one arc into its neighbour
        v = check(g)
        self.assertTrue(any(x.kind == "clearance" for x in v))

    def test_randomised_drc(self):
        rnd = random.Random(7)
        accepted = 0
        for _ in range(40):
            R = rnd.choice([2, 4, 6]); C = rnd.choice([2, 3, 4, 5])
            w = rnd.choice([0.5, 0.8, 1.2]); gap = rnd.choice([0.2, 0.3])
            kl = ["F.Cu"] + ["In%d.Cu" % i for i in range(1, R - 1)] + ["B.Cu"]
            p = CoilParams(d_out=rnd.choice([100, 160]), n_turns=rnd.choice([2, 3, 4.5]), d_in=None,
                           turn_gap=rnd.choice([4, 6, 8]), n_cols=C, strand_width=w, strand_gap=gap,
                           stackup=Stackup(R, 0.035, [0.2] * (R - 1), kl), steps_per_turn=rnd.choice([8, 12, 16]),
                           clockwise=rnd.random() < 0.5, stagger_alternate_turns=rnd.random() < 0.3)
            try:
                g = generate(p)
            except ValueError:
                continue
            accepted += 1
            self.assertEqual(check(g), [], msg=p.to_json())
        self.assertGreater(accepted, 10)

    def test_invalid_parameters_rejected(self):
        with self.assertRaises(ValueError):
            generate(CoilParams(strand_gap=0.05))
        with self.assertRaises(ValueError):
            generate(CoilParams(steps_per_turn=200))

    def test_auto_steps(self):
        g = generate(CoilParams(n_turns=3, steps_per_turn=-1))
        self.assertGreater(len(g.zones), 0)
        self.assertEqual(check(g), [])

    def test_omit_ring(self):
        g = generate(CoilParams(omit_rings=[1]))
        self.assertEqual(g.n_strands, 12)
        self.assertEqual(check(g), [])

    def test_segment_distance(self):
        P0 = np.array([[0, 0], [0, 0], [0, 0]], float); P1 = np.array([[1, 0], [1, 0], [1, 1]], float)
        Q0 = np.array([[0, 1], [2, 0], [1, 0]], float); Q1 = np.array([[1, 1], [3, 0], [0, 1]], float)
        d, _ = seg_seg_distance(P0, P1, Q0, Q1)
        np.testing.assert_allclose(d, [1, 1, 0], atol=1e-12)

    def test_neighbors_numpy_matches(self):
        pts = np.random.default_rng(3).uniform(0, 5, (800, 3))
        a = {tuple(x) for x in query_pairs(pts, 0.6)}
        b = {tuple(x) for x in query_pairs(pts, 0.6, force_numpy=True)}
        self.assertEqual(a, b)


class TestExports(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.g = generate(CoilParams())

    def test_board_roundtrip(self):
        txt = board_text(self.g)
        tree = parse(txt)[0]
        self.assertEqual(tree[0], "kicad_pcb")
        kinds = [x[0] for x in tree if isinstance(x, list)]
        self.assertEqual(kinds.count("via"), len(self.g.vias))
        self.assertEqual(kinds.count("arc"), len(self.g.arcs))
        self.assertEqual(kinds.count("footprint"), 2)
        vias = [x for x in tree if isinstance(x, list) and x[0] == "via"]
        pairs = {tuple(y[1:]) for v in vias for y in v if isinstance(y, list) and y[0] == "layers"}
        self.assertEqual(pairs, {("F.Cu", "In1.Cu"), ("In1.Cu", "In2.Cu"), ("In2.Cu", "B.Cu")})
        txt_items = [x for x in tree if isinstance(x, list) and x[0] == "gr_text"]
        rec = json.loads(txt_items[0][1][len("LITZGEN "):])
        self.assertEqual(CoilParams.from_dict(rec).to_dict(), self.g.params.to_dict())

    def test_footprint_parses(self):
        tree = parse(terminal_footprint_text(self.g.terminals[0]))[0]
        self.assertEqual(tree[0], "footprint")

    def test_fasthenry(self):
        t = fasthenry_text(self.g)
        lines = t.splitlines()
        self.assertTrue(any(x.startswith(".external") for x in lines))
        self.assertEqual(sum(1 for x in lines if x.startswith(".equiv")), 2)
        nodes = {x.split()[0] for x in lines if x.startswith("N")}
        for x in lines:
            if x.startswith("E"):
                a, b = x.split()[1:3]
                self.assertIn(a, nodes); self.assertIn(b, nodes)

    def test_openems_compiles(self):
        compile(openems_script(self.g), "openems_script", "exec")


class TestSolverPhysics(unittest.TestCase):
    def test_elliptic_against_numeric_neumann(self):
        for r1, r2, dz in [(0.05, 0.051, 0.0), (0.05, 0.08, 0.002), (0.03, 0.03, 0.0005)]:
            self.assertAlmostEqual(float(ring_mutual(np.array(r1), 0, np.array(r2), dz)) / ring_ref(r1, r2, dz), 1, places=4)

    def test_self_partial_matches_grover(self):
        l, w, t = 0.05, 0.8e-3, 0.07e-3
        grover = MU0 * l / (2 * np.pi) * (np.log(2 * l / (w + t)) + 0.5 + (w + t) / (3 * l))
        self.assertAlmostEqual(float(self_partial(np.array(l), np.array(0.2235 * (w + t)))) / grover, 1, delta=2e-3)

    def test_gmd_square(self):
        one = np.array([0.0]); u = np.array([1.0])
        self.assertAlmostEqual(float(np.exp(log_gmd_rect(one, u, one, u, one, u, one, u))[0]), 0.44705, places=4)

    def _loop(self, R, z, n, strand):
        th = np.linspace(0, 2 * np.pi, n + 1)
        P = np.stack([R * np.cos(th), R * np.sin(th), np.full_like(th, z)], 1)
        return P[:-1], P[1:], np.full(n, strand), np.full(n, 0.8e-3), np.full(n, 0.2235 * 0.87e-3), np.zeros(n, bool)

    def test_peec_circular_loop(self):
        a = self._loop(0.05, 0.0, 180, 0)
        b = self._loop(0.05, 0.001, 180, 1)
        m = SegmentMesh(*[np.concatenate([x, y]) for x, y in zip(a, b)])
        L = strand_inductance(m, 2)
        L_self = MU0 * 0.05 * (np.log(8 * 0.05 / (0.2235 * 0.87e-3)) - 2)
        self.assertAlmostEqual(L[0, 0] / L_self, 1, delta=2e-3)
        # strip-strip mutual: width-averaged ring formula
        x = (np.arange(30) + 0.5) / 30 * 0.8e-3 - 0.4e-3
        ref = np.mean(ring_mutual(0.05 + x[:, None], 0, 0.05 + x[None, :], 0.001))
        self.assertAlmostEqual(L[0, 1] / ref, 1, delta=3e-3)

    def test_eddy_dc_limit(self):
        cs = CrossSection(np.array([0.05]), np.array([0.0]), np.array([0.8e-3]), np.array([70e-6]))
        r = solve_cross_section(cs, 1.0, RHO, 10, 3)
        self.assertAlmostEqual(r.Zc[0, 0].real / (RHO * 0.05 / (0.8e-3 * 70e-6)), 1, places=4)

    def test_eddy_dowell_coax_foils(self):
        f = 6.78e6
        delta = math.sqrt(RHO / (math.pi * f * MU0))
        for D, tol in [(1.0, 0.01), (2.0, 0.03)]:
            h, H, g = D * delta, 0.1, 0.1e-3
            cs = CrossSection(np.array([0.05, 0.05 + h + g]), np.zeros(2), np.array([h, h]), np.array([H, H]))
            res = solve_cross_section(cs, f, RHO, nw=10, nt=80)
            I = np.array([1.0, -1.0])
            P = 0.5 * np.real(I @ res.Zc @ I)
            Rdc = RHO * 0.05 / (h * H) + RHO * (0.05 + h + g) / (h * H)
            F = D * (math.sinh(2 * D) + math.sin(2 * D)) / (math.cosh(2 * D) - math.cos(2 * D))
            self.assertAlmostEqual((2 * P / Rdc) / F, 1, delta=tol)


class TestEndToEndFast(unittest.TestCase):
    def test_small_coil_simulation(self):
        from litzgen.core.solver.solve import simulate
        p = CoilParams(d_out=80, n_turns=2, d_in=None, turn_gap=5, steps_per_turn=12)
        g = generate(p)
        res = simulate(g, SimParams.preset("fast"), keep_cross_section=True)
        self.assertGreater(res.Q, 10)
        self.assertAlmostEqual(abs(np.sum(res.strand_currents)), 1.0, places=9)
        self.assertGreater(res.R, res.R_dc)
        # two-turn spiral ~ within 15 % of the concentric-ring model
        self.assertAlmostEqual(res.L / res.L_axisym, 1.0, delta=0.15)
        from litzgen.core.report import html_report
        self.assertIn("<svg", html_report(g, res, 0))


if __name__ == "__main__":
    unittest.main()


class TestPcbnewWriterWithFake(unittest.TestCase):
    def test_place_regenerate_remove(self):
        self._place_regenerate_remove(kicad10_groups=False)

    def test_place_regenerate_remove_kicad10_groups(self):
        # KiCad 10: GetParentGroup() returns an EDA_GROUP without m_Uuid
        self._place_regenerate_remove(kicad10_groups=True)

    def _place_regenerate_remove(self, kicad10_groups):
        import importlib
        sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
        fake = importlib.import_module("fake_pcbnew")
        fake.EDA_GROUP_PARENTS = kicad10_groups
        sys.modules["pcbnew"] = fake
        try:
            W = importlib.import_module("litzgen.kicad.pcbnew_writer")
            importlib.reload(W)
            b = fake.BOARD()
            g = generate(CoilParams())
            n = W.place_coil(b, g)
            self.assertEqual(n["vias"], len(g.vias))
            self.assertEqual(n["terminals"], 2)
            self.assertEqual(b.GetCopperLayerCount(), 4)
            vias = [i for i in b.items if isinstance(i, fake.PCB_VIA)]
            self.assertEqual({v.pair for v in vias}, {(0, 1), (1, 2), (2, 31)})
            self.assertTrue(all(v.vtype == fake.VIATYPE_BLIND_BURIED for v in vias))
            rec = W.read_params_from_board(b, "LitzCoil")
            self.assertEqual(rec, g.params.to_dict())
            total = len(b.items)
            W.place_coil(b, g, replace=True)            # regenerate replaces, does not duplicate
            self.assertEqual(len(b.items), total)
            self.assertEqual(W.list_coils(b), ["LitzCoil"])
            removed = W.remove_coil(b, "LitzCoil")
            self.assertEqual(removed, total - 1)
            self.assertEqual(b.items, [])
        finally:
            fake.EDA_GROUP_PARENTS = False
            del sys.modules["pcbnew"]


class TestPcmMetadata(unittest.TestCase):
    """metadata.json must satisfy KiCad's PCM schema (https://go.kicad.org/pcm/schemas/v1),
    or Plugin and Content Manager refuses the package ("Unable to parse package metadata")."""

    def test_schema_patterns(self):
        import re
        root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        meta = json.load(open(os.path.join(root, "metadata.json"), encoding="utf-8"))
        for key in ("name", "description", "description_full", "identifier", "type",
                    "author", "license", "resources", "versions"):
            self.assertIn(key, meta)
        self.assertRegex(meta["identifier"], r"^[a-zA-Z][-a-zA-Z0-9.]{0,98}[a-zA-Z0-9]$")
        self.assertTrue(meta["tags"])
        self.assertEqual(len(meta["tags"]), len(set(meta["tags"])))
        for tag in meta["tags"]:
            self.assertRegex(tag, r"^[a-z][-a-z0-9]{0,48}[a-z0-9]$")
        from litzgen import __version__
        self.assertRegex(__version__, r"^\d{1,4}(\.\d{1,4}(\.\d{1,6})?)?$")
