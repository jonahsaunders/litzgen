# LitzGen design notes

LitzGen is a KiCad action plugin and a command-line tool. It generates concentric multi-bundle PCB Litz coils of the kind described by Kale & Wicht (WPTCE 2026), for any number of layers and columns, then checks and simulates them. For the paper's 160 mm, 5-turn, 4 × 4 coil it produces 7,600 arcs and 944 blind/buried vias with zero clearance violations. It predicts 3.475 µH (measured: 3.44 µH), 345 mΩ ESR (the paper's HFSS gives 383 mΩ) and Q = 429 at 6.78 MHz.

## Summary

- **Generates** the coil as native KiCad tracks, arcs, adjacent-layer vias and terminal bars, all in one group that can be regenerated. A lane-based step sequence rotates every concentric ring without deadlock, and all spacing comes from your design rules.
- **Checks** strand-to-strand clearance itself. Every strand is on the same net, so KiCad's DRC can't see a short between two strands.
- **Simulates** with two models: a 3D partial-inductance model for L and strand current sharing, and an axisymmetric eddy-current model of the full cross-section for AC resistance. It reports L, ESR, Q, per-strand currents, a loss split and a current-density map.
- **Exports** FastHenry models (strand level) and openEMS models (bundle level, for self-resonance).

The model also found something worth checking on the bench. At 6.78 MHz the inner 2 × 2 bundle carries only about 4 % of the current, and with 0.8 mm-wide strands, transposition does not lower ESR. See [Validation](#validation).

## Architecture

All physics and geometry live in a pure-Python core (numpy only) that never imports KiCad. Thin front-ends turn its output into board items. That keeps the solver testable without KiCad, and lets the same core serve the SWIG action plugin now and the IPC API later.

| Layer | Module | Job |
| --- | --- | --- |
| Core | `litzgen.core.params` | Coil, strand, stackup and via parameters (defaults = the paper's Table I) |
| Core | `litzgen.core.transposition` | Slot matrix → concentric rings → permutation and via schedule for each step |
| Core | `litzgen.core.geometry` | Maps (angle, lateral slot, layer) paths onto the spiral; emits tracks, arcs, vias, terminals |
| Core | `litzgen.core.drc` | Its own clearance checker, needed because every strand shares one net |
| Core | `litzgen.core.solver` | 3D PEEC for inductance and current sharing; axisymmetric eddy solve for AC resistance |
| Core | `litzgen.core.export` | `.kicad_pcb` S-expression writer, FastHenry `.inp`, openEMS script |
| Core | `litzgen.core.report` | Self-contained HTML report (inline SVG) |
| Front-end | `litzgen.kicad` | pcbnew ActionPlugin + wx dialog: generate, regenerate, simulate |
| Front-end | `litzgen.cli` | Generate / simulate / sweep / export without KiCad |

Data flow: parameters → transposition schedule → strand paths → (DRC, board items, solver mesh) → report.

## Geometry and transposition

The generator treats the bundle cross-section as an R × C matrix of slots. It peels that matrix into concentric rings and rotates every ring by one slot at each transposition step. For the paper's 4 × 4 matrix this gives the 12-strand outer bundle and the 2 × 2 inner bundle; a 6 × 6 matrix gives rings of 20, 12 and 4.

**Path model.** Each strand is a list of pieces in (angle, lateral slot, layer) space: dwells, jogs and via hops. Only the last step maps them onto the Archimedean spiral, so the KiCad writer, the DRC and the solver all read the same description.

**Why lanes are needed.** A completely full ring can't rotate one move at a time: every strand waits for the slot ahead of it to empty, so the moves deadlock. Each step therefore borrows a lane one pitch outside the ring. Ring 0 uses the gap between turns; ring k uses the column slots that ring k−1 has just vacated.

![One transposition zone, layer by layer](img/transposition-zone.png)

**Step sequence (one zone).**

1. X (exit): strands that must change layer jog sideways into the lane, while the rows shift sideways at the same time. On every layer the moving groups spread apart, so these parallel jogs never clash.
2. V (vias): strands in a lane change layer through adjacent-layer vias. The outboard lane goes bottom-first and the inboard lane top-first, so each via lands in an empty slot.
3. R (re-entry): lane strands jog back into the vacated column slots.

Ring k places 2(R−2k−1) vias per step, so the 4 × 4 matrix uses 6 + 2 = 8 vias per step. Every via spans one dielectric: blind F.Cu–In1.Cu and In2.Cu–B.Cu, buried In1.Cu–In2.Cu. A test checks that the plan is a clash-free rotation for every matrix from 2 × 2 to 8 × 8.

**Zone sizing comes from the design rules.** Parallel jogs sit p·cos(α+ψ) apart, where ψ is the spiral's own pitch angle. The jog angle α must therefore satisfy p·cos(α+ψ) ≥ (w + w_jog)/2 + clearance. Consecutive vias in a lane share a layer, and the track end caps there set the via pitch to w + clearance. Jogs are drawn as three-point arcs so the chord doesn't sag into the clearance.

| Paper coil (defaults) | Value |
| --- | --- |
| Zones | 118 (24 per turn, none within 6° of a terminal) |
| Jog / via pitch / zone length | 2.94 mm / 0.95 mm / 8.3 mm |
| Vias | 944: 236 blind F.Cu–In1.Cu, 472 buried, 236 blind In2.Cu–B.Cu |
| Strand length spread | 1859.6–1864.6 mm (0.27 %) |
| Built-in DRC | 0 violations |

The paper reports 1140 vias. It doesn't fully specify its routing, so the count here is simply what the step sequence above produces.

**Nets and DRC.** By default all strands share one net, which is electrically correct, but it means KiCad's DRC can't see a short between two strands. `litzgen.core.drc` treats every strand as its own conductor and checks copper clearance on each layer (tracks as capsules, via pads as discs) plus hole-to-hole spacing. A randomised sweep of 300 parameter sets (2–8 layers, 2–6 columns, staggered and necked jogs, hollow bundles) accepted 181 of them, and all 181 were clean. The other 119 were rejected up front with a specific reason. Turn on `per_strand_nets` if you also want KiCad's own DRC to check it.

## Simulation

The solver splits the problem the same way the physics splits. Inductance and strand current sharing come from a 3D partial-inductance model of the actual strand paths. AC resistance comes from an eddy-current solve of the actual cross-section, wherever the strands happen to be. The two meet in a solve of the network of parallel strands.

**1. 3D PEEC inductance (per strand).** Every strand is cut into straight segments of at most 3 mm, and every via becomes a vertical barrel. The paper coil gives about 14,900 segments. The Neumann integral is evaluated at three accuracy levels depending on distance:

- far pairs use the midpoint rule;
- mid-range pairs use 4-point Gauss–Legendre;
- near pairs use 8-point Gauss–Legendre averaged over 2 points across the trace width, or the reduced kernel for pieces of the same strand.

Self terms use the closed form of the reduced kernel, which reproduces Grover's formula for a rectangular bar. Only the 16 × 16 strand matrix is accumulated, so memory stays small.

**2. Axisymmetric eddy-current resistance.** At a given azimuth, every strand of every turn is treated as a rectangular ring conductor split into 10 × 3 cosine-graded filaments (the skin depth at 6.78 MHz is 25 µm). Per radian:

```math
Z_f = \operatorname{diag}\!\left(\frac{\rho r_i}{A_i}\right) + \frac{j\omega}{2\pi} M, \qquad Z_c = \left(B^{\mathsf T} Z_f^{-1} B\right)^{-1}
```

M uses the elliptic-integral formula for coaxial rings. Near pairs use the exact geometric mean distance between the two rectangles, which keeps elongated filaments accurate. Re(Z_c) covers skin effect, proximity effect and screening across all 80 strand conductors.

**3. Assembly along the path.** The azimuth is sampled every 0.5°. At each sample, the position of every strand in every turn is read from its path, lanes included. Each distinct pattern of occupied positions is solved once at 3 representative angles and interpolated in between; the paper coil has 33 such patterns. Re(Z_c) is then added into the strand resistance matrix. Its off-diagonal terms are the proximity coupling between strands, which is exactly what transposition is meant to even out.

**4. Network.** The strands are connected in parallel between the terminal bars:

```math
\left(R_{\text{eddy}} + R_{\text{via}} + j\omega L_{3D}\right) I = V\,\mathbf{1}, \qquad Z_{eq} = \frac{1}{\mathbf{1}^{\mathsf T} Z^{-1} \mathbf{1}}
```

Each via adds the resistance of its plated barrel, with a skin correction from Dowell's foil factor (0.83 mΩ per via here).

| Quality preset | Filaments per strand | Angles | Paper coil runtime (2 cores) | ESR |
| --- | --- | --- | --- | --- |
| fast | 8 × 3 | 2 | 32 s | 349 mΩ |
| standard | 10 × 3 | 3 | ≈ 4 min | 345 mΩ |
| fine | 14 × 4 | 4 | 16 min | 343 mΩ |

## Validation

The solver matches analytic references to within 0.2 % for inductance and about 3 % for eddy loss. For the paper's coil it lands within 1 % of the measured inductance and 10 % of the HFSS ESR.

### Analytic checks

| Check | Reference | Deviation |
| --- | --- | --- |
| Ring self-inductance, 180-segment polygon | µ0R[ln(8R/GMD) − 2] | −0.03 % |
| Strip-to-strip mutual inductance | width-averaged elliptic ring formula | ≤ 0.17 % (72 segments), ≤ 0.05 % (360) |
| Elliptic K, E by AGM | scipy / numeric Neumann integral | < 1e−15 / < 1e−4 |
| Self-GMD of a square / thin strip | 0.44705a / 0.2235w | exact to 5 digits |
| Eddy solver, DC limit | ρr/A | 2e−5 |
| Eddy solver, coaxial foils Δ = 1, 2, 3 | Dowell (m = 1) | +0.4 %, −2.1 %, −3.1 % (−0.8 % at Δ = 3 with a 1000:1 foil) |

### Against the paper (Table I coil, standard quality)

| Quantity at 6.78 MHz | Paper, calculated | Paper, HFSS | Paper, measured | LitzGen |
| --- | --- | --- | --- | --- |
| Inductance (µH) | 3.25 | 3.33 | 3.44 | 3.475 |
| ESR (mΩ) | – | 383 | 425 | 345 |
| Q | – | 369 | 344.5 | 429 |
| ESR, 6.4 → 7.2 MHz | – | – | 347 → 552 (+59 %) | 335 → 356 (+6 %) |

Table I in the paper doesn't add up. With d_out = 160 mm, 5 turns and 3.8 mm bundles, d_in = 69 mm implies a 4.54 mm gap between turns, while a 5.5 mm gap implies d_in = 59.4 mm. The defaults use d_in = 69 mm, which gives L = 3.48 µH, close to the measured value. The 5.5 mm gap gives 3.15 µH instead.

The measured ESR rises 59 % across a 12 % frequency band, but skin and proximity effects alone account for only about 6 %. That points to capacitive or self-resonance effects, which a quasi-static model leaves out.

### What the model says about the design

| Variant (same envelope) | Strands | Vias | L (µH) | ESR (mΩ) | Q | Inner-bundle current |
| --- | --- | --- | --- | --- | --- | --- |
| Dual-bundle, as in the paper | 16 | 944 | 3.475 | 345 | 429 | 4.0 % |
| Hollow: outer 12-strand ring only | 12 | 708 | 3.476 | 306 | 484 | – |
| Untransposed 4 × 4 | 16 | 0 | 3.474 | 304 | 487 | 2.6 % |
| Two-layer braid, 4 per layer, F.Cu/B.Cu | 8 | 236 | 3.519 | 327 | 459 | – |

Both results below are model predictions. Check them with FastHenry or a measurement before relying on them.

1. **The inner bundle carries almost no current.** With equal currents, each inner strand would link about 2 % more flux than an outer one (3.55 vs 3.48 µH). At 6.78 MHz that is roughly 3 Ω of extra reactance against about 5 Ω of strand resistance, so current crowds into the outer ring.
2. **With 0.8 mm strands, transposition does not lower ESR at 6.78 MHz.** Each strand is 32 skin depths wide, so eddy currents inside each strand dominate. The same solver does reproduce the textbook Litz gain with 0.12 × 0.018 mm strands: 405 mΩ transposed vs 503 mΩ untransposed at 6.78 MHz, and equal at 200 kHz.

## KiCad integration

| Item | How it is created | Notes |
| --- | --- | --- |
| Dwell and jog copper | `PCB_ARC` (3-point arcs, ≤ 10°) | follows the spiral, no chord sag |
| Vias | `PCB_VIA`, blind/buried, `SetLayerPair` | via-type enum resolved at run time for different KiCad versions |
| Terminal bars | footprint written as `.kicad_mod`, loaded with `FootprintLoad` | THT custom pad on all copper layers + up to 12 stitching holes |
| Parameter record | `PCB_TEXT` on Cmts.User | lets you reopen, edit and regenerate the coil |
| Group | `PCB_GROUP` named `LitzCoil` | move, delete or regenerate the coil as a unit |

**Exports.** FastHenry (`.inp`): every strand becomes a chain of rectangular bars split into nwinc × nhinc filaments. Strands are joined with `.equiv` at each terminal and driven through one `.external` port. openEMS: a bundle-level FDTD model with one conducting sheet per layer, for self-resonance and capacitance. Full strand detail isn't practical in FDTD at this frequency.

## Limitations and roadmap

| Area | Current state | Next step |
| --- | --- | --- |
| pcbnew front-end | Tested only against a stand-in for the pcbnew module | Generate the paper coil in KiCad 9 and run DRC, then again with *One net per strand* |
| KiCad IPC API (kipy) | Not written; the core doesn't depend on pcbnew | Add a thin kipy adapter for KiCad 10 |
| Parasitic capacitance / SRF | Not in the built-in solver; openEMS script provided but not run | Add a 2D electrostatic turn-to-turn model and a lumped SRF estimate |
| Eddy model | Locally axisymmetric; lanes snapped to slots; jog length handled with a correction factor | Compare against a FastHenry run of the exported `.inp` |
| Inner terminal | Terminal bar only; no lead-out (all layers are used) | Optional bridge footprint or a fifth layer |
| Inner-bundle starvation | 4 % of the current instead of 25 % | Add a step that swaps strands between rings so every strand visits every slot |

Two practical settings already exist: `steps_per_turn = -1` picks the densest transposition that fits, and `omit_rings = [1]` builds the hollow bundle.
