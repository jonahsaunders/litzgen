# Changelog

## 1.0.2 (2026-09-25)

- Fix KiCad Plugin and Content Manager rejecting the package ("Unable to parse package metadata"): the `wireless power` tag is now `wireless-power`, as PCM tags may not contain spaces

## 1.0.1 (2026-09-25)

- Fix repository links to point at jonahsaunders/litzgen

## 1.0.0 (2026-09-25)

First public release.

- Concentric multi-bundle PCB Litz coil generator (any layer and column count) using a lane-based step sequence that rotates each ring without deadlock, with adjacent-layer blind and buried vias
- Built-in strand-to-strand and hole-to-hole DRC
- Solver: 3D PEEC inductance and current sharing, axisymmetric eddy-current AC resistance, via losses, HTML report
- KiCad 8/9 action plugin (pcbnew) with a parameter dialog, regeneration from stored parameters and a background simulation
- CLI: `params`, `generate`, `simulate`, `sweep`, `export` (KiCad board, FastHenry, openEMS)
