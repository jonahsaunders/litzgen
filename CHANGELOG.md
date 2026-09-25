# Changelog

## Unreleased

- Redesigned the KiCad dialog following Apple's Human Interface Guidelines: parameters grouped into tabs with short labels, units and placeholders; one default button (Place Coil / Update Coil) on the trailing edge; confirmations before Remove Coil, Reset to Paper Defaults and placing a coil with clearance violations; field validation that points at the invalid field; plain-language status line with the activity log collapsed; simulation progress with Cancel and an Open Report button instead of opening the browser automatically; accessible names on every field; Return and Escape work
- Fix updating and removing a coil in KiCad 10: `GetParentGroup()` now returns an `EDA_GROUP` (no `m_Uuid`), and removing items with `BOARD.Remove()` corrupted pcbnew's SWIG state; items are now deleted with `BOARD.Delete()`
- Fix a crash risk when the dialog was closed while a simulation was running

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
