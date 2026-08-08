# pyVTTrac documentation

pyVTTrac tracks coherent features (clouds, eddies, any translating pattern)
through sequences of image-like data by repeated PIV-style template
matching, in a Lagrangian (PTV-style) manner. See the top-level
[`README.md`](../README.md) for a project overview, installation, and a
30-second example.

- **[Quickstart](quickstart.md)** — a runnable, plottable example
  (`examples/quickstart.py`) walked through step by step.
- **[API reference](api.md)** — every parameter of `track()`/`Tracker`,
  `TrackResult`'s fields and methods, `Grid`, `seed_grid()`, and the
  `Status` code table.
- **[Migrating from v1](migration-v1-to-v2.md)** — parameter mapping,
  before/after examples, and the numerical differences from v1 you should
  know about before comparing results across versions.

## Provenance

pyVTTrac's Fortran computational core is a direct port of
[`VTTrac.jl`](https://github.com/tsukada-cs/VTTrac.jl) v2.0.0, which is
itself based on [`VTTrac`](https://github.com/thorinouchi/VTTrac) by Takeshi
Horinouchi. `VTTrac.jl` also serves as the reference implementation this
package's test suite is validated against (`tests/test_golden.py`,
`tools/gen_golden.jl`).
