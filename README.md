# pyVTTrac

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://tsukada-cs.github.io/pyVTTrac/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://tsukada-cs.github.io/pyVTTrac/dev)
[![Build Status](https://github.com/tsukada-cs/pyVTTrac/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/tsukada-cs/pyVTTrac/actions/workflows/CI.yml?query=branch%3Amain)
[![PyPI](https://img.shields.io/pypi/v/pyVTTrac.svg)](https://pypi.org/project/pyVTTrac/) -->
![Python](https://img.shields.io/badge/python-3.10%2B-blue)
![Platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20macOS-lightgrey)
![License](https://img.shields.io/badge/license-BSD--2--Clause-green)

**Velocimetry by Template Tracking** — a NumPy-first Python package with a
native Fortran computational core (no Julia, no compiled-language runtime to
install).

pyVTTrac conducts simple PIV-style (particle image velocimetry) template
matching, applied repeatedly in a Lagrangian manner (as in PTV, particle
tracking velocimetry) over a configurable number of steps. It's built for
tracking coherent features — clouds, eddies, any translating pattern — in
sequences of image-like 2-D data. Both forward and backward tracking are
supported.

> **Upgrading from v1?** The Julia backend and the old `VTT`/`setup`/`trac`
> API were replaced in v2.0.0. See
> [`docs/migration-v1-to-v2.md`](docs/migration-v1-to-v2.md) — a handful of
> algorithm bugs were fixed along the way, so results can legitimately change.

## 30-second quickstart

```bash
pip install pyVTTrac
```

```python
import numpy as np
import pyvttrac as vt

# z: (nt, ny, nx) image-like data, e.g. satellite imagery, PIV frames, ...
z = np.load("images.npy").astype(np.float32)

# a regular grid of seed points to track
x0, y0 = vt.seed_grid(z.shape, spacing=8, margin=10)

res = vt.track(
    z, x0, y0, t0=0,
    template=(7, 7),              # (ny, nx) template size
    search_velocity=(2.0, 2.0),   # (vy, vx) search range, in pixels/frame
    nsteps=5,
)

print(f"tracked {res.ok.sum()} / {res.ok.size} points")
print(f"mean velocity: vx={np.nanmean(res.vx):.3f}, vy={np.nanmean(res.vy):.3f}")

# Need an xarray.Dataset or a tidy long-format table instead?
ds = res.to_xarray()       # requires: pip install "pyVTTrac[xarray]"
df = res.to_dataframe()    # requires: pip install "pyVTTrac[pandas]"
```

See [`docs/quickstart.md`](docs/quickstart.md) for a runnable, plottable
version of this example (`examples/quickstart.py`), and
[`docs/api.md`](docs/api.md) for the full parameter reference and `Status`
code table.

## How it works

- **Scoring**: template matching by sliding cross-correlation
  (`method="xcor"`, the default — `cov(x',y')/sig(x)/sig(y)`) or normalized
  covariance (`method="ncov"` — `cov(x',y')/sig(x)^2`), where `x` is the
  template and `y` is the slid target sub-image.
- **Subgrid refinement**: an optional 5-point paraboloid (`subgrid="paraboloid"`,
  the default) or Gaussian (`subgrid="gaussian"`) fit around the score peak,
  or disable it entirely (`subgrid=None`) for integer-pixel tracking.
  Non-integer seed positions are read via bilinear interpolation either way.
- **Multi-step tracking**: each seed is tracked for `nsteps` steps of `step`
  frames each (`step` can be negative, for backward tracking), carrying the
  previous step's velocity forward as the next search center.
  `vt.concat_bidirectional()` stitches a forward and a backward run from the
  same starting point into one trajectory.
- **Screening**: `min_score` (first step vs. subsequent steps), an optional
  `max_velocity_change` trajectory-consistency check, an optional
  `min_contrast` template check, and an optional `min_peak_prominence`
  interior-peak check.
- **Missing data**: `mask` (boolean, `True` = ignore) and/or `missing_value`
  (a sentinel in `z`); `NaN` in `z` is always treated as missing as well.
- **Coordinates**: positions/velocities are index-based (grid spacing = 1) by
  default; pass a `pyvttrac.Grid` to work in physical units instead.

## Performance

The Fortran core is called once per `track()` invocation (not once per
tracking step), and OpenMP-parallelizes across seed points. On a benchmark (400×400×12 field, 7396 seeds, 5 steps),
sequential tracking is faster than the Julia backend it replaces, and scales
close to linearly with the `workers` parameter:

```text
workers=1   0.145 s   (Julia v2.0.0 backend: 0.244 s)
workers=8   0.021 s
workers=-1  0.014 s   (all cores)
```

Reproduce with `python tools/bench.py`.

## Installation

```bash
pip install pyVTTrac
```

Requires Python 3.10+ on **Linux or macOS** (Windows is not supported).
Optional extras:

```bash
pip install "pyVTTrac[xarray]"   # for TrackResult.to_xarray()
pip install "pyVTTrac[pandas]"   # for TrackResult.to_dataframe()
```

Building from source needs a Fortran compiler (`gfortran` or equivalent);
OpenMP is used automatically if available, but the package builds and runs
fine without it. On macOS: `brew install gcc`.

### Running the tests

```bash
git clone https://github.com/tsukada-cs/pyVTTrac.git
cd pyVTTrac
pip install -e ".[test]"
pytest
```

## Related packages

- [`VTTrac.jl`](https://github.com/tsukada-cs/VTTrac.jl) by Taiga Tsukada —
  the Julia implementation this package was originally built on top of (up
  to v1), and whose v2.0.0 algorithm this package's Fortran core is a direct
  port of. It also serves as this package's reference implementation for
  golden-data testing (`tools/gen_golden.jl`).
- [`VTTrac`](https://github.com/thorinouchi/VTTrac) by Takeshi Horinouchi —
  the original implementation `VTTrac.jl` (and, transitively, pyVTTrac) is
  based on.

## References

- Horinouchi, T., S. Tsujino, M. Hayashi, U. Shimada, W. Yanase, A. Wada, and
  H. Yamada, 2023: Stationary and Transient Asymmetric Features in Tropical
  Cyclone Eye with Wavenumber-1 Instability: Case Study for Typhoon Haishen
  (2020) with Atmospheric Motion Vectors from 30-Second Imaging. *Monthly
  Weather Review*, 151, 253–273, <https://doi.org/10.1175/MWR-D-22-0179.1>.
- Tsukada, T., T. Horinouchi, and S. Tsujino, 2024: Wind Distribution in the
  Eye of Tropical Cyclone Revealed by a Novel Atmospheric Motion Vector
  Derivation. *JGR Atmospheres*, 129, e2023JD040585,
  <https://doi.org/10.1029/2023JD040585>.

## License

BSD 2-Clause License — see [`LICENSE`](LICENSE).
