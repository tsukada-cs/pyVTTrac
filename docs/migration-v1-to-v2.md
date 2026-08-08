# Migrating from pyVTTrac v1 to v2

pyVTTrac v2.0.0 is a complete rewrite. The Julia-based backend (`juliacall` +
the `VTTrac.jl` git submodule) has been replaced with a native Fortran
computational core, and the old `VTT` / `setup` / `trac` object API has been
replaced with a single function, `pyvttrac.track()`. **The v1 API has been
removed entirely** — there is no compatibility shim — so this guide exists to
make the switch mechanical.

If you only read one section, read [Numerical differences from
v1](#numerical-differences-from-v1---read-this-first): a handful of v1 bugs
were fixed in the reference algorithm this port is based on
(`VTTrac.jl` v2.0.0), so **results can legitimately change** even though
your code otherwise migrates one-to-one.

## Before / after

```python
# --- v1 (Julia-backed) ---
from pyVTTrac import VTTrac
vtt = VTTrac.VTT(z, t)
vtt.setup(nsx=5, nsy=5, vxhw=1.8, vyhw=1.8, ntrac=19)
ds = vtt.trac(tid0, x0, y0, out_subimage=True, out_score_ary=True)
print(ds.vx.mean())

# --- v2 ---
import numpy as np
import pyvttrac as vt
res = vt.track(z, x0, y0, t0=tid0, time=t,
               template=(5, 5), search_velocity=(1.8, 1.8),
               nsteps=19, diagnostics=True)
print(np.nanmean(res.vx))
# If you still want an xarray.Dataset (requires the optional xarray extra):
ds = res.to_xarray()
```

## Installing v2

- Julia, `juliacall`, and the `VTTrac.jl` git submodule are **no longer
  needed**. `git clone --recurse-submodules` is no longer required either.
- `pip install pyVTTrac` is now enough on its own.
- First-call latency drops from **~12 seconds (Julia cold start) to
  effectively zero**.
- Supported platforms are **Linux and macOS**; Windows is not supported (v1
  wasn't tested on Windows either, since it depended on `juliacall`).
- If you use `TrackResult.to_xarray()`, install the `xarray` extra:
  `pip install "pyVTTrac[xarray]"`. `to_dataframe()` (new in v2) needs
  `pip install "pyVTTrac[pandas]"`.

## Parameter mapping (v1 → v2)

| v1 | v2 | Notes |
|---|---|---|
| `from pyVTTrac import VTTrac` | `import pyvttrac as vt` | Import name is now lowercase. Distribution name on PyPI is unchanged (`pyVTTrac`). |
| `VTT()` + `setup()` + `trac()` | `vt.track()` (one call) | Reusable config: `vt.Tracker(...)`, see below. |
| `nsx`, `nsy` | `template=(ny, nx)` | **Axis order flipped to (y, x)**, consistent with `z`'s `(nt, ny, nx)` shape and NumPy conventions. |
| `vxhw`, `vyhw` | `search_velocity=(vy, vx)` | Order flipped, as above. Mutually exclusive with `search_radius`; exactly one is required. |
| `ixhw`, `iyhw` | `search_radius=(iy, ix)` | Order flipped, as above. |
| `ntrac` | `nsteps` | |
| `itstep` | `step` | Negative = backward tracking, same as v1. |
| `Sth0`, `Sth1` | `min_score=(first, subsequent)` | A scalar applies to both. |
| `vxch`, `vych` | `max_velocity_change=(dvy, dvx)` | Only active when both components are given (positive), same as v1. |
| `Cth` | `min_contrast` | |
| `peak_inside_th` | `min_peak_prominence` | |
| `use_init_temp` | `fixed_template` | |
| `score_method` | `method` | Same values: `"xcor"`, `"ncov"`. |
| `subgrid`, `subgrid_gaus` (2 bools) | `subgrid="paraboloid"` / `"gaussian"` / `None` | The two v1 booleans collapse into one parameter. |
| `zmiss` | `missing_value` | NaN in `z` is *also* always treated as missing now, regardless of `missing_value` (new in v2). |
| `fmiss`, `imiss` | *(removed)* | v1 let you pick the missing-value sentinels (default `-999.0`/`-999`); v2 always uses `NaN`/`-1` (see below) and does not expose a choice. |
| `out_subimage`, `out_score_ary` | `diagnostics=True` | Both diagnostic outputs are now controlled by one flag: `TrackResult.templates` / `.score_grids`. |
| `vxg0`, `vyg0` | `first_guess=(vy0, vx0)` | Order flipped, as above. |
| `asxarray=True` | `res.to_xarray()` | Default output is now plain NumPy; xarray is opt-in and optional. |
| `set_grid_par()` + `setup_eq_grid()` + `trac_eq_grid()` | `vt.Grid(...)` passed to `track(..., grid=...)` | Three v1 methods collapse into one object; see [Grid](#physical-coordinates-grid) below. |
| `ucfact` | `Grid(unit_factor=...)` | |
| `calc_ixyhw_from_v()` | internal | Handled automatically when `search_velocity` is given. |
| Return: 10-element tuple | Return: `TrackResult` | Attribute access instead of tuple unpacking (`res.vx`, not `result[5]`). |
| Missing value `-999.0` / `-999` | `NaN` (float outputs), `-1` (`t_index`) | See below. `count`/`status` are never negative, so they need no sentinel. |
| `status` raw `int` | `pyvttrac.Status` (`IntEnum`) | **Values are unchanged**, so numeric comparisons (`status == 3`) keep working as-is. |
| `vtt.attrs["nt"]`, etc. | Read directly from your own arrays (`z.shape`, ...) | There is no stateful tracker object to hold attrs; `Tracker.config` exposes the tracking parameters if you need them. |

### Reusable parameters: `Tracker`

If you called `setup()` once and `trac()` many times with different data,
use `Tracker` instead of repeating keyword arguments on every `track()` call:

```python
tracker = vt.Tracker(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=5)
res_a = tracker.track(z_a, x0, y0)
res_b = tracker.track(z_b, x0, y0)
```

### Physical coordinates: `Grid`

v1's `set_grid_par()` + `setup_eq_grid()` + `trac_eq_grid()` combo is now a
single `Grid` object passed to `track(..., grid=...)`:

```python
grid = vt.Grid(x0=0.0, y0=0.0, dx=2.0, dy=2.0, unit_factor=1000.0)
res = vt.track(z, x_phys, y_phys, grid=grid,
               template=(7, 7), search_velocity=(50.0, 50.0), nsteps=5)
# x0/y0/search_velocity/max_velocity_change/first_guess are given in physical
# units; x/y/vx/vy come back in physical units too.
```

The conversion formulas are unchanged from v1's `trac_eq_grid` /
`setup_eq_grid`:
- index = (phys − origin) / spacing
- index velocity = phys velocity / (spacing × `unit_factor`)

## Numerical differences from v1 — read this first

v2.0.0 is a port of `VTTrac.jl` v2.0.0, which fixed several bugs present when
v1 pyVTTrac was calling the older `VTTrac.jl`. **The same input can produce
different output between v1 and v2** for these reasons:

1. **`method="ncov"` + a mask**: the normalization was `cov/std(x)`; it is
   now correctly `cov/var(x)`, matching the unmasked case. Scores change, so
   `min_score` (`Sth0`/`Sth1`) thresholding can now accept or reject
   differently than it used to.
2. **`max_velocity_change` rejection**: when the 2nd tracking step is
   rejected, v1 incorrectly invalidated the *initial* point; v2 correctly
   invalidates only the derived (step-1) result and keeps the initial point.
3. **`diagnostics` (v1 `out_subimage`)**: the final sub-image in the
   diagnostic output is now read from the correct final time index (v1 read
   it from the wrong one).
4. **`xcor` variance computation**: stabilized numerically. For data with a
   large mean offset relative to its variance, v1's incremental variance
   update could go spuriously negative and silently corrupt an entire
   subsequent row of scores — measured at up to **13% error** in that
   scenario. v2 recomputes the (demeaned) variance directly at each window
   instead.
5. **`fixed_template` (v1 `use_init_temp`) + a mask**: this combination
   raised a runtime error in v1. It works correctly in v2.
6. Various invalid inputs that used to crash or silently misbehave in v1 now
   raise a clear exception instead (see `docs/api.md` for the full
   validation list).

**If your results changed after migrating, it is very likely v1's bug being
fixed, not a new bug in v2.** If you depend on any of the above for
reproducing a specific historical result, pin `pyVTTrac<2` rather than
patching around it.

## Other things worth knowing

- `TrackResult.to_dataframe()` is new in v2: a long-format table (one row per
  seed × step) convenient for AMV-style analysis, without needing xarray.
- Diagnostics arrays (`templates`, `score_grids`) are only computed when you
  ask for them (`diagnostics=True`); they cost time and memory, so leave them
  off unless you need to inspect intermediate sub-images/score maps.
- `workers` (new in v2) controls the number of OpenMP threads used to
  parallelize across seed points: `None` (default, respects `OMP_NUM_THREADS`
  / OpenMP's default), a specific count, or `-1` for all cores.
