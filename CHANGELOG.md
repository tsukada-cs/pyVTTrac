# Changelog

## 2.1.0

Driven by openTCAMV's migration to the v2 API (see openTCAMV's own
CHANGELOG for the consumer-side details).

### Fixed

- **`Grid.dx`/`Grid.dy` sign**: previously silently coerced to `abs(dx)`/
  `abs(dy)`, which flipped `to_index_x`/`to_index_y` for a descending
  coordinate axis (e.g. latitude/`y` in many satellite products —
  `Grid.from_coords([10, 9, 8, 7], ...).to_index_x(9.0)` returned `-1.0`
  instead of `+1.0`). `dx`/`dy` are now kept signed; only zero is rejected.
  **If you were relying on the old silent-`abs()` behavior with an
  explicitly negative `dx`/`dy`, your index/physical mapping will now come
  out with the opposite sign** (the previously-correct case — positive
  `dx`/`dy`, or a `Grid` built via `from_coords`/`grid="auto"` from an
  ascending coordinate — is unaffected).
- **`max_velocity_change` silently disabled by a negative `Grid.dy`/`dx`**:
  the physical-to-index conversion of a positive `max_velocity_change` came
  out negative whenever `Grid.dy`/`dx` was negative, and the `dvy > 0 and
  dvx > 0` gate then silently treated the check as unset. Now compared by
  magnitude (`abs()`), as intended.

### Added

- `Tracker.track(..., step=...)`: per-call override of the `Tracker`'s
  configured `step`, without needing a second `Tracker` for the opposite
  direction (e.g. one `Tracker` for both forward and backward tracking).
- `search_radius_from_velocity(search_velocity, dt, grid=None)`: the
  `ceil(abs(v_index * dt)) + 1` pixel-radius formula `track()` uses
  internally when `search_radius` isn't given, exposed for callers that need
  it ahead of time (e.g. to size a buffer, or to derive a radius from a
  reference `dt` other than the data's own mean spacing).
- `diagnostics` now also accepts `"templates"` or `"score_grids"` (or a
  tuple of the two) to populate only one of `TrackResult.templates` /
  `.score_grids`, avoiding the other (often large) allocation. `True`/
  `False` behave as before (both / neither).
- `TrackResult.step`: the `step` a result was tracked with. `to_xarray()`'s
  `step`/`step_v` coordinates are now scaled (and sign-flipped for backward
  tracking) by it, instead of always assuming `step=1` — a backward-tracked
  result's `to_xarray()` axes now carry the correct sign/magnitude.

## 2.0.0

### Breaking changes

- **Backend rewrite**: the Julia backend (`juliacall` + the `VTTrac.jl` git
  submodule) has been replaced with a native Fortran computational core
  (built via meson-python + f2py). Julia is no longer a runtime or build
  dependency of pyVTTrac; it remains only as the reference implementation
  used to generate this package's golden test data (`tools/gen_golden.jl`,
  developer-only, not run in CI).
- **API rewrite**: the `VTTrac.VTT` / `.setup()` / `.trac()` object API is
  removed. Use `pyvttrac.track()` (one call) or `pyvttrac.Tracker` (reusable
  config). Import name changes from `pyVTTrac` (mixed case) to `pyvttrac`
  (lowercase); the PyPI distribution name is unchanged (`pyVTTrac`).
- **`(x, y)` axis order**: `nsx`/`nsy`, `vxhw`/`vyhw`, `ixhw`/`iyhw`, and
  `vxg0`/`vyg0` are replaced by `(ny, nx)`-ordered tuples: `template`,
  `search_velocity`, `search_radius`, `first_guess`.
- **Missing-value convention**: invalid float outputs are now `NaN` (was
  `-999.0`); invalid `t_index` entries are `-1` (was `-999`). `NaN` in `z` is
  now always treated as missing data, in addition to (and equivalent to) an
  explicit `missing_value`.
- **Return type**: `track()` returns a `TrackResult` dataclass (attribute
  access) instead of a 10-element tuple. `to_xarray()` produces the closest
  equivalent to v1's `asxarray=True` `Dataset`; `xarray` is now an optional
  dependency (`pip install "pyVTTrac[xarray]"`).
- **Windows is not supported.** Linux and macOS only.
- **Python >= 3.10** is now required.

See `docs/migration-v1-to-v2.md` for the full parameter mapping and
before/after examples.

### New

- `pyvttrac.Grid`: unifies v1's `set_grid_par()` / `setup_eq_grid()` /
  `trac_eq_grid()` / `calc_ixyhw_from_v_eq_grid()` into a single object for
  physical-coordinate tracking.
- `pyvttrac.seed_grid()`: convenience helper to generate a regular grid of
  seed points.
- `TrackResult.to_dataframe()`: long-format table (one row per seed × step),
  convenient for AMV-style analysis. Requires the optional `pandas`
  dependency (`pip install "pyVTTrac[pandas]"`).
- `workers` parameter: controls the number of OpenMP threads parallelizing
  across seed points (`None` = OpenMP default, `-1` = all cores, or a
  specific count).
- Explicit input validation with clear exceptions in place of v1 inputs that
  used to crash or silently misbehave.
- First-call latency drops from ~12 s (Julia cold start) to effectively zero,
  since there is no longer a Julia runtime to initialize. Sequential
  tracking throughput is also faster than the Julia backend it replaces
  (measured ~40% faster under the benchmark), and
  scales close to linearly with `workers` via OpenMP.

### Fixed (inherited from VTTrac.jl v2.0.0 — results can differ from v1)

1. `method="ncov"` with a mask: normalization corrected from `cov/std(x)` to
   `cov/var(x)`, matching the unmasked case. Scores, and therefore
   `min_score` accept/reject decisions, can change.
2. `max_velocity_change` rejection: v1 incorrectly invalidated the *initial*
   tracking point on rejection; only the derived point is invalidated now.
3. `diagnostics` (v1 `out_subimage`): the final sub-image is now read from
   the correct final time index (was read from the wrong one in v1).
4. `xcor`/`ncov` variance computation: numerically stabilized. v1's
   incremental variance update could go spuriously negative for data with a
   large mean offset relative to its variance, silently corrupting
   subsequent scores (measured up to 13% error in that scenario).
5. `fixed_template` (v1 `use_init_temp`) combined with a mask: raised a
   runtime error in v1; works correctly now.

### Removed

- `VTTrac.VTT`, `.setup()`, `.trac()`, `.set_grid_par()`, `.setup_eq_grid()`,
  `.trac_eq_grid()`, `.calc_ixyhw_from_v()`, `.calc_ixyhw_from_v_eq_grid()`.
- The `VTTrac.jl` git submodule and the `juliacall` runtime dependency.
- `setup.py` / `setup.cfg` (replaced by `pyproject.toml` + meson-python).

---

## 1.0.0 and earlier

See git history prior to the v2.0.0 rewrite.
