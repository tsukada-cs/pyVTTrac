# API reference

For a worked example, see [`docs/quickstart.md`](quickstart.md). For the
mapping from the removed v1 API, see
[`docs/migration-v1-to-v2.md`](migration-v1-to-v2.md).

## `pyvttrac.track(z, x0, y0, t0=0, *, ...)`

Track template patches through a `(nt, ny, nx)` image sequence. Returns a
[`TrackResult`](#trackresult).

### Data arguments

| Parameter | Type | Default | Description |
|---|---|---|---|
| `z` | array-like, `(nt, ny, nx)` | required | Image-like data; at least 2 time steps. Converted to `float32`. |
| `x0`, `y0` | array-like, any shared shape | required | Seed positions, index units (or physical units if `grid` is given). Non-integer values are read via bilinear interpolation. |
| `t0` | scalar or array matching `x0`/`y0` | `0` | Starting time index for each seed. |
| `time` | array-like, `(nt,)` | `None` (→ `0, 1, 2, ...`) | Physical/absolute time coordinate for `z`'s first axis; used for `search_velocity`/`max_velocity_change` unit conversion and as the `dt` between steps. |
| `mask` | array-like bool, same shape as `z` | `None` | `True` = ignore this pixel when scoring. A mask that is all-`False` is treated the same as no mask. |
| `first_guess` | `(vy0, vx0)`, each broadcastable to `x0`'s shape | `None` (→ zero) | Initial velocity guess to search around, for the first step only. |
| `missing_value` | scalar | `None` | A sentinel value in `z` treated as missing. `NaN` in `z` is *always* treated as missing as well, regardless of this setting. |
| `grid` | `Grid`, `"auto"`, or `None` | `None` | If given, `x0`/`y0`/velocities are interpreted (and returned) in physical units — see [`Grid`](#grid). `"auto"` infers a `Grid` from an `xarray.DataArray`'s coordinates (requires `xarray`; coordinates must be equally spaced). |
| `diagnostics` | `bool`, `"templates"`, `"score_grids"`, or a tuple of the latter two | `False` | Which of `TrackResult.templates`/`.score_grids` to populate. `False` = neither (default), `True` = both. A single name or `(name, ...)` populates only those, avoiding the (often large) allocation for the side you don't need. |

### Tracking-configuration arguments

These are the same ones `Tracker(...)` accepts, for reuse across multiple
`track()` calls (see [`Tracker`](#tracker)).

| Parameter | Type | Default | Description |
|---|---|---|---|
| `template` | `(ny, nx)` | required | Template sub-image size. |
| `search_velocity` | `(vy, vx)` | `None` | Search range as a velocity (index units, or physical units under `grid`); the actual search radius is derived from this and `time`. Mutually exclusive with `search_radius`; exactly one is required. |
| `search_radius` | `(iy, ix)`, pixels | `None` | Search range as a pixel radius directly. Mutually exclusive with `search_velocity`. |
| `nsteps` | int | `2` | Number of tracking steps per seed. |
| `method` | `"xcor"` \| `"ncov"` | `"xcor"` | Scoring: cross-correlation coefficient, or covariance normalized by the template's own variance. |
| `subgrid` | `"paraboloid"` \| `"gaussian"` \| `None` | `"paraboloid"` | Subgrid peak refinement method, or `None` for integer-pixel tracking. |
| `step` | int | `1` | Index stride between tracked frames; negative = backward tracking. Must be nonzero. On `Tracker.track()`, also accepted as a per-call keyword (`tracker.track(z, x0, y0, step=-1)`) that overrides the `Tracker`'s configured `step` for that call only, without mutating `tracker.config` — useful for running forward and backward tracking from one `Tracker`. |
| `min_score` | float or `(first, subsequent)` | `(0.8, 0.7)` | Minimum score to accept a step; a scalar applies to both. |
| `max_velocity_change` | `(dvy, dvx)` | `None` | Reject a step whose velocity differs from the previous step's by more than this, in *both* components. Only active when both components are given and positive. |
| `min_contrast` | float | `None` | Reject a template whose `max - min` is below this. |
| `min_peak_prominence` | float | `None` | Reject a template with no sufficiently prominent interior peak/trough (needs `template` >= 3 in both dimensions). |
| `fixed_template` | bool | `False` | Keep matching against the *initial* template throughout the trajectory, instead of re-reading a fresh template at each step. |
| `min_samples` | int | `1` | Minimum number of jointly-visible (unmasked) pixels required to compute a score, when `mask` is given. |
| `workers` | int | `None` | OpenMP threads used to parallelize across seeds: `None` (OpenMP default / `OMP_NUM_THREADS`), `-1` (all cores), or a specific count. |

## `Tracker`

```python
tracker = vt.Tracker(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=5, ...)
res_a = tracker.track(z_a, x0, y0)
res_b = tracker.track(z_b, x0, y0)
```

Bundles the tracking-configuration arguments listed above (everything
`track()` accepts *except* the data arguments); `Tracker.track()` accepts
the same data arguments as `track()` (`z`, `x0`, `y0`, `t0`, `time`, `mask`,
`first_guess`, `missing_value`, `grid`, `diagnostics`). `tracker.config` is
the underlying `TrackingConfig`, if you need to inspect it.

## `TrackResult`

Returned by `track()` / `Tracker.track()`. Let `seed_shape` be the
broadcast shape of `x0`/`y0`.

| Field | Shape | dtype | Description |
|---|---|---|---|
| `count` | `seed_shape` | int | Number of valid trajectory points, including the initial one (`0` .. `nsteps+1`). |
| `status` | `seed_shape` | int (`Status` values) | Why tracking stopped for this seed (`Status.OK` if it completed all `nsteps`). |
| `t_index` | `(nsteps+1, *seed_shape)` | int | Time index at each point; `-1` once invalid. |
| `x`, `y` | `(nsteps+1, *seed_shape)` | float | Tracked positions; `NaN` once invalid. |
| `vx`, `vy` | `(nsteps, *seed_shape)` | float | Velocity between consecutive points; `NaN` where invalid. |
| `score` | `(nsteps, *seed_shape)` | float | Match score at each step; `NaN` where invalid. |
| `templates` | `(nsteps+1, ny, nx, *seed_shape)` or `None` | float | Sub-images along the trajectory, if requested via `diagnostics`. |
| `score_grids` | `(nsteps, 2·iy+1, 2·ix+1, *seed_shape)` or `None` | float | Full score arrays at each step, if requested via `diagnostics`. |
| `step` | scalar | int | The `step` this result was tracked with (config's or the per-call override). Determines the `step`/`step_v` coordinates in `to_xarray()`. |
| `grid` | scalar | `Grid` or `None` | The `Grid` used, if any (`None` if tracking was done in index units). |
| `seed_x`, `seed_y` | `seed_shape` | float or `None` | The seed positions as passed to `track()`, in whatever units `x0`/`y0` were given (physical units when `grid` is set). |
| `search_radius` | scalar | `(iy, ix)` or `None` | The pixel search radius actually used (whether given directly or derived from `search_velocity`). |
| `status_forward`, `status_backward` | `seed_shape` | int or `None` | Set only on a `concat_bidirectional()` result: each leg's own `status`. `None` on an ordinary result. |
| `step_offset` | scalar | int | Index, along the position axis, of the origin (step 0). `0` for an ordinary result; set by `concat_bidirectional()`. Feeds into `to_xarray()`'s `step`/`step_v` coordinates -- you shouldn't need to touch this directly. |

Properties/methods:

- `.ok` — bool array, `seed_shape`: `status == Status.OK`.
- `.to_xarray(*, units=None, global_attrs=None)` — an `xarray.Dataset` with
  labeled dimensions/coordinates and **CF-1.13 / ACDD-1.3 metadata**:
  a `Conventions` global attribute, `long_name`/`units` on every variable
  (`"1"`, i.e. dimensionless, by default), `status`'s CF `flag_values`/
  `flag_meanings` (derived from the `Status` enum), and `_FillValue` set via
  each variable's `.encoding` (not `.attrs`, which `Dataset.to_netcdf()`
  rejects) so the `NaN`/`-1` sentinels round-trip through netCDF correctly.
  `title`/`summary`/`keywords`/`source`/`history`/`date_created` are filled
  with generic-but-accurate defaults. If `self.grid` is set and has
  `x_units`/`y_units`/`velocity_units`, those become the default `units` for
  `x`/`y` and `vx`/`vy` respectively (an explicit `units=` entry still wins).
  If `self.seed_x`/`self.seed_y` are set, they're added as `seed_x`/`seed_y`
  coordinates -- a 1-D axis per dimension (`("seed_1",)`/`("seed_0",)`) when
  they come from a separable `seed_grid()`-style grid, otherwise a single
  auxiliary coordinate spanning all seed dimensions.
  What it still *can't* know — deployment-specific ACDD discovery attributes
  (`institution`, `creator_name`, `license`, `id`, ...) — you supply via
  `global_attrs={"institution": ..., ...}`; both `units=` and `global_attrs=`
  are merged over (and can override) the auto-generated defaults.
  The `step`/`step_v` coordinates are scaled by `self.step` and shifted by
  `self.step_offset` (`0` except on a `concat_bidirectional()` result):
  `step = (arange(nsteps+1) - step_offset) * step`,
  `step_v = (arange(nsteps) - step_offset) * step + 0.5 * sign(step)` --
  for an ordinary (non-combined) result this is unchanged from before
  (`step_offset` is always `0`), including the sign-flip for backward
  tracking. Requires the optional `xarray` dependency; raises a clear
  `ImportError` (mentioning `pip install "pyVTTrac[xarray]"`) if it isn't
  installed.
- `.to_dataframe()` — a long-format `pandas.DataFrame`, one row per
  `(seed, step)`, with columns `seed, step, count, status, t_index, x, y,
  vx, vy, score` (`vx`/`vy`/`score` are `NaN` on the `step == 0` row, since
  no step has been taken yet). Requires the optional `pandas` dependency;
  raises a clear `ImportError` if it isn't installed.

## `Grid`

```python
grid = vt.Grid(x0=0.0, y0=0.0, dx=2.0, dy=2.0, unit_factor=1.0)
```

Maps between physical and index coordinates, for `track(..., grid=grid)`:

| Field | Description |
|---|---|
| `x0`, `y0` | Physical-coordinate origin (where index 0 sits). |
| `dx`, `dy` | Physical spacing per index step. **Signed** — kept as given (only zero is rejected). Use a negative `dy` (or `dx`) for a descending coordinate axis (e.g. latitude in many satellite products); `to_index_*`/`to_phys_*` round-trip correctly either way. Code that needs a magnitude (e.g. a pixel search radius) must take `abs()` itself. |
| `unit_factor` | Unit-conversion factor applied to velocities only (e.g. seconds-per-timestep, if `dx`/`dy` are in meters and you want velocities in m/s). Default `1.0`. |
| `x_units`, `y_units`, `velocity_units` | Optional unit strings (e.g. `"m"`, `"degrees_east"`, `"m s-1"`). Pure metadata — never affect any conversion or validation here; used only as `TrackResult.to_xarray()`'s default `units` for `x`/`y` and `vx`/`vy`. |

Conversions: `index = (phys - origin) / spacing`;
`index_velocity = phys_velocity / (spacing * unit_factor)` (and their
inverses). When `grid` is given, `x0`/`y0`/`search_velocity`/
`max_velocity_change`/`first_guess` are interpreted in physical units on the
way in, and `x`/`y`/`vx`/`vy` come back in physical units.

`Grid.from_coords(x_coord, y_coord, *, x_units=None, y_units=None,
velocity_units=None)` infers a `Grid` from two 1-D, equally-spaced
coordinate arrays (raises `ValueError` if they aren't equally spaced),
optionally tagged with the same unit metadata as the `Grid` constructor;
this is what `grid="auto"` uses under the hood for an
`xarray.DataArray`.

## `search_radius_from_velocity(search_velocity, dt, grid=None)`

```python
iy, ix = vt.search_radius_from_velocity((2.0, 2.0), dt=1.5)
iy, ix = vt.search_radius_from_velocity((50.0, 50.0), dt=1.5, grid=grid)  # physical units
```

Pixel search radius `(iy, ix)` from a `(vy, vx)` search velocity and a
reference time interval `dt`, using the same formula `track()` uses
internally when `search_radius` isn't given explicitly:
`ceil(abs(v_index * dt)) + 1`. Pass `grid` if `search_velocity` is in
physical units; omit it for index units. Useful when you need the derived
search radius ahead of time (e.g. to size a diagnostics buffer, or to pick a
`dt` other than the data's own mean spacing).

## `velocity_from_search_radius(search_radius, dt, grid=None)`

```python
vy, vx = vt.velocity_from_search_radius((6, 9), dt=1.5)
```

The (approximate) inverse of `search_radius_from_velocity()`: the `(vy, vx)`
velocity magnitude a given pixel search radius `(iy, ix)` is guaranteed to
cover, for a reference time interval `dt` (`(radius - 1) / abs(dt)`, the
same convention v1 pyVTTrac used in `set_ixyhw_directly`). Not an exact
inverse, since `search_radius_from_velocity()`'s `ceil()` loses information
— the round trip `search_radius_from_velocity` → `velocity_from_search_radius`
returns a value `>=` the original. Always non-negative, and unaffected by
`dt`'s sign: a search radius alone carries no direction. Pass `grid` to get
the result in physical units. Useful for recovering an equivalent search
velocity to log/report (v1 exposed `vxhw`/`vyhw` directly; v2's
`search_radius`-first design doesn't retain them).

## `concat_bidirectional(forward, backward, *, drop_shared_origin=True)`

```python
tracker = vt.Tracker(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=5)
fwd = tracker.track(z, x0, y0, t0=10, step=1)
bwd = tracker.track(z, x0, y0, t0=10, step=-1)
combined = vt.concat_bidirectional(fwd, bwd)
```

Combines a forward- and backward-tracked `TrackResult` from a **shared
origin** (same `t0`/`x0`/`y0`) into one `TrackResult` spanning both
directions. Only this shared-origin case is supported — `forward.step > 0`,
`backward.step < 0`, equal magnitude, matching seed shapes, and
`forward.t_index[0] == backward.t_index[0]`; anything else raises
`ValueError`. Seeds whose forward/backward legs start at different times
must be combined by the caller.

The combined position axis (`t_index`/`x`/`y`/`templates`) runs from
`-len(backward.vx)` to `+len(forward.vx)` steps (`to_xarray()`'s `step`/
`step_v` coordinates reflect this automatically, via `step_offset`).
Velocities are **not** sign-flipped when reversing the backward leg — `vx =
(xw - xcur) / dt` already comes out correctly signed when `dt < 0`, so only
the array order changes.

- `status` is `Status.OK` only if both legs completed OK; otherwise
  whichever leg's status is not OK (`forward`'s, if both failed). Each
  leg's own status survives separately as `status_forward`/
  `status_backward`.
- `count` is the number of valid (`t_index != -1`) points on the combined
  position axis.
- `grid`/`seed_x`/`seed_y`/`search_radius` carry over only when `forward`
  and `backward` agree; otherwise `None`.
- `diagnostics` (`templates`/`score_grids`) must be present on both legs or
  neither, with matching non-step-axis shapes.
- `drop_shared_origin=True` (default) records the shared starting point
  once; `False` keeps both legs' own copy of it (one array element longer).

## `seed_grid(shape, spacing, margin=0)`

```python
x0, y0 = vt.seed_grid(z.shape, spacing=8, margin=10)
```

A regular grid of seed points over the last two dimensions of `shape`
(i.e. `(ny, nx)`), in index units. `spacing`/`margin` are a scalar (applied
to both axes) or `(y, x)` pair. Returns 2-D `x0`, `y0` arrays.

## `Status`

`IntEnum` describing why a seed's tracking stopped (`status` field of
`TrackResult`, and the return value of `count` reflects how far it got).
Values match `VTTrac.jl`'s `STATUS_*` constants (and v1 pyVTTrac's raw
status integers), so numeric comparisons keep working across versions.

| Value | Name | Meaning |
|---|---|---|
| 0 | `OK` | Completed successfully (or the check hasn't failed yet). |
| 1 | `TID_START_OUT_OF_RANGE` | This step's starting time index is out of range. |
| 2 | `TEMPLATE_READ_FAILED` | Template sub-image couldn't be read (out of bounds, missing data, or — with a mask — too few visible pixels per `min_samples`). |
| 3 | `LOW_CONTRAST` | Template failed the `min_contrast` check. |
| 4 | `PEAK_NOT_INSIDE_TEMPLATE` | Template failed the `min_peak_prominence` check. |
| 5 | `TID_END_OUT_OF_RANGE` | This step's ending time index is out of range. |
| 6 | `SCORE_COMPUTATION_FAILED` | The sliding score couldn't be computed (search window out of bounds, missing data, or — with a mask — too few jointly-visible pixels everywhere). |
| 7 | `PEAK_NOT_FOUND` | No valid score peak, or the peak sits on the edge of the search window. |
| 8 | `SCORE_BELOW_THRESHOLD` | Peak score below `min_score` (first component for step 1, second for later steps). |
| 9 | `VELOCITY_CHANGE_TOO_LARGE` | Velocity change along the trajectory exceeded `max_velocity_change`. |

## Input validation

`track()`/`Tracker()`/`Grid()` raise `ValueError` (rather than crashing or
silently misbehaving) for, among others:

- `z` not 3-dimensional, or fewer than 2 time steps.
- `mask` shape not matching `z`.
- `template` components < 1, or < 3 when `min_peak_prominence` is set.
- `method` not `"xcor"`/`"ncov"`; `subgrid` not a valid choice.
- Neither, or both, of `search_velocity`/`search_radius` given.
- `x0`/`y0` shape mismatch; `t0` or `time` shape mismatch.
- `step == 0`; `nsteps < 1`; `min_samples < 1`.
- `max_velocity_change` components <= 0.
- `Grid.dx`/`Grid.dy` == 0.
