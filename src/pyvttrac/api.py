from __future__ import annotations

import math
from typing import Optional

import numpy as np

from . import _core
from .config import TrackingConfig
from .grid import Grid
from .result import TrackResult

_METHOD_CODE = {"xcor": 0, "ncov": 1}
_FMISS = -999.0
_IMISS = -999
_AUTO_ZMISS = np.float32(-3.4028235e38)  # near -FLT_MAX; used only when z has NaN but no missing_value was given


class Tracker:
    """A reusable set of tracking parameters; call `.track()` with different data."""

    def __init__(
        self,
        template,
        *,
        search_velocity=None,
        search_radius=None,
        nsteps=2,
        method="xcor",
        subgrid="paraboloid",
        step=1,
        min_score=(0.8, 0.7),
        max_velocity_change=None,
        min_contrast=None,
        min_peak_prominence=None,
        fixed_template=False,
        min_samples=1,
        workers=None,
    ):
        self.config = TrackingConfig(
            template=template,
            search_velocity=search_velocity,
            search_radius=search_radius,
            nsteps=nsteps,
            method=method,
            subgrid=subgrid,
            step=step,
            min_score=min_score,
            max_velocity_change=max_velocity_change,
            min_contrast=min_contrast,
            min_peak_prominence=min_peak_prominence,
            fixed_template=fixed_template,
            min_samples=min_samples,
            workers=workers,
        )

    def track(
        self,
        z,
        x0,
        y0,
        t0=0,
        *,
        step=None,
        time=None,
        mask=None,
        first_guess=None,
        missing_value=None,
        grid=None,
        diagnostics=False,
    ) -> TrackResult:
        """`step=None` (default) uses the `step` this `Tracker` was configured
        with; passing `step` here overrides it for this call only, without
        mutating `self.config`. Useful for running forward (`step=+1`) and
        backward (`step=-1`) tracking from a single `Tracker`."""
        return _track_impl(
            self.config,
            z,
            x0,
            y0,
            t0,
            step=step,
            time=time,
            mask=mask,
            first_guess=first_guess,
            missing_value=missing_value,
            grid=grid,
            diagnostics=diagnostics,
        )


def track(
    z,
    x0,
    y0,
    t0=0,
    *,
    template,
    search_velocity=None,
    search_radius=None,
    nsteps=2,
    time=None,
    mask=None,
    method="xcor",
    subgrid="paraboloid",
    step=1,
    min_score=(0.8, 0.7),
    max_velocity_change=None,
    min_contrast=None,
    min_peak_prominence=None,
    fixed_template=False,
    min_samples=1,
    first_guess=None,
    missing_value=None,
    grid=None,
    diagnostics=False,
    workers=None,
) -> TrackResult:
    """Track template patches through a `(nt, ny, nx)` image sequence.

    See docs/api.md for the full parameter reference; docs/migration-v1-to-v2.md
    for the mapping from the removed v1 `VTT`/`setup`/`trac` API.
    """
    config = TrackingConfig(
        template=template,
        search_velocity=search_velocity,
        search_radius=search_radius,
        nsteps=nsteps,
        method=method,
        subgrid=subgrid,
        step=step,
        min_score=min_score,
        max_velocity_change=max_velocity_change,
        min_contrast=min_contrast,
        min_peak_prominence=min_peak_prominence,
        fixed_template=fixed_template,
        min_samples=min_samples,
        workers=workers,
    )
    return _track_impl(
        config,
        z,
        x0,
        y0,
        t0,
        step=None,
        time=time,
        mask=mask,
        first_guess=first_guess,
        missing_value=missing_value,
        grid=grid,
        diagnostics=diagnostics,
    )


def search_radius_from_velocity(search_velocity, dt, grid=None) -> "tuple[int, int]":
    """Pixel search radius `(iy, ix)` from a `(vy, vx)` search velocity and a
    reference time interval `dt`.

    If `grid` is given, `search_velocity` is in physical units; otherwise
    it's in index units. Same formula `track()` uses internally when
    `search_radius` isn't given explicitly: `ceil(abs(v_index * dt)) + 1`.
    """
    vy, vx = search_velocity
    if grid is not None:
        vy = grid.velocity_to_index_y(vy)
        vx = grid.velocity_to_index_x(vx)
    iy = math.ceil(abs(float(vy) * dt)) + 1
    ix = math.ceil(abs(float(vx) * dt)) + 1
    return iy, ix


def _parse_diagnostics(diagnostics) -> "tuple[bool, bool]":
    """Returns `(want_templates, want_score_grids)`.

    `diagnostics` accepts: `False` (neither, default), `True` (both),
    `"templates"`, `"score_grids"`, or a tuple/list containing either or
    both of those two strings.
    """
    if diagnostics is False:
        return False, False
    if diagnostics is True:
        return True, True
    if isinstance(diagnostics, str):
        names = (diagnostics,)
    else:
        names = tuple(diagnostics)
    valid = {"templates", "score_grids"}
    unknown = set(names) - valid
    if unknown:
        raise ValueError(
            f"`diagnostics` contains unknown value(s) {sorted(unknown)!r}; "
            f"expected a subset of {sorted(valid)!r}, True, or False"
        )
    return "templates" in names, "score_grids" in names


def _resolve_auto_grid(z, grid):
    if not (isinstance(grid, str) and grid == "auto"):
        return grid
    try:
        import xarray as xr
    except ImportError as e:
        raise ImportError('grid="auto" requires the optional "xarray" dependency') from e
    if not isinstance(z, xr.DataArray):
        raise ValueError('grid="auto" requires `z` to be an xarray.DataArray')
    if z.ndim != 3:
        raise ValueError("`z` must be 3-dimensional (time, y, x)")
    y_dim, x_dim = z.dims[-2], z.dims[-1]
    return Grid.from_coords(z.coords[x_dim].values, z.coords[y_dim].values)


def _track_impl(
    config: TrackingConfig,
    z,
    x0,
    y0,
    t0,
    *,
    step,
    time,
    mask,
    first_guess,
    missing_value,
    grid,
    diagnostics: bool,
) -> TrackResult:
    effective_step = config.step if step is None else step
    if effective_step == 0:
        raise ValueError("`step` must not be 0")

    grid = _resolve_auto_grid(z, grid)

    z = np.asarray(z)
    if z.ndim != 3:
        raise ValueError("`z` must be 3-dimensional (nt, ny, nx)")
    nt, ny, nx = z.shape
    if nt < 2:
        raise ValueError("`z` must have at least 2 time steps (z.shape[0] >= 2)")

    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.shape != z.shape:
            raise ValueError("`mask` must have the same shape as `z`")

    x0 = np.asarray(x0, dtype=np.float64)
    y0 = np.asarray(y0, dtype=np.float64)
    if x0.shape != y0.shape:
        raise ValueError("`x0` and `y0` must have the same shape")
    seed_shape = x0.shape
    n = x0.size

    if np.ndim(t0) == 0:
        tid0 = np.full(seed_shape, int(t0), dtype=np.int64)
    else:
        tid0 = np.asarray(t0, dtype=np.int64)
        if tid0.shape != seed_shape:
            raise ValueError("`t0` must be a scalar or match the shape of `x0`/`y0`")

    if first_guess is None:
        vyg = np.zeros(seed_shape, dtype=np.float64)
        vxg = np.zeros(seed_shape, dtype=np.float64)
    else:
        vyg0, vxg0 = first_guess
        vyg = np.broadcast_to(np.asarray(vyg0, dtype=np.float64), seed_shape).copy()
        vxg = np.broadcast_to(np.asarray(vxg0, dtype=np.float64), seed_shape).copy()

    if time is None:
        t = np.arange(nt, dtype=np.float64)
    else:
        t = np.asarray(time, dtype=np.float64)
        if t.shape != (nt,):
            raise ValueError("`time` must have shape (nt,)")

    use_grid = isinstance(grid, Grid)
    if use_grid:
        x0_idx = grid.to_index_x(x0)
        y0_idx = grid.to_index_y(y0)
        vyg = grid.velocity_to_index_y(vyg)
        vxg = grid.velocity_to_index_x(vxg)
        svy, svx = (
            (grid.velocity_to_index_y(config.search_velocity[0]),
             grid.velocity_to_index_x(config.search_velocity[1]))
            if config.search_velocity is not None
            else (None, None)
        )
        dvy, dvx = (
            (float(grid.velocity_to_index_y(config.max_velocity_change[0])),
             float(grid.velocity_to_index_x(config.max_velocity_change[1])))
            if config.max_velocity_change is not None
            else (None, None)
        )
    else:
        x0_idx, y0_idx = x0, y0
        svy, svx = config.search_velocity if config.search_velocity is not None else (None, None)
        dvy, dvx = (
            config.max_velocity_change if config.max_velocity_change is not None else (None, None)
        )

    nsy, nsx = config.template

    if config.search_radius is not None:
        iyhw, ixhw = config.search_radius
    else:
        iyhw, ixhw = search_radius_from_velocity(
            (svy, svx), dt=effective_step * (t[-1] - t[0]) / (nt - 1)
        )

    use_zmiss = missing_value is not None
    zmiss = np.float32(missing_value) if use_zmiss else np.float32(0.0)

    use_mask = mask is not None and bool(mask.any())
    if use_mask:
        visible_f = np.ascontiguousarray(~mask).transpose(2, 1, 0)
    else:
        # Fortran always expects a (nx,ny,nt)-shaped array (never indexed when
        # use_mask is False); a small allocation is the price of a single,
        # simple Fortran signature.
        visible_f = np.ones((nx, ny, nt), dtype=np.bool_, order="F")

    zf = np.ascontiguousarray(z, dtype=np.float32).transpose(2, 1, 0)  # zero-copy: (nx,ny,nt), F-contig
    nan_mask = np.isnan(zf)
    if nan_mask.any():
        # NaN is always treated as missing, in addition to (and made
        # equivalent to) any explicit `missing_value`: both collapse onto
        # the same Fortran-side zmiss sentinel, so the two are indistinguishable
        # downstream. This necessarily breaks the zero-copy fast path (we
        # must not mutate the caller's `z` in place).
        if not use_zmiss:
            zmiss = _AUTO_ZMISS
            use_zmiss = True
        zf = zf.copy(order="F")
        zf[nan_mask] = zmiss

    tid0_f = tid0.reshape(-1) + 1  # 0-based -> Fortran 1-based
    x0_f = x0_idx.reshape(-1).astype(np.float64) + 1.0
    y0_f = y0_idx.reshape(-1).astype(np.float64) + 1.0
    vxg_f = np.ascontiguousarray(vxg.reshape(-1), dtype=np.float64)
    vyg_f = np.ascontiguousarray(vyg.reshape(-1), dtype=np.float64)

    method_code = _METHOD_CODE[config.method]
    # `max_velocity_change` is a magnitude; take abs() here rather than
    # requiring dvy/dvx > 0, since a Grid with negative dy/dx (a descending
    # coordinate axis) flips the sign of the physical->index conversion
    # above even though the physical value was validated positive.
    use_vch = dvy is not None and dvx is not None
    vxch = abs(dvx) if dvx is not None else -999.0
    vych = abs(dvy) if dvy is not None else -999.0
    use_peak_th = config.min_peak_prominence is not None and config.min_peak_prominence > 0
    peak_th = np.float32(config.min_peak_prominence) if use_peak_th else np.float32(-1.0)
    use_cth = config.min_contrast is not None and config.min_contrast > 0
    cth = np.float32(config.min_contrast) if use_cth else np.float32(-1.0)

    sth0, sth1 = config.min_score
    nsteps = config.nsteps

    want_templates, want_score_grids = _parse_diagnostics(diagnostics)
    if want_templates:
        # Steps never reached (tracking stopped earlier) are left at this
        # fill value, matching VTTrac.jl's `fill(o.zmiss, ...)` initial
        # allocation exactly, so golden-data comparisons don't need to
        # special-case unwritten regions.
        zss_fill = zmiss if use_zmiss else np.float32(0.0)
        zss_buf = np.full((nsx, nsy, nsteps + 1, n), zss_fill, dtype=np.float32, order="F")
    else:
        # Fortran always expects a 4-D array (never indexed when want_zss is
        # False); a minimal dummy keeps the call signature uniform.
        zss_buf = np.zeros((1, 1, 1, 1), dtype=np.float32, order="F")
    if want_score_grids:
        scr_buf = np.full(
            (2 * ixhw + 1, 2 * iyhw + 1, nsteps, n), np.float32(_FMISS), dtype=np.float32, order="F"
        )
    else:
        scr_buf = np.zeros((1, 1, 1, 1), dtype=np.float32, order="F")

    nthreads = 0 if config.workers is None else int(config.workers)

    cnt, status, tid, xo, yo, vxo, vyo, score = _core.pyvttrac_core.track(
        zf, t, use_zmiss, zmiss, use_mask, visible_f, config.min_samples,
        nsx, nsy, ixhw, iyhw, effective_step, nsteps,
        tid0_f, x0_f, y0_f, vxg_f, vyg_f,
        method_code, config.use_subgrid, config.use_gaussian_subgrid,
        sth0, sth1, use_vch, vxch, vych, use_peak_th, peak_th, use_cth, cth,
        config.fixed_template, _FMISS, _IMISS, nthreads,
        want_templates, zss_buf, want_score_grids, scr_buf,
    )

    t_index = np.where(tid == _IMISS, -1, tid - 1).astype(np.int64)
    x_out = np.where(xo == _FMISS, np.nan, xo - 1.0)
    y_out = np.where(yo == _FMISS, np.nan, yo - 1.0)
    vx_out = np.where(vxo == _FMISS, np.nan, vxo)
    vy_out = np.where(vyo == _FMISS, np.nan, vyo)
    score_out = np.where(score == _FMISS, np.nan, score)

    if use_grid:
        x_out = grid.to_phys_x(x_out)
        y_out = grid.to_phys_y(y_out)
        vx_out = grid.velocity_to_phys_x(vx_out)
        vy_out = grid.velocity_to_phys_y(vy_out)

    count_out = cnt.reshape(seed_shape)
    status_out = status.reshape(seed_shape)
    t_index = t_index.reshape(nsteps + 1, *seed_shape)
    x_out = x_out.reshape(nsteps + 1, *seed_shape)
    y_out = y_out.reshape(nsteps + 1, *seed_shape)
    vx_out = vx_out.reshape(nsteps, *seed_shape)
    vy_out = vy_out.reshape(nsteps, *seed_shape)
    score_out = score_out.reshape(nsteps, *seed_shape)

    templates = None
    score_grids = None
    if want_templates:
        templates = np.ascontiguousarray(zss_buf.transpose(2, 1, 0, 3)).reshape(
            nsteps + 1, nsy, nsx, *seed_shape
        )
    if want_score_grids:
        score_grids = np.ascontiguousarray(scr_buf.transpose(2, 1, 0, 3)).reshape(
            nsteps, 2 * iyhw + 1, 2 * ixhw + 1, *seed_shape
        )

    return TrackResult(
        count=count_out,
        status=status_out,
        t_index=t_index,
        x=x_out,
        y=y_out,
        vx=vx_out,
        vy=vy_out,
        score=score_out,
        templates=templates,
        score_grids=score_grids,
        step=effective_step,
    )
