from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from .status import Status


@dataclass
class TrackResult:
    """Result of `track()` / `Tracker.track()`.

    Shapes (seed_shape = the broadcast shape of the `x0`/`y0` seed arrays):
        count, status:            seed_shape
        t_index, x, y:            (nsteps+1, *seed_shape)
        vx, vy, score:            (nsteps,   *seed_shape)
        templates (diagnostics):  (nsteps+1, ny_sub, nx_sub, *seed_shape)
        score_grids (diagnostics):(nsteps,   ny_search, nx_search, *seed_shape)

    Invalid entries are NaN for float arrays and -1 for `t_index` (invalid
    `count`/`status` are never negative, so no sentinel is needed there).
    """

    count: np.ndarray
    status: np.ndarray
    t_index: np.ndarray
    x: np.ndarray
    y: np.ndarray
    vx: np.ndarray
    vy: np.ndarray
    score: np.ndarray
    templates: Optional[np.ndarray] = None
    score_grids: Optional[np.ndarray] = None
    step: int = 1

    @property
    def ok(self) -> np.ndarray:
        return self.status == Status.OK

    def to_xarray(self, *, units: Optional[dict] = None, global_attrs: Optional[dict] = None):
        """Convert to an `xarray.Dataset` with CF-1.13 / ACDD-1.3 metadata.

        The dimensional layout is unchanged from previous versions (`step`/
        `step_v` + one dim per seed axis) -- this only adds standards-compliant
        attributes on top, so existing code indexing the returned `Dataset`
        keeps working.

        What's filled in automatically:
        - Global: `Conventions`, a generic-but-accurate `title`/`summary`/
          `keywords`, `source`, `history`, `date_created`.
        - Per variable: `long_name`, `units` (default `"1"`, i.e.
          dimensionless/index units -- see `units` below), and `status`'s
          CF `flag_values`/`flag_meanings` (derived from the `Status` enum).
        - `_FillValue` is set via each variable's `.encoding` (not `.attrs`,
          which `Dataset.to_netcdf()` requires) so NaN/-1 sentinels round-trip
          through netCDF correctly.

        What it cannot know, and you should supply for real ACDD compliance:
        - Discovery attributes tied to *your* deployment (`institution`,
          `creator_name`, `creator_email`, `license`, `id`,
          `naming_authority`, `project`, ...) -- pass via `global_attrs`.
        - Physical units for `x`/`y`/`vx`/`vy`/`templates` when you tracked in
          physical coordinates (via `grid=`) or `z` has known units --
          `TrackResult` doesn't retain the `Grid` used, so pass e.g.
          `units={"x": "m", "y": "m", "vx": "m s-1", "vy": "m s-1"}`.
          `count`/`status`/`t_index`/`score`/`step`/`step_v` are always
          dimensionless and shouldn't need overriding.

        Parameters
        ----------
        units : dict, optional
            Maps variable/coordinate name -> CF/UDUNITS unit string, merged
            over the `"1"` (dimensionless) defaults.
        global_attrs : dict, optional
            Extra/overriding dataset-level attributes (e.g. ACDD discovery
            attributes), merged over the auto-generated ones.
        """
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "TrackResult.to_xarray() requires the optional 'xarray' dependency. "
                'Install it with `pip install "pyVTTrac[xarray]"`.'
            ) from e

        from datetime import datetime, timezone

        from . import __version__

        seed_shape = self.count.shape
        seed_dims = [f"seed_{i}" for i in range(len(seed_shape))]
        nsteps = self.vx.shape[0]

        data_vars = dict(
            count=(seed_dims, self.count),
            status=(seed_dims, self.status),
            t_index=(["step", *seed_dims], self.t_index),
            x=(["step", *seed_dims], self.x),
            y=(["step", *seed_dims], self.y),
            vx=(["step_v", *seed_dims], self.vx),
            vy=(["step_v", *seed_dims], self.vy),
            score=(["step_v", *seed_dims], self.score),
        )
        coords = dict(
            step=np.arange(nsteps + 1) * self.step,
            step_v=np.arange(nsteps) * self.step + 0.5 * np.sign(self.step),
        )
        if self.templates is not None:
            data_vars["templates"] = (["step", "ny_sub", "nx_sub", *seed_dims], self.templates)
        if self.score_grids is not None:
            data_vars["score_grids"] = (
                ["step_v", "ny_search", "nx_search", *seed_dims],
                self.score_grids,
            )
        ds = xr.Dataset(data_vars=data_vars, coords=coords)

        long_names = {
            "count": "number of valid trajectory points, including the initial position",
            "status": "tracking status code",
            "t_index": "time index into the tracked data's time axis",
            "x": "x position",
            "y": "y position",
            "vx": "x-component velocity between consecutive tracked points",
            "vy": "y-component velocity between consecutive tracked points",
            "score": "template-matching score at the tracked position",
            "templates": "template sub-image at each tracking step",
            "score_grids": "template-matching score over the full search window",
            "step": "tracking step index relative to the initial position, scaled by `step`",
            "step_v": "tracking step index at which each velocity/score applies",
        }
        default_units = {name: "1" for name in long_names}  # dimensionless/index units by default
        merged_units = {**default_units, **(units or {})}

        for name in ds.variables:
            if name in long_names:
                ds[name].attrs["long_name"] = long_names[name]
            if name in merged_units and name != "status":  # a flag variable has no `units`
                ds[name].attrs["units"] = merged_units[name]

        ds["status"].attrs["flag_values"] = np.array([s.value for s in Status], dtype=np.int32)
        ds["status"].attrs["flag_meanings"] = " ".join(s.name for s in Status)

        ds["t_index"].encoding["_FillValue"] = -1
        for name in ("x", "y", "vx", "vy", "score", "templates", "score_grids"):
            if name in ds:
                ds[name].encoding["_FillValue"] = np.nan

        now = datetime.now(timezone.utc).replace(microsecond=0).isoformat()
        ds.attrs.update(
            {
                "Conventions": "CF-1.13, ACDD-1.3",
                "title": "Template tracking result (pyvttrac)",
                "summary": (
                    "Trajectories of tracked template positions and velocities produced by "
                    "pyvttrac.track(), a PIV/PTV-style velocimetry-by-template-tracking algorithm."
                ),
                "keywords": "velocimetry, template tracking, PIV, PTV, trajectory",
                "source": "pyvttrac.track()",
                "history": f"{now}: created by pyvttrac {__version__}",
                "date_created": now,
            }
        )
        if global_attrs:
            ds.attrs.update(global_attrs)
        return ds

    def to_dataframe(self):
        """Long-format DataFrame: one row per (seed, step)."""
        try:
            import pandas as pd
        except ImportError as e:
            raise ImportError(
                "TrackResult.to_dataframe() requires the optional 'pandas' dependency. "
                'Install it with `pip install "pyVTTrac[pandas]"`.'
            ) from e

        nsteps = self.vx.shape[0]
        count = self.count.reshape(-1)
        status = self.status.reshape(-1)
        t_index = self.t_index.reshape(nsteps + 1, -1)
        x = self.x.reshape(nsteps + 1, -1)
        y = self.y.reshape(nsteps + 1, -1)
        vx = self.vx.reshape(nsteps, -1)
        vy = self.vy.reshape(nsteps, -1)
        score = self.score.reshape(nsteps, -1)

        n_seed = count.shape[0]
        seed_idx = np.repeat(np.arange(n_seed), nsteps + 1)
        step_idx = np.tile(np.arange(nsteps + 1), n_seed)

        vx_col = np.full(n_seed * (nsteps + 1), np.nan)
        vy_col = np.full(n_seed * (nsteps + 1), np.nan)
        score_col = np.full(n_seed * (nsteps + 1), np.nan)
        row = 0
        for m in range(n_seed):
            for step in range(nsteps + 1):
                if step >= 1:
                    vx_col[row] = vx[step - 1, m]
                    vy_col[row] = vy[step - 1, m]
                    score_col[row] = score[step - 1, m]
                row += 1

        return pd.DataFrame(
            dict(
                seed=seed_idx,
                step=step_idx,
                count=count[seed_idx],
                status=status[seed_idx],
                t_index=t_index[step_idx, seed_idx],
                x=x[step_idx, seed_idx],
                y=y[step_idx, seed_idx],
                vx=vx_col,
                vy=vy_col,
                score=score_col,
            )
        )
