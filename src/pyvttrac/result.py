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

    def to_xarray(self):
        try:
            import xarray as xr
        except ImportError as e:
            raise ImportError(
                "TrackResult.to_xarray() requires the optional 'xarray' dependency. "
                'Install it with `pip install "pyVTTrac[xarray]"`.'
            ) from e

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
        return xr.Dataset(data_vars=data_vars, coords=coords)

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
