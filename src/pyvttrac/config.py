from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple, Union

import numpy as np

Pair = Tuple[float, float]


def _float_pair(value, name: str) -> Pair:
    if np.isscalar(value):
        return (float(value), float(value))
    try:
        a, b = value
    except (TypeError, ValueError) as e:
        raise ValueError(f"`{name}` must be a scalar or a 2-tuple") from e
    return (float(a), float(b))


def _int_pair(value, name: str) -> Tuple[int, int]:
    if np.isscalar(value):
        return (int(value), int(value))
    try:
        a, b = value
    except (TypeError, ValueError) as e:
        raise ValueError(f"`{name}` must be a scalar or a 2-tuple") from e
    return (int(a), int(b))


@dataclass(frozen=True)
class TrackingConfig:
    """Reusable tracking parameters (everything `Tracker` bundles up-front).

    Data-specific arguments (z, x0/y0/t0, time, mask, first_guess,
    missing_value, grid, diagnostics) are passed per-call to `track()`
    instead -- see api.py.
    """

    template: Tuple[int, int]  # (ny, nx)
    search_velocity: Optional[Pair] = None  # (vy, vx), index units unless a Grid is used
    search_radius: Optional[Tuple[int, int]] = None  # (iy, ix), pixels
    nsteps: int = 2
    method: str = "xcor"
    subgrid: Optional[str] = "paraboloid"  # "paraboloid" | "gaussian" | None
    step: int = 1
    min_score: Union[float, Pair] = (0.8, 0.7)
    max_velocity_change: Optional[Pair] = None  # (dvy, dvx)
    min_contrast: Optional[float] = None
    min_peak_prominence: Optional[float] = None
    fixed_template: bool = False
    min_samples: int = 1
    workers: Optional[int] = None

    def __post_init__(self):
        object.__setattr__(self, "template", _int_pair(self.template, "template"))
        object.__setattr__(self, "min_score", _float_pair(self.min_score, "min_score"))
        if self.search_velocity is not None:
            object.__setattr__(
                self, "search_velocity", _float_pair(self.search_velocity, "search_velocity")
            )
        if self.search_radius is not None:
            object.__setattr__(
                self, "search_radius", _int_pair(self.search_radius, "search_radius")
            )
        if self.max_velocity_change is not None:
            object.__setattr__(
                self, "max_velocity_change", _float_pair(self.max_velocity_change, "max_velocity_change")
            )

        nsy, nsx = self.template
        if nsy < 1 or nsx < 1:
            raise ValueError("`template` sizes must each be >= 1")
        if self.min_peak_prominence is not None and (nsy < 3 or nsx < 3):
            raise ValueError(
                "`template` sizes must each be >= 3 when `min_peak_prominence` is set"
            )
        if self.method not in ("xcor", "ncov"):
            raise ValueError('`method` must be "xcor" or "ncov"')
        if self.subgrid not in ("paraboloid", "gaussian", None):
            raise ValueError('`subgrid` must be "paraboloid", "gaussian", or None')
        if (self.search_velocity is None) == (self.search_radius is None):
            raise ValueError(
                "exactly one of `search_velocity` or `search_radius` must be specified"
            )
        if self.step == 0:
            raise ValueError("`step` must not be 0")
        if self.nsteps < 1:
            raise ValueError("`nsteps` must be >= 1")
        if self.min_samples < 1:
            raise ValueError("`min_samples` must be >= 1")
        if self.max_velocity_change is not None:
            dvy, dvx = self.max_velocity_change
            if dvy <= 0 or dvx <= 0:
                raise ValueError("`max_velocity_change` components must be positive")

    @property
    def use_subgrid(self) -> bool:
        return self.subgrid is not None

    @property
    def use_gaussian_subgrid(self) -> bool:
        return self.subgrid == "gaussian"
