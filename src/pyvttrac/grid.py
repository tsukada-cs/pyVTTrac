from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Grid:
    """Physical-coordinate <-> index-coordinate mapping for `track(..., grid=...)`.

    Replaces v1's `set_grid_par` / `setup_eq_grid` / `trac_eq_grid` /
    `calc_ixyhw_from_v_eq_grid`. Conversion formulas match those v1 methods:
    - index = (phys - origin) / spacing
    - index velocity = phys velocity / (spacing * unit_factor)
    """

    x0: float
    y0: float
    dx: float
    dy: float
    unit_factor: float = 1.0

    def __post_init__(self):
        if self.dx == 0 or self.dy == 0:
            raise ValueError("Grid.dx and Grid.dy must be non-zero")
        object.__setattr__(self, "dx", abs(self.dx))
        object.__setattr__(self, "dy", abs(self.dy))

    def to_index_x(self, x_phys):
        return (np.asarray(x_phys, dtype=np.float64) - self.x0) / self.dx

    def to_index_y(self, y_phys):
        return (np.asarray(y_phys, dtype=np.float64) - self.y0) / self.dy

    def to_phys_x(self, x_index):
        return np.asarray(x_index, dtype=np.float64) * self.dx + self.x0

    def to_phys_y(self, y_index):
        return np.asarray(y_index, dtype=np.float64) * self.dy + self.y0

    def velocity_to_index_x(self, vx_phys):
        return np.asarray(vx_phys, dtype=np.float64) / (self.dx * self.unit_factor)

    def velocity_to_index_y(self, vy_phys):
        return np.asarray(vy_phys, dtype=np.float64) / (self.dy * self.unit_factor)

    def velocity_to_phys_x(self, vx_index):
        return np.asarray(vx_index, dtype=np.float64) * self.dx * self.unit_factor

    def velocity_to_phys_y(self, vy_index):
        return np.asarray(vy_index, dtype=np.float64) * self.dy * self.unit_factor

    @staticmethod
    def from_coords(x_coord, y_coord) -> "Grid":
        """Infer a Grid from 1-D, equally-spaced coordinate arrays (used by `grid="auto"`)."""
        x_coord = np.asarray(x_coord, dtype=np.float64)
        y_coord = np.asarray(y_coord, dtype=np.float64)
        dx = _uniform_step(x_coord, "x")
        dy = _uniform_step(y_coord, "y")
        return Grid(x0=float(x_coord[0]), y0=float(y_coord[0]), dx=dx, dy=dy)


def _uniform_step(coord: np.ndarray, axis_name: str) -> float:
    if coord.size < 2:
        raise ValueError(f"{axis_name} coordinate must have at least 2 points to infer a Grid")
    diffs = np.diff(coord)
    step = diffs[0]
    if step == 0 or not np.allclose(diffs, step, rtol=1e-6):
        raise ValueError(
            f"{axis_name} coordinate is not equally spaced; pass an explicit `Grid` "
            'instead of grid="auto" (non-equally-spaced grids are outside the algorithm\'s scope)'
        )
    return float(step)
