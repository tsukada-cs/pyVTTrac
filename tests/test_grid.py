"""Grid: physical <-> index coordinate round-trips, unit_factor, and
non-equally-spaced-coordinate rejection.
"""

import numpy as np
import pytest

import pyvttrac as vt
from pyvttrac.grid import Grid
from pyvttrac.status import Status


def test_index_to_phys_roundtrip():
    grid = Grid(x0=100.0, y0=-50.0, dx=2.0, dy=3.0)
    x_index = np.array([0.0, 1.5, 10.0, -2.0])
    y_index = np.array([0.0, -1.5, 5.0, 20.0])

    x_phys = grid.to_phys_x(x_index)
    y_phys = grid.to_phys_y(y_index)
    np.testing.assert_allclose(grid.to_index_x(x_phys), x_index)
    np.testing.assert_allclose(grid.to_index_y(y_phys), y_index)


def test_phys_to_index_roundtrip():
    grid = Grid(x0=0.0, y0=0.0, dx=0.5, dy=1.25)
    x_phys = np.array([0.0, 10.0, -3.25, 100.0])
    y_phys = np.array([0.0, 5.0, -1.25, 50.0])

    x_index = grid.to_index_x(x_phys)
    y_index = grid.to_index_y(y_phys)
    np.testing.assert_allclose(grid.to_phys_x(x_index), x_phys)
    np.testing.assert_allclose(grid.to_phys_y(y_index), y_phys)


def test_dx_dy_keep_their_sign():
    grid = Grid(x0=0.0, y0=0.0, dx=-2.0, dy=-3.0)
    assert grid.dx == -2.0
    assert grid.dy == -3.0


def test_dx_dy_reject_zero():
    with pytest.raises(ValueError):
        Grid(x0=0.0, y0=0.0, dx=0.0, dy=1.0)
    with pytest.raises(ValueError):
        Grid(x0=0.0, y0=0.0, dx=1.0, dy=0.0)


def test_roundtrip_with_descending_coordinate():
    # A descending y axis (e.g. latitude) has negative dy; index<->phys must
    # still round-trip, and to_index_y must return the correct sign (this is
    # the exact regression covered by Grid.from_coords([10,9,8,7], ...)).
    grid = Grid.from_coords(x_coord=np.array([0.0, 1.0, 2.0, 3.0]), y_coord=np.array([10.0, 9.0, 8.0, 7.0]))
    assert grid.dy == pytest.approx(-1.0)
    assert grid.to_index_y(9.0) == pytest.approx(1.0)

    y_index = np.array([0.0, 1.0, 2.0, 3.0])
    y_phys = grid.to_phys_y(y_index)
    np.testing.assert_allclose(y_phys, [10.0, 9.0, 8.0, 7.0])
    np.testing.assert_allclose(grid.to_index_y(y_phys), y_index)


def test_max_velocity_change_is_effective_with_negative_dy():
    # A descending y axis must not silently disable max_velocity_change
    # (this is the `dvy > 0` bug: a negative dy makes the physical->index
    # conversion of a positive max_velocity_change come out negative).
    z = _wave_field()
    grid = Grid(x0=0.0, y0=59.0, dx=1.0, dy=-1.0)
    x0, y0 = grid.to_phys_x(np.array([15.0])), grid.to_phys_y(np.array([15.0]))

    res_free = vt.track(
        z, x0, y0, t0=0, template=(7, 7), search_velocity=(6.0, 6.0), nsteps=2, grid=grid,
    )
    res_capped = vt.track(
        z, x0, y0, t0=0, template=(7, 7), search_velocity=(6.0, 6.0), nsteps=2, grid=grid,
        max_velocity_change=(0.01, 0.01),
    )
    assert res_capped.status[0] == Status.VELOCITY_CHANGE_TOO_LARGE
    assert res_free.status[0] != Status.VELOCITY_CHANGE_TOO_LARGE


def _wave_field(nt=8, ny=60, nx=60, k=2 * np.pi / 12, cx=1.5, cy=1.5):
    t = np.arange(nt, dtype=np.float64)
    tg, yg, xg = np.meshgrid(t, np.arange(ny), np.arange(nx), indexing="ij")
    z = np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))
    return z.astype(np.float32)


def test_velocity_roundtrip_with_unit_factor():
    grid = Grid(x0=0.0, y0=0.0, dx=2.0, dy=2.0, unit_factor=1000.0)
    vx_phys = np.array([0.0, 50.0, -25.0])
    vy_phys = np.array([0.0, -10.0, 40.0])

    vx_index = grid.velocity_to_index_x(vx_phys)
    vy_index = grid.velocity_to_index_y(vy_phys)
    np.testing.assert_allclose(grid.velocity_to_phys_x(vx_index), vx_phys)
    np.testing.assert_allclose(grid.velocity_to_phys_y(vy_index), vy_phys)

    # unit_factor actually changes the conversion (index space is smaller by
    # `unit_factor` for a given physical velocity)
    grid_no_uf = Grid(x0=0.0, y0=0.0, dx=2.0, dy=2.0, unit_factor=1.0)
    assert grid.velocity_to_index_x(50.0) == pytest.approx(grid_no_uf.velocity_to_index_x(50.0) / 1000.0)


def test_unit_factor_default_is_one():
    grid = Grid(x0=0.0, y0=0.0, dx=2.0, dy=2.0)
    vx_phys = 10.0
    np.testing.assert_allclose(grid.velocity_to_index_x(vx_phys), vx_phys / 2.0)


def test_from_coords_infers_origin_and_spacing():
    x_coord = np.array([10.0, 12.0, 14.0, 16.0])
    y_coord = np.array([-5.0, -3.0, -1.0])
    grid = Grid.from_coords(x_coord, y_coord)
    assert grid.x0 == 10.0
    assert grid.y0 == -5.0
    assert grid.dx == pytest.approx(2.0)
    assert grid.dy == pytest.approx(2.0)


def test_from_coords_rejects_non_equally_spaced():
    x_coord = np.array([0.0, 1.0, 2.0, 5.0])  # not equally spaced
    y_coord = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError):
        Grid.from_coords(x_coord, y_coord)


def test_from_coords_requires_at_least_two_points():
    with pytest.raises(ValueError):
        Grid.from_coords(np.array([1.0]), np.array([0.0, 1.0]))
