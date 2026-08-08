"""Grid: physical <-> index coordinate round-trips, unit_factor, and
non-equally-spaced-coordinate rejection.
"""

import numpy as np
import pytest

from pyvttrac.grid import Grid


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


def test_dx_dy_are_stored_as_absolute_value():
    grid = Grid(x0=0.0, y0=0.0, dx=-2.0, dy=-3.0)
    assert grid.dx == 2.0
    assert grid.dy == 3.0


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
