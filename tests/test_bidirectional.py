"""concat_bidirectional(): combining a forward- and backward-tracked
TrackResult from a shared origin into one trajectory.
"""

import numpy as np
import pytest

import pyvttrac as vt
from pyvttrac.status import Status


def _wave_field(nt=12, ny=60, nx=60, k=2 * np.pi / 12, cx=1.2, cy=1.2):
    t = np.arange(nt, dtype=np.float64)
    tg, yg, xg = np.meshgrid(t, np.arange(ny), np.arange(nx), indexing="ij")
    z = np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))
    return z.astype(np.float32)


def _forward_backward(t0=5, nsteps_fwd=3, nsteps_bwd=3, diagnostics=False, x0=None, y0=None):
    z = _wave_field()
    if x0 is None:
        x0 = np.array([30.0])
    if y0 is None:
        y0 = np.array([30.0])
    kw = dict(template=(7, 7), search_velocity=(2.0, 2.0), min_score=(0.3, 0.3))
    tracker_f = vt.Tracker(nsteps=nsteps_fwd, **kw)
    tracker_b = vt.Tracker(nsteps=nsteps_bwd, **kw)
    f = tracker_f.track(z, x0, y0, t0=t0, step=1, diagnostics=diagnostics)
    b = tracker_b.track(z, x0, y0, t0=t0, step=-1, diagnostics=diagnostics)
    return f, b


def test_vx_sign_does_not_flip_across_the_combined_trajectory():
    # The critical regression this function exists to get right: reversing
    # the backward leg must NOT sign-flip its velocities.
    f, b = _forward_backward()
    c = vt.concat_bidirectional(f, b)
    assert np.all(np.sign(c.vx) == np.sign(c.vx[0]))
    assert np.all(np.sign(c.vy) == np.sign(c.vy[0]))
    # and it should match the (positive) forward direction, not be flipped
    assert np.sign(c.vx[0]) == np.sign(f.vx[0])


def test_position_is_monotonic_for_uniformly_translating_field():
    f, b = _forward_backward()
    c = vt.concat_bidirectional(f, b)
    dx = np.diff(c.x[:, 0])
    assert np.all(dx > 0)  # uniformly rightward-moving wave


def test_t_index_is_monotonically_increasing():
    f, b = _forward_backward()
    c = vt.concat_bidirectional(f, b)
    assert np.all(np.diff(c.t_index[:, 0]) == 1)


def test_shapes_drop_shared_origin_true():
    nf, nb = 3, 4
    f, b = _forward_backward(nsteps_fwd=nf, nsteps_bwd=nb)
    c = vt.concat_bidirectional(f, b, drop_shared_origin=True)
    assert c.x.shape == (nb + nf + 1, *f.count.shape)
    assert c.t_index.shape == (nb + nf + 1, *f.count.shape)
    assert c.vx.shape == (nb + nf, *f.count.shape)
    assert c.score.shape == (nb + nf, *f.count.shape)


def test_shapes_drop_shared_origin_false():
    nf, nb = 3, 4
    f, b = _forward_backward(nsteps_fwd=nf, nsteps_bwd=nb)
    c_drop = vt.concat_bidirectional(f, b, drop_shared_origin=True)
    c_keep = vt.concat_bidirectional(f, b, drop_shared_origin=False)
    assert c_keep.x.shape[0] == c_drop.x.shape[0] + 1
    # the origin appears twice, back-to-back, at indices nb (backward's own
    # copy) and nb+1 (forward's own copy) in the kept version
    np.testing.assert_array_equal(c_keep.x[nb], c_keep.x[nb + 1])
    np.testing.assert_array_equal(c_keep.t_index[nb], c_keep.t_index[nb + 1])
    # velocity axis length is unaffected by drop_shared_origin
    assert c_keep.vx.shape == c_drop.vx.shape


def test_count_status_and_step():
    nf, nb = 3, 4
    f, b = _forward_backward(nsteps_fwd=nf, nsteps_bwd=nb)
    c = vt.concat_bidirectional(f, b)
    assert c.count[0] == nb + nf + 1
    assert c.status[0] == Status.OK
    assert c.status_forward[0] == f.status[0]
    assert c.status_backward[0] == b.status[0]
    assert c.step == f.step
    assert c.step > 0


def test_status_combination_rules():
    # forward fails (too-strict min_score), backward succeeds
    z = _wave_field()
    x0, y0 = np.array([30.0]), np.array([30.0])
    kw = dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=3)
    f_fail = vt.track(z, x0, y0, t0=5, step=1, min_score=(0.999, 0.999), **kw)
    b_ok = vt.track(z, x0, y0, t0=5, step=-1, min_score=(0.3, 0.3), **kw)
    assert f_fail.status[0] != Status.OK
    assert b_ok.status[0] == Status.OK

    c = vt.concat_bidirectional(f_fail, b_ok)
    # forward's non-OK status wins when only one side failed
    assert c.status[0] == f_fail.status[0]

    # both fail: forward's status wins (documented tie-break)
    f_fail2 = vt.track(z, x0, y0, t0=5, step=1, min_score=(0.999, 0.999), **kw)
    b_fail2 = vt.track(z, x0, y0, t0=5, step=-1, min_score=(0.998, 0.998), **kw)
    assert f_fail2.status[0] != Status.OK
    assert b_fail2.status[0] != Status.OK
    c2 = vt.concat_bidirectional(f_fail2, b_fail2)
    assert c2.status[0] == f_fail2.status[0]


def test_requires_opposite_step_signs():
    f, b = _forward_backward()
    with pytest.raises(ValueError):
        vt.concat_bidirectional(b, f)  # swapped: forward.step < 0
    with pytest.raises(ValueError):
        vt.concat_bidirectional(f, f)  # both positive


def test_requires_equal_step_magnitude():
    z = _wave_field()
    x0, y0 = np.array([30.0]), np.array([30.0])
    kw = dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=2)
    f = vt.track(z, x0, y0, t0=5, step=1, **kw)
    b = vt.track(z, x0, y0, t0=5, step=-2, **kw)
    with pytest.raises(ValueError):
        vt.concat_bidirectional(f, b)


def test_requires_matching_seed_shape():
    z = _wave_field()
    kw = dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=2)
    f = vt.track(z, np.array([30.0, 31.0]), np.array([30.0, 31.0]), t0=5, step=1, **kw)
    b = vt.track(z, np.array([30.0]), np.array([30.0]), t0=5, step=-1, **kw)
    with pytest.raises(ValueError):
        vt.concat_bidirectional(f, b)


def test_requires_shared_origin():
    z = _wave_field()
    kw = dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=2)
    f = vt.track(z, np.array([30.0]), np.array([30.0]), t0=5, step=1, **kw)
    b = vt.track(z, np.array([30.0]), np.array([30.0]), t0=6, step=-1, **kw)  # different t0
    with pytest.raises(ValueError):
        vt.concat_bidirectional(f, b)


def test_diagnostics_both_or_neither_required():
    f, b = _forward_backward(diagnostics=True)
    f_no_diag, b_no_diag = _forward_backward(diagnostics=False)
    with pytest.raises(ValueError):
        vt.concat_bidirectional(f, b_no_diag)
    with pytest.raises(ValueError):
        vt.concat_bidirectional(f_no_diag, b)


def test_diagnostics_shapes_when_both_present():
    nf, nb = 3, 4
    f, b = _forward_backward(nsteps_fwd=nf, nsteps_bwd=nb, diagnostics=True)
    c = vt.concat_bidirectional(f, b)
    assert c.templates.shape == (nb + nf + 1, *f.templates.shape[1:-1], *f.count.shape)
    assert c.score_grids.shape == (nb + nf, *f.score_grids.shape[1:-1], *f.count.shape)


def test_grid_and_provenance_carried_over_when_they_agree():
    grid = vt.Grid(x0=0.0, y0=0.0, dx=2.0, dy=2.0)
    z = _wave_field()
    x0_idx, y0_idx = np.array([30.0]), np.array([30.0])
    x0, y0 = grid.to_phys_x(x0_idx), grid.to_phys_y(y0_idx)
    kw = dict(template=(7, 7), search_velocity=(4.0, 4.0), nsteps=2, grid=grid)
    f = vt.track(z, x0, y0, t0=5, step=1, **kw)
    b = vt.track(z, x0, y0, t0=5, step=-1, **kw)

    c = vt.concat_bidirectional(f, b)
    assert c.grid == grid
    np.testing.assert_array_equal(c.seed_x, f.seed_x)
    np.testing.assert_array_equal(c.seed_y, f.seed_y)
    assert c.search_radius == f.search_radius


def test_provenance_is_none_when_forward_and_backward_disagree():
    z = _wave_field()
    x0, y0 = np.array([30.0]), np.array([30.0])
    f = vt.track(z, x0, y0, t0=5, template=(7, 7), search_radius=(3, 3), nsteps=2, step=1)
    b = vt.track(z, x0, y0, t0=5, template=(7, 7), search_radius=(4, 4), nsteps=2, step=-1)
    c = vt.concat_bidirectional(f, b)
    assert c.search_radius is None


def test_to_xarray_on_combined_result():
    pytest.importorskip("xarray")
    nf, nb = 3, 4
    f, b = _forward_backward(nsteps_fwd=nf, nsteps_bwd=nb)
    c = vt.concat_bidirectional(f, b)
    ds = c.to_xarray()
    np.testing.assert_array_equal(ds["step"].values, np.arange(-nb, nf + 1))
    np.testing.assert_allclose(ds["step_v"].values, np.arange(-nb, nf) + 0.5)
