"""Python API behavior: thread-safety/determinism, diagnostics not perturbing primary outputs, 
NaN<->missing_value equivalence, arbitrary seed shapes, Tracker vs track() parity, and to_xarray()/to_dataframe().
"""

import numpy as np
import pytest

import pyvttrac as vt
from pyvttrac.status import Status


def _wave_field(nt=10, ny=60, nx=60, k=2 * np.pi / 12, cx=1.2, cy=1.2):
    t = np.arange(nt, dtype=np.float64)
    tg, yg, xg = np.meshgrid(t, np.arange(ny), np.arange(nx), indexing="ij")
    z = np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))
    return z.astype(np.float32)


def _common_kwargs():
    return dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=4, min_score=(0.3, 0.3))


def test_workers_1_vs_8_bit_identical():
    z = _wave_field()
    x0, y0 = vt.seed_grid(z.shape, spacing=3, margin=8)
    res1 = vt.track(z, x0, y0, t0=0, workers=1, **_common_kwargs())
    res8 = vt.track(z, x0, y0, t0=0, workers=8, **_common_kwargs())

    np.testing.assert_array_equal(res1.count, res8.count)
    np.testing.assert_array_equal(res1.status, res8.status)
    np.testing.assert_array_equal(res1.t_index, res8.t_index)
    for name in ("x", "y", "vx", "vy", "score"):
        a, b = getattr(res1, name), getattr(res8, name)
        assert np.array_equal(a, b, equal_nan=True), f"{name} differs between workers=1 and workers=8"


def test_diagnostics_does_not_change_primary_outputs():
    z = _wave_field()
    x0, y0 = vt.seed_grid(z.shape, spacing=5, margin=8)
    res_plain = vt.track(z, x0, y0, t0=0, diagnostics=False, **_common_kwargs())
    res_diag = vt.track(z, x0, y0, t0=0, diagnostics=True, **_common_kwargs())

    np.testing.assert_array_equal(res_plain.count, res_diag.count)
    np.testing.assert_array_equal(res_plain.status, res_diag.status)
    np.testing.assert_array_equal(res_plain.t_index, res_diag.t_index)
    for name in ("x", "y", "vx", "vy", "score"):
        np.testing.assert_array_equal(getattr(res_plain, name), getattr(res_diag, name))

    assert res_plain.templates is None
    assert res_plain.score_grids is None
    assert res_diag.templates is not None
    assert res_diag.score_grids is not None


def test_nan_and_missing_value_are_equivalent():
    z = _wave_field(nt=6, ny=30, nx=30)
    # A 7x7 template centered on the seed will read this point on the first step.
    z_sentinel = z.copy()
    z_sentinel[0, 15, 15] = -999.0
    z_nan = z.copy()
    z_nan[0, 15, 15] = np.nan

    kwargs = dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=2)
    res_sentinel = vt.track(
        z_sentinel, np.array([15.0]), np.array([15.0]), t0=0, missing_value=-999.0, **kwargs
    )
    res_nan = vt.track(z_nan, np.array([15.0]), np.array([15.0]), t0=0, **kwargs)

    assert res_sentinel.status[0] == Status.TEMPLATE_READ_FAILED
    assert res_nan.status[0] == Status.TEMPLATE_READ_FAILED
    np.testing.assert_array_equal(res_sentinel.count, res_nan.count)
    np.testing.assert_array_equal(res_sentinel.status, res_nan.status)


def test_nan_does_not_mutate_input_array():
    z = _wave_field(nt=6, ny=20, nx=20)
    z[0, 10, 10] = np.nan
    z_before = z.copy()
    vt.track(
        z, np.array([10.0]), np.array([10.0]), t0=0,
        template=(5, 5), search_velocity=(1.0, 1.0), nsteps=1,
    )
    np.testing.assert_array_equal(z, z_before)
    assert np.isnan(z[0, 10, 10])


@pytest.mark.parametrize("seed_shape", [(), (5,), (3, 4)])
def test_seed_shape_variants(seed_shape):
    z = _wave_field()
    rng = np.random.RandomState(1)
    x0 = rng.uniform(15, 45, size=seed_shape)
    y0 = rng.uniform(15, 45, size=seed_shape)
    res = vt.track(z, x0, y0, t0=0, **_common_kwargs())

    nsteps = 4
    assert res.count.shape == seed_shape
    assert res.status.shape == seed_shape
    assert res.t_index.shape == (nsteps + 1, *seed_shape)
    assert res.x.shape == (nsteps + 1, *seed_shape)
    assert res.y.shape == (nsteps + 1, *seed_shape)
    assert res.vx.shape == (nsteps, *seed_shape)
    assert res.vy.shape == (nsteps, *seed_shape)
    assert res.score.shape == (nsteps, *seed_shape)


def test_tracker_and_track_agree():
    z = _wave_field()
    x0, y0 = vt.seed_grid(z.shape, spacing=5, margin=8)
    kwargs = _common_kwargs()

    tracker = vt.Tracker(**kwargs)
    res_tracker = tracker.track(z, x0, y0, t0=0)
    res_direct = vt.track(z, x0, y0, t0=0, **kwargs)

    np.testing.assert_array_equal(res_tracker.count, res_direct.count)
    np.testing.assert_array_equal(res_tracker.status, res_direct.status)
    for name in ("t_index", "x", "y", "vx", "vy", "score"):
        np.testing.assert_array_equal(getattr(res_tracker, name), getattr(res_direct, name))


def test_to_dataframe_shape_and_columns():
    pytest.importorskip("pandas")
    z = _wave_field()
    x0, y0 = vt.seed_grid(z.shape, spacing=8, margin=8)
    res = vt.track(z, x0, y0, t0=0, **_common_kwargs())
    df = res.to_dataframe()

    nsteps = 4
    n_seed = x0.size
    assert len(df) == n_seed * (nsteps + 1)
    for col in ("seed", "step", "count", "status", "t_index", "x", "y", "vx", "vy", "score"):
        assert col in df.columns
    # step==0 rows have no vx/vy/score yet (no step has been taken)
    assert df.loc[df["step"] == 0, "vx"].isna().all()


def test_to_xarray_shape_and_dims():
    pytest.importorskip("xarray")
    z = _wave_field()
    x0, y0 = vt.seed_grid(z.shape, spacing=8, margin=8)
    res = vt.track(z, x0, y0, t0=0, diagnostics=True, **_common_kwargs())
    ds = res.to_xarray()

    nsteps = 4
    assert ds["x"].shape == (nsteps + 1, *x0.shape)
    assert ds["vx"].shape == (nsteps, *x0.shape)
    assert "templates" in ds
    assert "score_grids" in ds


def test_to_xarray_missing_dependency_raises_informative_error(monkeypatch):
    import builtins

    z = _wave_field(nt=6, ny=20, nx=20)
    res = vt.track(
        z, np.array([10.0]), np.array([10.0]), t0=0,
        template=(5, 5), search_velocity=(1.0, 1.0), nsteps=1,
    )

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "xarray":
            raise ImportError("simulated: xarray not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="xarray"):
        res.to_xarray()


def test_to_dataframe_missing_dependency_raises_informative_error(monkeypatch):
    import builtins

    z = _wave_field(nt=6, ny=20, nx=20)
    res = vt.track(
        z, np.array([10.0]), np.array([10.0]), t0=0,
        template=(5, 5), search_velocity=(1.0, 1.0), nsteps=1,
    )

    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name == "pandas":
            raise ImportError("simulated: pandas not installed")
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="pandas"):
        res.to_dataframe()
