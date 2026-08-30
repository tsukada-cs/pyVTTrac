"""min_score_margin: ridge screening on the score surface.

On quasi-1-D texture the correlation surface is a ridge running to the edge
of the search window, so the along-ridge displacement is ill-determined even
though the peak score is high. The margin check rejects exactly that case.
"""
import numpy as np
import pytest

import pyvttrac as vt
from pyvttrac import Status


def _make_sequence(pattern, vx, nt=4):
    """Translate `pattern` by vx pixels per step (frozen field, wrap-free interior)."""
    ny, nx = pattern.shape
    z = np.empty((nt, ny, nx), dtype=np.float32)
    xx = np.arange(nx)
    for t in range(nt):
        for iy in range(ny):
            z[t, iy] = np.interp(xx - vx * t, xx, pattern[iy])
    return z


def _seeds():
    x0, y0 = np.meshgrid(np.arange(24, 40, 4, dtype=float), np.arange(24, 40, 4, dtype=float))
    return x0, y0


@pytest.fixture()
def streak_field():
    """Nearly 1-D texture: strong variation across x, a whisper across y, so
    the integer peak sits inside the window but the ridge reaches its edge."""
    ny = nx = 64
    yy, xx = np.mgrid[0:ny, 0:nx].astype(np.float32)
    pattern = np.sin(2 * np.pi * xx / 10.0) + 0.02 * np.sin(2 * np.pi * yy / 40.0)
    return _make_sequence(pattern, vx=1.0)


@pytest.fixture()
def isotropic_field():
    """Smoothed white noise: isotropic texture whose autocorrelation decays
    well inside the search window, i.e. the score peak is compact."""
    rng = np.random.default_rng(0)
    ny = nx = 64
    noise = rng.normal(0.0, 1.0, (ny, nx)).astype(np.float32)
    kernel = np.exp(-0.5 * (np.arange(-4, 5) / 1.2) ** 2).astype(np.float32)
    pattern = np.apply_along_axis(lambda r: np.convolve(r, kernel, "same"), 1, noise)
    pattern = np.apply_along_axis(lambda c: np.convolve(c, kernel, "same"), 0, pattern)
    return _make_sequence(pattern, vx=1.0)


def test_streaks_rejected_with_margin(streak_field):
    x0, y0 = _seeds()
    kwargs = dict(template=(9, 9), search_radius=(3, 3), nsteps=2,
                  min_score=(0.5, 0.5), subgrid="paraboloid")
    plain = vt.track(streak_field, x0, y0, t0=0, **kwargs)
    screened = vt.track(streak_field, x0, y0, t0=0, min_score_margin=0.05, **kwargs)

    # without the margin the ridge tracks "successfully"
    assert np.isfinite(plain.vx).all()
    # with it, every seed is rejected for exactly the ridge reason
    assert (screened.status == Status.SCORE_MARGIN_TOO_SMALL).all()
    assert np.isnan(screened.vx).all()


def test_isotropic_texture_unaffected_by_margin(isotropic_field):
    x0, y0 = _seeds()
    kwargs = dict(template=(9, 9), search_radius=(3, 3), nsteps=2,
                  min_score=(0.5, 0.5), subgrid="paraboloid")
    plain = vt.track(isotropic_field, x0, y0, t0=0, **kwargs)
    screened = vt.track(isotropic_field, x0, y0, t0=0, min_score_margin=0.05, **kwargs)

    assert (plain.status == Status.OK).all()
    np.testing.assert_array_equal(screened.status, plain.status)
    np.testing.assert_allclose(screened.vx, plain.vx, rtol=0, atol=0)
    # linear-interpolation advection plus subgrid refinement is only good to
    # ~0.1 px/step here; the point is that the field tracks sensibly at all
    np.testing.assert_allclose(np.asarray(plain.vx)[np.isfinite(plain.vx)], 1.0, atol=0.15)


def test_margin_disabled_matches_default(isotropic_field, streak_field):
    x0, y0 = _seeds()
    kwargs = dict(template=(9, 9), search_radius=(3, 3), nsteps=2, min_score=(0.5, 0.5))
    for z in (isotropic_field, streak_field):
        a = vt.track(z, x0, y0, t0=0, **kwargs)
        b = vt.track(z, x0, y0, t0=0, min_score_margin=None, **kwargs)
        np.testing.assert_array_equal(a.status, b.status)
        np.testing.assert_allclose(a.vx, b.vx, equal_nan=True)


def test_negative_margin_rejected():
    with pytest.raises(ValueError, match="min_score_margin"):
        vt.Tracker(template=(9, 9), search_radius=(3, 3), min_score_margin=-0.1)
