"""Input validation: every check must raise a clear
exception rather than crash, corrupt memory, or silently misbehave.
"""

import numpy as np
import pytest

import pyvttrac as vt
from pyvttrac.config import TrackingConfig
from pyvttrac.grid import Grid
from pyvttrac.status import Status


def _small_z(nt=6, ny=20, nx=20):
    rng = np.random.RandomState(0)
    return rng.randn(nt, ny, nx).astype(np.float32)


# ---------------------------------------------------------------------------
# TrackingConfig / Tracker construction-time validation
# ---------------------------------------------------------------------------


def test_template_size_must_be_positive():
    with pytest.raises(ValueError):
        TrackingConfig(template=(0, 5), search_velocity=(1.0, 1.0))
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, -1), search_velocity=(1.0, 1.0))


def test_min_peak_prominence_requires_template_at_least_3():
    with pytest.raises(ValueError):
        TrackingConfig(template=(2, 2), search_velocity=(1.0, 1.0), min_peak_prominence=0.03)
    # must not raise when template is big enough
    TrackingConfig(template=(3, 3), search_velocity=(1.0, 1.0), min_peak_prominence=0.03)


def test_method_must_be_xcor_or_ncov():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), method="bogus")


def test_subgrid_must_be_valid_choice():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), subgrid="bogus")
    # valid choices must not raise
    for s in ("paraboloid", "gaussian", None):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), subgrid=s)


def test_search_velocity_and_search_radius_are_mutually_exclusive_and_required():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5))  # neither given
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), search_radius=(3, 3))  # both given


def test_step_must_not_be_zero():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), step=0)


def test_nsteps_must_be_at_least_1():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), nsteps=0)


def test_min_samples_must_be_at_least_1():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), min_samples=0)


def test_max_velocity_change_components_must_be_positive():
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), max_velocity_change=(0.0, 1.0))
    with pytest.raises(ValueError):
        TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), max_velocity_change=(1.0, -0.5))


def test_min_score_scalar_is_broadcast_to_pair():
    cfg = TrackingConfig(template=(5, 5), search_velocity=(1.0, 1.0), min_score=0.75)
    assert cfg.min_score == (0.75, 0.75)


# ---------------------------------------------------------------------------
# track() / Tracker.track() data-shape validation
# ---------------------------------------------------------------------------


def test_z_must_be_3d():
    z2d = _small_z()[0]
    with pytest.raises(ValueError):
        vt.track(z2d, np.array([5.0]), np.array([5.0]), template=(5, 5), search_velocity=(1.0, 1.0))


def test_z_must_have_at_least_2_time_steps():
    z1 = _small_z(nt=1)
    with pytest.raises(ValueError):
        vt.track(z1, np.array([5.0]), np.array([5.0]), template=(5, 5), search_velocity=(1.0, 1.0))


def test_mask_shape_must_match_z():
    z = _small_z()
    bad_mask = np.zeros((z.shape[0] - 1, *z.shape[1:]), dtype=bool)
    with pytest.raises(ValueError):
        vt.track(
            z, np.array([5.0]), np.array([5.0]), template=(5, 5), search_velocity=(1.0, 1.0),
            mask=bad_mask,
        )


def test_x0_y0_shape_mismatch():
    z = _small_z()
    with pytest.raises(ValueError):
        vt.track(
            z, np.array([5.0, 6.0]), np.array([5.0]),
            template=(5, 5), search_velocity=(1.0, 1.0),
        )


def test_t0_shape_mismatch():
    z = _small_z()
    with pytest.raises(ValueError):
        vt.track(
            z, np.array([5.0, 6.0]), np.array([5.0, 6.0]), t0=np.array([0, 1, 2]),
            template=(5, 5), search_velocity=(1.0, 1.0),
        )


def test_time_shape_mismatch():
    z = _small_z()
    with pytest.raises(ValueError):
        vt.track(
            z, np.array([5.0]), np.array([5.0]), time=np.arange(z.shape[0] - 1, dtype=np.float64),
            template=(5, 5), search_velocity=(1.0, 1.0),
        )


# ---------------------------------------------------------------------------
# Grid validation
# ---------------------------------------------------------------------------


def test_grid_dx_dy_must_be_nonzero():
    with pytest.raises(ValueError):
        Grid(x0=0.0, y0=0.0, dx=0.0, dy=1.0)
    with pytest.raises(ValueError):
        Grid(x0=0.0, y0=0.0, dx=1.0, dy=0.0)


# ---------------------------------------------------------------------------
# Template-read boundary handling: a seed whose subgrid-interpolation stencil
# leans toward the domain's low edge (fractional part < 0, so the bilinear
# stencil needs one extra column/row *below* the template's own footprint)
# must be rejected like any other out-of-range template, not read past the
# array's start.
# ---------------------------------------------------------------------------


def _edge_marker_z(nt=4, ny=20, nx=20):
    z = np.zeros((nt, ny, nx), dtype=np.float32)
    z[:, :, -1] = 1000.0  # distinctive value at the far (right/bottom) edge
    return z


def test_low_edge_subpixel_seed_x_is_rejected_not_corrupted():
    z = _edge_marker_z()
    # 0-based x=1.9 -> fractional part < 0 after rounding -> stencil reaches
    # for the column left of the template's own left edge.
    res = vt.track(
        z, np.array([1.9]), np.array([10.0]), t0=0,
        template=(5, 5), search_radius=(2, 2), nsteps=1,
        first_guess=(np.array([0.0]), np.array([5.0])),
        min_score=(-2.0, -2.0), diagnostics="templates",
    )
    assert res.status[0] == Status.TEMPLATE_READ_FAILED
    # The far-edge marker must never appear in this template: if it does,
    # the read wrapped around the array instead of being rejected.
    assert not np.any(res.templates[0, :, :, 0] == 1000.0)


def test_low_edge_subpixel_seed_y_is_rejected_not_corrupted():
    z = _edge_marker_z()
    res = vt.track(
        z, np.array([10.0]), np.array([1.9]), t0=0,
        template=(5, 5), search_radius=(2, 2), nsteps=1,
        first_guess=(np.array([5.0]), np.array([0.0])),
        min_score=(-2.0, -2.0), diagnostics="templates",
    )
    assert res.status[0] == Status.TEMPLATE_READ_FAILED
    assert not np.any(res.templates[0, :, :, 0] == 1000.0)


def test_low_edge_subpixel_seed_is_rejected_with_mask_too():
    # read_template_visible() has the same stencil-vs-bounds-check shape as
    # read_template(); exercise it via a mask that doesn't otherwise reject.
    z = _edge_marker_z()
    mask = np.zeros(z.shape, dtype=bool)
    res = vt.track(
        z, np.array([1.9]), np.array([10.0]), t0=0,
        template=(5, 5), search_radius=(2, 2), nsteps=1,
        first_guess=(np.array([0.0]), np.array([5.0])),
        min_score=(-2.0, -2.0), mask=mask, min_samples=1,
    )
    assert res.status[0] == Status.TEMPLATE_READ_FAILED


def test_integer_seed_at_same_edge_is_unaffected():
    # Sanity check: an integer-valued seed at the same nominal position
    # (no subpixel stencil, isx == 0) is a legitimate in-bounds read and
    # must keep succeeding -- the boundary fix must not overcorrect.
    z = _edge_marker_z()
    res = vt.track(
        z, np.array([2.0]), np.array([10.0]), t0=0,
        template=(5, 5), search_radius=(2, 2), nsteps=1,
        first_guess=(np.array([0.0]), np.array([5.0])),
        min_score=(-2.0, -2.0),
    )
    assert res.status[0] != Status.TEMPLATE_READ_FAILED
