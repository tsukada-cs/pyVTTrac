"""Loose performance smoke test: catches gross regressions 
(e.g. an accidental O(n^2), or OpenMP silently not engaging)
without being sensitive to CI-runner variance. For the strict section-11
benchmark gate against VTTrac.jl, see tools/bench.py (run manually, not in CI).
"""

import time

import numpy as np

import pyvttrac as vt


def _field(nt=10, ny=120, nx=120):
    k = 2 * np.pi / 12
    cx = cy = 1.2
    t = np.arange(nt, dtype=np.float64)
    tg, yg, xg = np.meshgrid(t, np.arange(ny), np.arange(nx), indexing="ij")
    z = np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))
    return z.astype(np.float32)


def test_tracking_completes_within_generous_time_budget():
    z = _field()
    x0, y0 = vt.seed_grid(z.shape, spacing=4, margin=10)  # ~700 seeds

    t0 = time.perf_counter()
    res = vt.track(
        z, x0, y0, t0=0,
        template=(7, 7), search_velocity=(2.0, 2.0), nsteps=5, workers=1,
    )
    elapsed = time.perf_counter() - t0

    # Should take well under a second on any real machine; 15s leaves large
    # headroom for slow/shared CI runners while still catching a severe
    # regression (e.g. an accidentally-quadratic search).
    assert elapsed < 15.0, f"tracking took {elapsed:.2f}s, expected well under 15s"
    assert res.ok.mean() > 0.5  # sanity: the field/params should mostly track


def test_openmp_workers_do_not_slow_things_down():
    """Requesting more workers must never make things *slower* -- a cheap,
    non-flaky proxy for "OpenMP is wired up correctly" that doesn't assert a
    specific speedup ratio (shared/throttled CI runners can't guarantee one).
    """
    z = _field()
    x0, y0 = vt.seed_grid(z.shape, spacing=4, margin=10)
    kwargs = dict(template=(7, 7), search_velocity=(2.0, 2.0), nsteps=5)

    def best_of(workers, repeats=3):
        best = None
        for _ in range(repeats):
            t0 = time.perf_counter()
            vt.track(z, x0, y0, t0=0, workers=workers, **kwargs)
            dt = time.perf_counter() - t0
            best = dt if best is None else min(best, dt)
        return best

    t_serial = best_of(1)
    t_parallel = best_of(-1)  # all cores
    # Generous slack (2x) to absorb noise/oversubscription on shared CI hosts.
    assert t_parallel < t_serial * 2.0 + 0.05
