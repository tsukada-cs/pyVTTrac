"""Performance gate.

Not part of the pytest suite (timing comparisons against another language's
implementation don't belong in CI-gating tests -- see tests/test_performance.py
for the loose, CI-safe smoke test). Run manually:

    python tools/bench.py

Expected (measured on the reference dev machine):
    Julia v2.0.0 (sequential):       0.2443 s  (6.605 us/(template*step))
    Fortran prototype (sequential):  0.2036 s  (5.505 us/(template*step))
    Fortran prototype (OpenMP x8):   --        (0.786 us/(template*step))

Gate: this package's `workers=1` time must be <= Julia's 0.2443 s, and
OpenMP must show a significant speedup as `workers` increases.
"""

import time

import numpy as np

import pyvttrac as vt


def build_field():
    nt, ny, nx = 12, 400, 400
    k = 2 * np.pi / 12
    cx = cy = 1.2
    t = np.arange(nt, dtype=np.float64)
    tg, yg, xg = np.meshgrid(t, np.arange(ny), np.arange(nx), indexing="ij")
    z = np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))
    return z.astype(np.float32)


def build_seeds():
    xs = np.arange(30, 371, 4, dtype=np.float64)
    ys = np.arange(30, 371, 4, dtype=np.float64)
    return np.meshgrid(xs, ys)


def run(z, x0, y0, nsteps, workers, repeats=5):
    best = None
    res = None
    for _ in range(repeats):
        t0 = time.perf_counter()
        res = vt.track(
            z, x0, y0, t0=0,
            template=(7, 7), search_radius=(3, 3), nsteps=nsteps,
            workers=workers,
        )
        dt = time.perf_counter() - t0
        best = dt if best is None else min(best, dt)
    return best, res


def main():
    z = build_field()
    x0, y0 = build_seeds()
    nsteps = 5
    n_points = x0.size

    print(f"field: {z.shape}, seeds: {n_points}, nsteps: {nsteps}")
    print(f"{'workers':>8}  {'time(s)':>9}  {'us/(tmpl*step)':>16}  {'valid':>7}  {'mean_vx':>8}")

    baseline = None
    for workers in (1, 2, 4, 8, -1):
        dt, res = run(z, x0, y0, nsteps, workers)
        if workers == 1:
            baseline = dt
        us_per = dt / (n_points * nsteps) * 1e6
        valid = int(res.count.sum())
        print(f"{workers!s:>8}  {dt:>9.4f}  {us_per:>16.3f}  {valid:>7}  {np.nanmean(res.vx):>8.4f}")

    julia_baseline_s = 0.2443
    print(f"\nsequential (workers=1): {baseline:.4f} s  (Julia v2.0.0 baseline: {julia_baseline_s} s)")
    if baseline <= julia_baseline_s:
        print(f"PASS: {baseline:.4f} s <= {julia_baseline_s} s")
    else:
        print(f"FAIL: {baseline:.4f} s >  {julia_baseline_s} s")


if __name__ == "__main__":
    main()
