# Quickstart

This walks through [`examples/quickstart.py`](../examples/quickstart.py),
a self-contained script over synthetic data (no external files needed).
Run it directly, or step through it cell-by-cell (the `#%%` markers are
[Jupyter-style cells](https://code.visualstudio.com/docs/python/jupyter-support-py))
in VS Code or Spyder.

```bash
pip install pyVTTrac matplotlib
python examples/quickstart.py
```

## 1. Get some image-like data

pyVTTrac tracks features through a `(nt, ny, nx)` sequence of 2-D arrays —
satellite imagery, PIV frames, model output, anything where a coherent
pattern translates over time. The quickstart script synthesizes a traveling
wave so the "true" velocity is known in advance (`cx = cy = 1.2` pixels per
time step):

```python
import numpy as np

nt, ny, nx = 20, 80, 100
t = np.arange(nt, dtype=np.float64)
tg, yg, xg = t[:, None, None], np.arange(ny)[None, :, None], np.arange(nx)[None, None, :]
k = 2 * np.pi / 12
cx, cy = 1.2, 1.2
z = (np.sin(k * (xg - cx * tg)) * np.cos(k * (yg - cy * tg))).astype(np.float32)
```

`z` must have `dtype` convertible to `float32` and at least 2 time steps.

## 2. Choose seed points

`seed_grid()` lays out a regular grid of points to track, in index units
(pixels), leaving a margin so the template never starts off the edge of the
image:

```python
import pyvttrac as vt

x0, y0 = vt.seed_grid(z.shape, spacing=(1.0, 2.5), margin=(10.5, 7.5))
```

You can just as easily pass your own arbitrary-shaped `x0`/`y0` arrays
(a 1-D list of points of interest, a 2-D grid, whatever fits your problem);
`track()`'s output arrays follow whatever shape you give it.

## 3. Track

```python
res = vt.track(
    z, x0, y0, t0=0, time=t,
    template=(5, 5),              # (ny, nx) template size
    search_velocity=(1.8, 1.8),   # (vy, vx) search range, pixels/time-unit
    nsteps=nt - 1,
    diagnostics=True,             # keep intermediate templates/score maps
)
```

`template` should comfortably contain the feature you're tracking;
`search_velocity` should comfortably bound the true displacement between
steps (or use `search_radius=(iy, ix)` to specify the search window directly,
in pixels, instead of a velocity).

## 4. Read the result

`res` is a `TrackResult`:

```python
res.ok            # bool array, seed_shape: which seeds tracked to completion
res.status        # Status enum values, seed_shape: why a seed stopped, if it did
res.x, res.y      # (nsteps+1, *seed_shape): tracked positions (NaN once invalid)
res.vx, res.vy    # (nsteps,   *seed_shape): velocities between steps
res.score         # (nsteps,   *seed_shape): match quality at each step
res.t_index       # (nsteps+1, *seed_shape): time index at each step (-1 once invalid)
```

```python
print(f"tracked {res.ok.sum()} / {res.ok.size} points to completion")
print(f"mean vx = {np.nanmean(res.vx):.3f}  (true: {cx})")
```

Need a labeled `xarray.Dataset` or a tidy long-format table instead of raw
arrays?

```python
ds = res.to_xarray()      # pip install "pyVTTrac[xarray]"
df = res.to_dataframe()   # pip install "pyVTTrac[pandas]"
```

## Where to go next

- [`docs/api.md`](api.md) — full parameter reference and the `Status` code
  table (why a given seed's tracking stopped).
- [`docs/migration-v1-to-v2.md`](migration-v1-to-v2.md) — coming from v1?
