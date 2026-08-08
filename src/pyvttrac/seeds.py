from __future__ import annotations

import numpy as np


def seed_grid(shape, spacing, margin=0):
    """A regular grid of seed points inside a (..., ny, nx) image, in index units.

    Parameters
    ----------
    shape : tuple
        Shape of the data array to seed, e.g. `z.shape` (only the last two
        dimensions, ny and nx, are used).
    spacing : float or (float, float)
        Spacing between seed points, as a scalar (same in x and y) or (sy, sx).
    margin : float or (float, float), default 0
        Margin to leave empty at the image edges, as a scalar or (my, mx).

    Returns
    -------
    x0, y0 : np.ndarray
        2-D arrays of seed x/y positions (index units), suitable for `track()`.
    """
    ny, nx = shape[-2], shape[-1]
    if np.isscalar(spacing):
        sy = sx = float(spacing)
    else:
        sy, sx = spacing
    if np.isscalar(margin):
        my = mx = float(margin)
    else:
        my, mx = margin

    xs = np.arange(mx, nx - mx, sx, dtype=np.float64)
    ys = np.arange(my, ny - my, sy, dtype=np.float64)
    x0, y0 = np.meshgrid(xs, ys)
    return x0, y0
