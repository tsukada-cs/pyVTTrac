from __future__ import annotations

import numpy as np

from .result import TrackResult
from .status import Status


def concat_bidirectional(
    forward: TrackResult, backward: TrackResult, *, drop_shared_origin: bool = True
) -> TrackResult:
    """Combine a forward- and a backward-tracked `TrackResult` from a shared
    starting point into one trajectory spanning both directions.

    Only the shared-origin case is supported: `forward` and `backward` must
    have been tracked from the same `t0`/`x0`/`y0`, with `forward.step > 0`,
    `backward.step < 0`, and equal magnitude. Seeds whose forward and
    backward legs start at different times are outside this function's
    scope and must be combined by the caller.

    The combined position axis (`t_index`/`x`/`y`/`templates`) runs from
    `-len(backward.vx)` to `+len(forward.vx)` steps; the combined
    velocity/score axis (`vx`/`vy`/`score`/`score_grids`) is staggered by
    0.5 relative to that, same as within a single `TrackResult`. Velocities
    are **not** sign-flipped when reversing the backward leg -- `vx = (xw -
    xcur) / dt` already comes out as a correctly-signed physical velocity
    when `dt < 0` (backward tracking), so only the array order is reversed.

    `status` is `Status.OK` only when both legs completed OK; otherwise it's
    whichever leg's status is not OK (`forward`'s, if both failed). Each
    leg's own status is kept separately in `status_forward`/
    `status_backward`. The combined `count` is the number of valid
    (`t_index != -1`) points on the combined position axis.

    `grid`/`seed_x`/`seed_y`/`search_radius` are carried over only when
    `forward` and `backward` agree on them; otherwise `None`.

    Parameters
    ----------
    forward, backward : TrackResult
    drop_shared_origin : bool, default True
        The two legs' starting point is identical by construction. If
        `True` (default), it appears once in the combined result. If
        `False`, it appears twice (once from each leg) and the combined
        position axis is one element longer.
    """
    if not (forward.step > 0 and backward.step < 0):
        raise ValueError(
            "concat_bidirectional() requires forward.step > 0 and backward.step < 0 "
            f"(got forward.step={forward.step!r}, backward.step={backward.step!r})"
        )
    if abs(forward.step) != abs(backward.step):
        raise ValueError(
            "concat_bidirectional() requires |forward.step| == |backward.step| "
            f"(got {forward.step!r} and {backward.step!r})"
        )
    if forward.count.shape != backward.count.shape:
        raise ValueError(
            "forward and backward must have the same seed shape "
            f"(got {forward.count.shape!r} and {backward.count.shape!r})"
        )
    if not np.array_equal(forward.t_index[0], backward.t_index[0]):
        raise ValueError(
            "concat_bidirectional() only supports the shared-origin case: "
            "forward.t_index[0] and backward.t_index[0] must be identical "
            "(both legs tracked from the same t0/x0/y0). Seeds whose "
            "forward/backward legs start at different times must be "
            "combined by the caller."
        )

    has_templates = forward.templates is not None
    if has_templates != (backward.templates is not None):
        raise ValueError("forward and backward must either both have `templates` or neither")
    if has_templates and forward.templates.shape[1:] != backward.templates.shape[1:]:
        raise ValueError(
            "forward.templates and backward.templates have incompatible shapes "
            f"(got {forward.templates.shape!r} and {backward.templates.shape!r})"
        )

    has_score_grids = forward.score_grids is not None
    if has_score_grids != (backward.score_grids is not None):
        raise ValueError("forward and backward must either both have `score_grids` or neither")
    if has_score_grids and forward.score_grids.shape[1:] != backward.score_grids.shape[1:]:
        raise ValueError(
            "forward.score_grids and backward.score_grids have incompatible shapes "
            f"(got {forward.score_grids.shape!r} and {backward.score_grids.shape!r})"
        )

    def concat_position(fwd_arr, bwd_arr):
        parts = [bwd_arr[1:][::-1]]
        if not drop_shared_origin:
            parts.append(bwd_arr[0:1])
        parts.append(fwd_arr)
        return np.concatenate(parts, axis=0)

    def concat_velocity(fwd_arr, bwd_arr):
        # No sign flip: see the docstring note on backward velocity sign.
        return np.concatenate([bwd_arr[::-1], fwd_arr], axis=0)

    t_index = concat_position(forward.t_index, backward.t_index)
    x = concat_position(forward.x, backward.x)
    y = concat_position(forward.y, backward.y)
    vx = concat_velocity(forward.vx, backward.vx)
    vy = concat_velocity(forward.vy, backward.vy)
    score = concat_velocity(forward.score, backward.score)
    templates = concat_position(forward.templates, backward.templates) if has_templates else None
    score_grids = (
        concat_velocity(forward.score_grids, backward.score_grids) if has_score_grids else None
    )

    count = np.sum(t_index != -1, axis=0)

    both_ok = (forward.status == Status.OK) & (backward.status == Status.OK)
    status = np.where(
        both_ok, Status.OK, np.where(forward.status != Status.OK, forward.status, backward.status)
    )

    grid = forward.grid if forward.grid == backward.grid else None
    seed_x = (
        forward.seed_x
        if forward.seed_x is not None
        and backward.seed_x is not None
        and np.array_equal(forward.seed_x, backward.seed_x)
        else None
    )
    seed_y = (
        forward.seed_y
        if forward.seed_y is not None
        and backward.seed_y is not None
        and np.array_equal(forward.seed_y, backward.seed_y)
        else None
    )
    search_radius = (
        forward.search_radius if forward.search_radius == backward.search_radius else None
    )

    return TrackResult(
        count=count,
        status=status,
        t_index=t_index,
        x=x,
        y=y,
        vx=vx,
        vy=vy,
        score=score,
        templates=templates,
        score_grids=score_grids,
        step=forward.step,
        step_offset=backward.vx.shape[0],
        grid=grid,
        seed_x=seed_x,
        seed_y=seed_y,
        search_radius=search_radius,
        status_forward=forward.status,
        status_backward=backward.status,
    )
