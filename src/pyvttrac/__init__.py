from importlib.metadata import PackageNotFoundError, version as _version

from .api import Tracker, search_radius_from_velocity, track, velocity_from_search_radius
from .bidirectional import concat_bidirectional
from .grid import Grid
from .result import TrackResult
from .seeds import seed_grid
from .status import Status

try:
    # Single source of truth is meson.build's `project(..., version: ...)`;
    # meson-python exposes it as this distribution's installed metadata (see
    # `dynamic = ["version"]` in pyproject.toml), so it's read back here
    # rather than duplicated as a literal.
    __version__ = _version("pyVTTrac")
except PackageNotFoundError:
    # Not installed (e.g. running from a source checkout without `pip
    # install -e .`) -- there's no metadata to read yet.
    __version__ = "unknown"

__all__ = [
    "track",
    "Tracker",
    "Grid",
    "TrackResult",
    "Status",
    "seed_grid",
    "search_radius_from_velocity",
    "velocity_from_search_radius",
    "concat_bidirectional",
    "__version__",
]
