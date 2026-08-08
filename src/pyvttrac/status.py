from enum import IntEnum


class Status(IntEnum):
    """Per-seed tracking outcome. Values match VTTrac.jl's STATUS_* constants."""

    OK = 0
    TID_START_OUT_OF_RANGE = 1
    TEMPLATE_READ_FAILED = 2
    LOW_CONTRAST = 3
    PEAK_NOT_INSIDE_TEMPLATE = 4
    TID_END_OUT_OF_RANGE = 5
    SCORE_COMPUTATION_FAILED = 6
    PEAK_NOT_FOUND = 7
    SCORE_BELOW_THRESHOLD = 8
    VELOCITY_CHANGE_TOO_LARGE = 9
