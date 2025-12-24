# -*- coding: utf-8 -*-
"""Time helpers for LPERFECT."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import datetime helpers.
from datetime import datetime, timezone  # import datetime import datetime, timezone

# Import numpy for datetime64 parsing.
import numpy as np  # import numpy as np

# Import CF schema constants.
from .cf_schema import RAIN_TIME_UNITS  # import .cf_schema import RAIN_TIME_UNITS


def utc_now_iso() -> str:  # define function utc_now_iso
    """Return current UTC time as an ISO-8601 string with 'Z'."""  # execute statement
    # Get current time in UTC.
    now = datetime.now(timezone.utc)  # set now
    # Convert to ISO string and force 'Z' suffix.
    return now.isoformat().replace("+00:00", "Z")  # return now.isoformat().replace("+00:00", "Z")


def parse_iso8601_to_datetime64(s: str | None) -> np.datetime64 | None:  # define function parse_iso8601_to_datetime64
    """Parse an ISO-8601 timestamp into numpy.datetime64.

    Notes
    -----
    - Accepts 'Z' suffix.
    - Returns None if s is None/empty.
    """
    # Handle null / empty.
    if not s:  # check condition not s:
        return None  # return None
    # Strip whitespace.
    ss = s.strip()  # set ss
    # Replace Z with explicit UTC offset for numpy.
    if ss.endswith("Z"):  # check condition ss.endswith("Z"):
        ss = ss[:-1] + "+00:00"  # set ss
    # Parse to aware datetime, convert to UTC, drop tzinfo (numpy has no tz support).
    dt = datetime.fromisoformat(ss)  # set dt
    dt_utc = dt.astimezone(timezone.utc).replace(tzinfo=None)  # set dt_utc
    # Convert to numpy datetime64 without emitting timezone warnings.
    return np.datetime64(dt_utc)  # return np.datetime64(dt_utc)


def datetime_to_hours_since_1900(dt: datetime) -> int:  # define function datetime_to_hours_since_1900
    """Convert datetime to integer hours since 1900-01-01 UTC."""  # execute statement
    base = datetime(1900, 1, 1, tzinfo=timezone.utc)  # set base
    hours = (dt - base).total_seconds() / 3600.0  # set hours
    return int(round(hours))  # return int(round(hours))


def rain_time_units() -> str:  # define function rain_time_units
    """Return the CF time units used by rain forcing."""  # execute statement
    return RAIN_TIME_UNITS  # return RAIN_TIME_UNITS
