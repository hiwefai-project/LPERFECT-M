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


def parse_iso8601_to_utc_datetime(value: str | None) -> datetime:  # define function parse_iso8601_to_utc_datetime
    """Parse ISO-8601 string into an aware UTC datetime (fallback: now)."""  # execute statement
    if not value:  # check condition not value:
        return datetime.now(timezone.utc)  # return datetime.now(timezone.utc)
    text = value.strip()  # set text
    if text.endswith("Z"):  # check condition text.endswith("Z"):
        text = text[:-1] + "+00:00"  # set text
    dt = datetime.fromisoformat(text)  # set dt
    if dt.tzinfo is None:  # check condition dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)  # set dt
    return dt.astimezone(timezone.utc)  # return dt.astimezone(timezone.utc)


def datetime_to_hours_since_1900(dt: datetime) -> float:  # define function datetime_to_hours_since_1900
    """Convert datetime to hours since 1900-01-01 UTC."""  # execute statement
    base = datetime(1900, 1, 1, tzinfo=timezone.utc)  # set base
    hours = (dt - base).total_seconds() / 3600.0  # set hours
    return float(hours)  # return float(hours)


def rain_time_units() -> str:  # define function rain_time_units
    """Return the CF time units used by rain forcing."""  # execute statement
    return RAIN_TIME_UNITS  # return RAIN_TIME_UNITS
