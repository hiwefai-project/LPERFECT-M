# -*- coding: utf-8 -*-
"""Configuration handling for LPERFECT.

The model is configured via:
1) A JSON configuration file (config.json).
2) Optional CLI overrides (handled in cli.py).
"""

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import JSON for reading configuration files.
import json  # import json

# Import typing primitives.
from typing import Any, Dict  # import typing import Any, Dict

# Import CF schema defaults.
from .cf_schema import (  # import .cf_schema
    CF_CONVENTIONS,
    RAIN_CRS_VAR,
    RAIN_LAT_VAR,
    RAIN_LON_VAR,
    RAIN_RATE_UNITS,
    RAIN_RATE_VAR,
    RAIN_TIME_UNITS,
    RAIN_TIME_VAR,
)


def default_config() -> Dict[str, Any]:  # define function default_config
    """Return a complete default configuration dictionary."""  # execute statement
    # Domain is NetCDF-only.
    return {  # return {
        "domain": {  # execute statement
            "mode": "netcdf",  # execute statement
            "domain_nc": "domain.nc",  # execute statement
            "varmap": {  # execute statement
                "dem": "dem",  # execute statement
                "d8": "d8",  # execute statement
                "cn": "cn",  # execute statement
                "channel_mask": "channel_mask",  # execute statement
                "x": "longitude",  # execute statement
                "y": "latitude",  # execute statement
            },  # execute statement
        },  # execute statement
        "model": {  # execute statement
            "start_time": "2025-12-21T00:00:00Z",  # execute statement
            "T_s": 7200,  # execute statement
            "dt_s": 5,  # execute statement
            "encoding": "esri",  # execute statement
            "ia_ratio": 0.2,  # execute statement
            "particle_vol_m3": 0.25,  # execute statement
            "travel_time_s": 5,  # execute statement
            "travel_time_channel_s": 1,  # execute statement
            "travel_time_mode": "fixed",  # execute statement
            "travel_time_auto": {  # execute statement
                "hillslope_velocity_ms": 0.5,  # execute statement
                "channel_velocity_ms": 1.5,  # execute statement
                "min_s": 0.25,  # execute statement
                "max_s": 3600.0,  # execute statement
            },  # execute statement
            "outflow_sink": True,  # execute statement
            "log_every": 10,  # execute statement
        },  # execute statement
        "rain": {  # execute statement
            "schema": {  # execute statement
                "time_var": RAIN_TIME_VAR,  # execute statement
                "lat_var": RAIN_LAT_VAR,  # execute statement
                "lon_var": RAIN_LON_VAR,  # execute statement
                "rain_var": RAIN_RATE_VAR,  # execute statement
                "crs_var": RAIN_CRS_VAR,  # execute statement
                "time_units": RAIN_TIME_UNITS,  # execute statement
                "rate_units": RAIN_RATE_UNITS,  # execute statement
                "require_cf": True,  # execute statement
                "require_time_dim": True,  # execute statement
            },  # execute statement
            "sources": {  # execute statement
                "rain": {  # execute statement
                    "kind": "netcdf",  # execute statement
                    "path": "rain_time_dependent.nc",  # execute statement
                    "var": RAIN_RATE_VAR,  # execute statement
                "time_var": RAIN_TIME_VAR,  # execute statement
                "select": "previous",  # execute statement
                "mode": "intensity_mmph",  # execute statement
                "weight": 1.0,  # execute statement
            },  # execute statement
            }  # execute statement
        },  # execute statement
        "risk": {  # execute statement
            "enabled": True,  # execute statement
            "balance": 0.55,  # execute statement
            "p_low": 5.0,  # execute statement
            "p_high": 95.0,  # execute statement
        },  # execute statement
        "restart": {  # execute statement
            "in": None,  # execute statement
            "out": "restart_state.nc",  # execute statement
            "every": 120,  # execute statement
            "strict_grid_check": True,  # execute statement
        },  # execute statement
        "output": {  # execute statement
            "out_netcdf": "flood_depth.nc",  # execute statement
            "Conventions": CF_CONVENTIONS,  # execute statement
            "title": "LPERFECT flood depth + hydrogeological risk index",  # execute statement
            "institution": "UniParthenope",  # execute statement
        },  # execute statement
        "compute": {  # execute statement
            "device": "cpu",  # execute statement
            "shared_memory": {  # execute statement
                "enabled": False,  # execute statement
                "workers": None,  # execute statement
                "min_particles_per_worker": 5000,  # execute statement
                "chunk_size": 65536,  # execute statement
            },  # execute statement
        },  # execute statement
    }  # execute statement


def load_json(path: str) -> Dict[str, Any]:  # define function load_json
    """Load a JSON file into a Python dictionary."""  # execute statement
    # Open the file with UTF-8 encoding.
    with open(path, "r", encoding="utf-8") as f:  # manage context open(path, "r", encoding="utf-8") as f:
        # Parse JSON into Python dict.
        return json.load(f)  # return json.load(f)


def deep_update(base: Dict[str, Any], other: Dict[str, Any]) -> Dict[str, Any]:  # define function deep_update
    """Recursively merge dict `other` into dict `base` (non-destructive)."""  # execute statement
    # Start from a shallow copy of base.
    out = dict(base)  # set out
    # Iterate keys from other.
    for k, v in other.items():  # loop over k, v in other.items():
        # If both sides are dicts, merge recursively.
        if isinstance(v, dict) and isinstance(out.get(k), dict):  # check condition isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)  # type: ignore[arg-type]  # execute statement
        else:  # fallback branch
            # Otherwise override.
            out[k] = v  # execute statement
    # Return merged dictionary.
    return out  # return out
