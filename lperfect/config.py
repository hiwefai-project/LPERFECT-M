# -*- coding: utf-8 -*-
"""Configuration handling for LPERFECT.

The model is configured via:
1) A JSON configuration file (config.json).
2) Optional CLI overrides (handled in cli.py).
"""

# Import JSON for reading configuration files.
import json

# Import typing primitives.
from typing import Any, Dict


def default_config() -> Dict[str, Any]:
    """Return a complete default configuration dictionary."""
    # Domain is NetCDF-only.
    return {
        "domain": {
            "mode": "netcdf",
            "domain_nc": "domain.nc",
            "varmap": {
                "dem": "dem",
                "d8": "d8",
                "cn": "cn",
                "channel_mask": "channel_mask",
                "x": "x",
                "y": "y",
            },
        },
        "model": {
            "start_time": "2025-12-21T00:00:00Z",
            "T_s": 7200,
            "dt_s": 5,
            "encoding": "esri",
            "ia_ratio": 0.2,
            "particle_vol_m3": 0.25,
            "travel_time_s": 5,
            "travel_time_channel_s": 1,
            "outflow_sink": True,
            "log_every": 10,
        },
        "rain": {
            "sources": {
                "radar": {
                    "kind": "netcdf",
                    "path": "radar_nowcast.nc",
                    "var": "rain_rate",
                    "time_var": "time",
                    "select": "nearest",
                    "mode": "intensity_mmph",
                    "weight": 0.6,
                },
                "station": {
                    "kind": "netcdf",
                    "path": "stations_nowcast.nc",
                    "var": "rain_rate",
                    "time_var": "time",
                    "select": "nearest",
                    "mode": "intensity_mmph",
                    "weight": 0.2,
                },
                "model": {
                    "kind": "netcdf",
                    "path": "wrf_forecast.nc",
                    "var": "rain_rate",
                    "time_var": "time",
                    "select": "nearest",
                    "mode": "intensity_mmph",
                    "weight": 0.2,
                },
            }
        },
        "risk": {
            "enabled": True,
            "balance": 0.55,
            "p_low": 5.0,
            "p_high": 95.0,
        },
        "restart": {
            "in": None,
            "out": "restart_state.nc",
            "every": 120,
            "strict_grid_check": True,
        },
        "output": {
            "out_netcdf": "flood_depth.nc",
            "Conventions": "CF-1.10",
            "title": "LPERFECT flood depth + hydrogeological risk index",
            "institution": "UniParthenope",
        },
        "compute": {
            "device": "cpu",
        },
    }


def load_json(path: str) -> Dict[str, Any]:
    """Load a JSON file into a Python dictionary."""
    # Open the file with UTF-8 encoding.
    with open(path, "r", encoding="utf-8") as f:
        # Parse JSON into Python dict.
        return json.load(f)


def deep_update(base: Dict[str, Any], other: Dict[str, Any]) -> Dict[str, Any]:
    """Recursively merge dict `other` into dict `base` (non-destructive)."""
    # Start from a shallow copy of base.
    out = dict(base)
    # Iterate keys from other.
    for k, v in other.items():
        # If both sides are dicts, merge recursively.
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_update(out[k], v)  # type: ignore[arg-type]
        else:
            # Otherwise override.
            out[k] = v
    # Return merged dictionary.
    return out
