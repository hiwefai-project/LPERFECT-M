#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_wrf_rain_to_rain_time_dependent.py

Convert a WRF-derived NetCDF (like wrf5_d02_*.nc) into a NetCDF compliant with a
simple time-dependent rain-rate CDL (cdl/rain_time_dependent.cdl).

Your source file (per the provided CDL) has:
- dimensions: time, latitude, longitude
- coordinates: time(time), latitude(latitude), longitude(longitude)
- rain accumulation for the hour in: RAIN_DELTA or DELTA_RAIN (mm)

We produce:
- rain_rate(time, latitude, longitude) as *liquid water equivalent precipitation rate*
  in **kg m-2 s-1**, which is CF-friendly (1 mm = 1 kg m-2).

Conversion:
- If DELTA_RAIN is "hourly cumulated rain" in mm over Δt = 1 hour:
    rain_rate_kg_m2_s = (DELTA_RAIN_mm / 3600.0)
- If you know your accumulation period is different, pass --accum-hours.

The output keeps:
- time/lat/lon coordinates
- FillValue handling (1e37 -> NaN -> FillValue)
- minimal CF-1.10 global attributes

Usage example:
  python convert_wrf_rain_to_rain_time_dependent.py \
    --in wrf5_d02_20251221Z1200.nc \
    --out rain_time_dependent_20251221Z1200.nc \
    --rain-var RAIN_DELTA

Notes:
- Your sample CDL shows DELTA_RAIN, while your message says RAIN_DELTA.
  This script supports BOTH; it auto-detects if --rain-var is not provided.
"""

from __future__ import annotations  # pedagogical: enables forward refs in type hints (py3.11+ ok)

import argparse  # pedagogical: standard CLI parsing
import logging   # pedagogical: structured logs for reproducible pipelines
from typing import Optional, Any  # pedagogical: type hints for readability

import numpy as np               # pedagogical: numerical arrays and NaN handling
import xarray as xr              # pedagogical: NetCDF I/O and labeled data


# pedagogical: module-level logger (configured in main())
LOG = logging.getLogger("wrf_rain_converter")


def _to_float(x: Any) -> Optional[float]:
    """Safely convert a scalar to float; return None on failure."""
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _infer_coord_name(ds: xr.Dataset, candidates: list[str]) -> Optional[str]:
    """
    Find a coordinate name in a dataset (works for both ds.coords and ds.variables).
    This makes the script robust to variants like lat/latitude/y and lon/longitude/x.
    """
    for c in candidates:
        if c in ds.coords or c in ds.variables:
            return c
    return None


def _get_fill_value(da: xr.DataArray) -> Optional[float]:
    """Read _FillValue or missing_value attribute if present."""
    for k in ("_FillValue", "missing_value"):
        if k in da.attrs:
            v = _to_float(da.attrs.get(k))
            if v is not None:
                return v
    return None


def _mask_fill(da: xr.DataArray) -> xr.DataArray:
    """
    Replace FillValue with NaN and ensure float dtype.
    This simplifies conversion math and avoids plotting/processing artifacts.
    """
    out = da.astype(float)
    fill = _get_fill_value(out)
    if fill is not None:
        out = out.where(out != fill)  # pedagogical: xarray 'where' keeps coords/dims intact
    return out


def convert(
    in_path: str,
    out_path: str,
    rain_var: Optional[str],
    accum_hours: float,
) -> None:
    """
    Perform the conversion and write the output NetCDF.
    """
    # 1) Open input dataset (xarray will pick netCDF4/h5netcdf backend automatically)
    LOG.info("Opening input: %s", in_path)
    ds = xr.open_dataset(in_path)

    # 2) Infer coordinate names (or rely on typical names)
    lat_name = _infer_coord_name(ds, ["latitude", "lat", "y"])
    lon_name = _infer_coord_name(ds, ["longitude", "lon", "x"])
    time_name = _infer_coord_name(ds, ["time"])

    if not lat_name or not lon_name or not time_name:
        raise ValueError(
            f"Could not infer coords. Found lat={lat_name}, lon={lon_name}, time={time_name}. "
            "Use standard names latitude/longitude/time in the input."
        )

    # 3) Choose the rain accumulation variable
    #    - user may pass --rain-var
    #    - else auto-detect common candidates
    if rain_var is None:
        for cand in ["RAIN_DELTA", "DELTA_RAIN", "HOURLY_RAIN", "RAIN"]:
            if cand in ds.variables:
                rain_var = cand
                break

    if rain_var is None or rain_var not in ds.variables:
        raise KeyError(
            f"Rain variable not found. Pass --rain-var. Available vars: {list(ds.variables)}"
        )

    LOG.info("Using rain accumulation variable: %s", rain_var)

    # 4) Extract the accumulation field and clean FillValues
    rain_acc_mm = _mask_fill(ds[rain_var])  # units mm (per your CDL)

    # 5) Make sure dimensions match (time, lat, lon)
    #    We don't hard-fail on order; we just ensure the dims exist.
    for req in (time_name, lat_name, lon_name):
        if req not in rain_acc_mm.dims:
            raise ValueError(
                f"Rain variable dims are {rain_acc_mm.dims}, expected to include {req}."
            )

    # 6) Convert accumulated mm over accum_hours to rate in kg m-2 s-1
    #    Pedagogical:
    #    - 1 mm of water over 1 m^2 equals 1 kg of water (density ~ 1000 kg/m^3).
    #    - So: mm / hour -> kg m-2 s-1 by dividing by (hours * 3600).
    seconds = float(accum_hours) * 3600.0
    rain_rate = rain_acc_mm / seconds

    # 7) Prepare output variable attributes consistent with CF conventions
    rain_rate.attrs = {
        "standard_name": "lwe_precipitation_rate",  # CF-friendly for liquid water equivalent precipitation rate
        "long_name": "rain_rate",
        "description": f"Rain rate derived from {rain_var} (accumulated over {accum_hours} h)",
        "units": "kg m-2 s-1",
        "source_variable": rain_var,
        "accumulation_period_hours": accum_hours,
    }

    # 8) Build output dataset with required coordinates and a single data variable
    out = xr.Dataset(
        data_vars={
            "rain_rate": rain_rate
        },
        coords={
            time_name: ds[time_name],
            lat_name: ds[lat_name],
            lon_name: ds[lon_name],
        },
        attrs={
            "Conventions": "CF-1.10",
            "title": "Time-dependent rain rate (converted from WRF-derived accumulation)",
            "source": "WRF post-processed output",
            "history": "Converted to rain_time_dependent format by convert_wrf_rain_to_rain_time_dependent.py",
        },
    )

    # 9) Normalize coordinate names to the target CDL expectations if needed
    #    Your other CDL uses exactly: time, latitude, longitude.
    rename_map = {}
    if time_name != "time":
        rename_map[time_name] = "time"
    if lat_name != "latitude":
        rename_map[lat_name] = "latitude"
    if lon_name != "longitude":
        rename_map[lon_name] = "longitude"
    if rename_map:
        out = out.rename(rename_map)

    # 10) Copy coordinate attributes (units/long_name) if present, else set defaults
    out["time"].attrs = dict(ds[time_name].attrs) if time_name in ds else {}
    out["latitude"].attrs = dict(ds[lat_name].attrs) if lat_name in ds else {}
    out["longitude"].attrs = dict(ds[lon_name].attrs) if lon_name in ds else {}

    # Ensure minimal CF-ish coordinate attrs
    out["latitude"].attrs.setdefault("units", "degrees_north")
    out["longitude"].attrs.setdefault("units", "degrees_east")
    out["time"].attrs.setdefault("units", ds[time_name].attrs.get("units", "hours since 1900-01-01 00:00:0.0"))

    # 11) Ensure a FillValue is set for the output
    #     We'll use the same FillValue convention as the input if available, otherwise 1e37.
    fill_out = _get_fill_value(ds[rain_var]) or 1.0e37
    out["rain_rate"].attrs["_FillValue"] = float(fill_out)

    # 12) Replace NaNs with FillValue on write (classic NetCDF behavior)
    #     xarray encoding controls how missing values are stored.
    encoding = {
        "rain_rate": {
            "_FillValue": float(fill_out),
            "dtype": "float32",
            "zlib": True,
            "complevel": 4,
        },
        "time": {"dtype": "int32"},
        "latitude": {"dtype": "float32"},
        "longitude": {"dtype": "float32"},
    }

    # 13) Write output
    LOG.info("Writing output: %s", out_path)
    out.to_netcdf(out_path, encoding=encoding)
    LOG.info("Done.")


def main() -> int:
    # CLI definition
    ap = argparse.ArgumentParser(
        description="Convert WRF-derived hourly rain accumulation NetCDF into CF rain_time_dependent NetCDF."
    )
    ap.add_argument("--in", dest="in_path", required=True, help="Input NetCDF (WRF-derived).")
    ap.add_argument("--out", dest="out_path", required=True, help="Output NetCDF (rain_time_dependent).")
    ap.add_argument(
        "--rain-var",
        default=None,
        help="Name of rain accumulation variable (e.g., RAIN_DELTA or DELTA_RAIN). If omitted, auto-detect.",
    )
    ap.add_argument(
        "--accum-hours",
        type=float,
        default=1.0,
        help="Accumulation period in hours for the input variable (default: 1.0).",
    )
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = ap.parse_args()

    # Logging config
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")

    # Run conversion
    convert(
        in_path=args.in_path,
        out_path=args.out_path,
        rain_var=args.rain_var,
        accum_hours=float(args.accum_hours),
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
