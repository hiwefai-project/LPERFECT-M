#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_wrf_rain_to_rain_time_dependent.py

Convert one or many WRF-derived NetCDF files (wrf5_d02_*.nc) into a single
time-series NetCDF compliant with a simple time-dependent rain-rate CDL
(cdl/rain_time_dependent.cdl).

Input files (per the provided CDL) typically have:
- dimensions: time, latitude, longitude
- coordinates: time(time), latitude(latitude), longitude(longitude)
- hourly accumulated rain in: RAIN_DELTA or DELTA_RAIN (mm)

Output file contains:
- rain_rate(time, latitude, longitude) in kg m-2 s-1

Key feature (extended):
- You can pass multiple files (or a glob) and the script will **merge time steps**
  into a single output NetCDF, sorted by time and de-duplicated.

Conversion:
- If accumulation is in mm over Δt hours:
    rain_rate_kg_m2_s = (accum_mm / (Δt * 3600))
  using 1 mm == 1 kg m-2 (liquid water equivalent).
"""

from __future__ import annotations

import argparse
import glob
import logging
from typing import Optional, Any, List, Tuple

import numpy as np
import xarray as xr

LOG = logging.getLogger("wrf_rain_converter")


def _to_float(x: Any) -> Optional[float]:
    """Safely convert a scalar to float; return None on failure."""
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _infer_coord_name(ds: xr.Dataset, candidates: List[str]) -> Optional[str]:
    """Return the first candidate that exists as a coord or variable in ds."""
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
    """Replace FillValue with NaN and ensure float dtype."""
    out = da.astype(float)
    fill = _get_fill_value(out)
    if fill is not None:
        out = out.where(out != fill)
    return out


def _expand_inputs(inputs: List[str]) -> List[str]:
    """
    Expand a list of file paths / globs into an explicit list of file paths.
    Keeps order stable but removes duplicates.
    """
    expanded: List[str] = []
    for item in inputs:
        matches = glob.glob(item)
        if matches:
            expanded.extend(matches)
        else:
            expanded.append(item)

    seen = set()
    out: List[str] = []
    for p in expanded:
        if p not in seen:
            out.append(p)
            seen.add(p)
    return out


def _assert_same_grid(
    base_lat: np.ndarray,
    base_lon: np.ndarray,
    lat: np.ndarray,
    lon: np.ndarray,
    path: str,
) -> None:
    """Ensure that multiple inputs are on the same grid (required for safe merge)."""
    if base_lat.shape != lat.shape or base_lon.shape != lon.shape:
        raise ValueError(
            f"Grid shape mismatch in {path}: lat {lat.shape} lon {lon.shape} "
            f"vs base lat {base_lat.shape} lon {base_lon.shape}"
        )
    if not np.allclose(base_lat, lat, rtol=0, atol=1e-10) or not np.allclose(base_lon, lon, rtol=0, atol=1e-10):
        raise ValueError(f"Grid coordinate mismatch in {path}: lat/lon values differ from base grid.")


def _convert_one_dataset(
    ds: xr.Dataset,
    rain_var: str,
    accum_hours: float,
    time_name: str,
    lat_name: str,
    lon_name: str,
) -> xr.Dataset:
    """Convert one dataset to an intermediate dataset with rain_rate."""
    rain_acc_mm = _mask_fill(ds[rain_var])

    for req in (time_name, lat_name, lon_name):
        if req not in rain_acc_mm.dims:
            raise ValueError(f"Rain variable dims {rain_acc_mm.dims} do not include required coord {req}.")

    seconds = float(accum_hours) * 3600.0
    rain_rate = rain_acc_mm / seconds

    rain_rate.attrs = {
        "standard_name": "lwe_precipitation_rate",
        "long_name": "rain_rate",
        "description": f"Rain rate derived from {rain_var} (accumulated over {accum_hours} h)",
        "units": "kg m-2 s-1",
        "source_variable": rain_var,
        "accumulation_period_hours": accum_hours,
    }

    return xr.Dataset(
        data_vars={"rain_rate": rain_rate},
        coords={time_name: ds[time_name], lat_name: ds[lat_name], lon_name: ds[lon_name]},
    )


def convert_many(
    in_paths: List[str],
    out_path: str,
    rain_var: Optional[str],
    accum_hours: float,
    sort_time: bool,
    dedupe_time: bool,
) -> None:
    """Convert multiple files and merge into a single time series output."""
    paths = _expand_inputs(in_paths)
    if not paths:
        raise ValueError("No input files found. Check --in arguments/globs.")

    LOG.info("Inputs (%d):", len(paths))
    for p in paths:
        LOG.info("  - %s", p)

    intermediate: List[xr.Dataset] = []

    base_lat = None
    base_lon = None

    chosen_rain_var: Optional[str] = rain_var
    fill_out: Optional[float] = None
    time_units: Optional[str] = None
    lat_attrs: dict = {}
    lon_attrs: dict = {}
    time_attrs: dict = {}

    for idx, path in enumerate(paths):
        LOG.info("Opening input %d/%d: %s", idx + 1, len(paths), path)
        ds = xr.open_dataset(path)

        lat_name = _infer_coord_name(ds, ["latitude", "lat", "y"])
        lon_name = _infer_coord_name(ds, ["longitude", "lon", "x"])
        time_name = _infer_coord_name(ds, ["time"])
        if not lat_name or not lon_name or not time_name:
            raise ValueError(f"Could not infer coords in {path}. Found lat={lat_name} lon={lon_name} time={time_name}")

        if chosen_rain_var is None:
            for cand in ["RAIN_DELTA", "DELTA_RAIN", "HOURLY_RAIN", "RAIN"]:
                if cand in ds.variables:
                    chosen_rain_var = cand
                    break
            if chosen_rain_var is None:
                raise KeyError(f"Could not auto-detect rain var in {path}. Vars: {list(ds.variables)}")

        if chosen_rain_var not in ds.variables:
            raise KeyError(f"Rain variable '{chosen_rain_var}' not found in {path}. Vars: {list(ds.variables)}")

        if fill_out is None:
            fill_out = _get_fill_value(ds[chosen_rain_var]) or 1.0e37
        if time_units is None:
            time_units = str(ds[time_name].attrs.get("units", "hours since 1900-01-01 00:00:0.0"))

        if idx == 0:
            base_lat = np.asarray(ds[lat_name].values, dtype=float)
            base_lon = np.asarray(ds[lon_name].values, dtype=float)
            lat_attrs = dict(ds[lat_name].attrs)
            lon_attrs = dict(ds[lon_name].attrs)
            time_attrs = dict(ds[time_name].attrs)
        else:
            _assert_same_grid(
                base_lat=base_lat,
                base_lon=base_lon,
                lat=np.asarray(ds[lat_name].values, dtype=float),
                lon=np.asarray(ds[lon_name].values, dtype=float),
                path=path,
            )

        inter = _convert_one_dataset(
            ds=ds,
            rain_var=chosen_rain_var,
            accum_hours=accum_hours,
            time_name=time_name,
            lat_name=lat_name,
            lon_name=lon_name,
        )

        rename_map = {}
        if time_name != "time":
            rename_map[time_name] = "time"
        if lat_name != "latitude":
            rename_map[lat_name] = "latitude"
        if lon_name != "longitude":
            rename_map[lon_name] = "longitude"
        if rename_map:
            inter = inter.rename(rename_map)

        intermediate.append(inter)

    assert chosen_rain_var is not None

    LOG.info("Concatenating %d datasets along time...", len(intermediate))
    merged = xr.concat(intermediate, dim="time")

    if sort_time:
        LOG.info("Sorting by time coordinate...")
        merged = merged.sortby("time")

    if dedupe_time:
        LOG.info("De-duplicating identical time values (keep first)...")
        tvals = np.asarray(merged["time"].values)
        _, first_idx = np.unique(tvals, return_index=True)
        merged = merged.isel(time=np.sort(first_idx))

    merged.attrs = {
        "Conventions": "CF-1.10",
        "title": "Time-dependent rain rate (converted from WRF-derived accumulation)",
        "source": "WRF post-processed output",
        "history": "Converted and merged by convert_wrf_rain_to_rain_time_dependent.py",
        "input_files": ", ".join(paths),
        "input_rain_variable": chosen_rain_var,
        "accumulation_period_hours": float(accum_hours),
    }

    merged["time"].attrs = time_attrs
    merged["latitude"].attrs = lat_attrs
    merged["longitude"].attrs = lon_attrs
    merged["latitude"].attrs.setdefault("units", "degrees_north")
    merged["longitude"].attrs.setdefault("units", "degrees_east")
    merged["time"].attrs.setdefault("units", time_units or "hours since 1900-01-01 00:00:0.0")

    fill_out = float(fill_out or 1.0e37)
    merged["rain_rate"].attrs["_FillValue"] = fill_out

    encoding = {
        "rain_rate": {"_FillValue": fill_out, "dtype": "float32", "zlib": True, "complevel": 4},
        "time": {"dtype": "int32"},
        "latitude": {"dtype": "float32"},
        "longitude": {"dtype": "float32"},
    }

    LOG.info("Writing output: %s", out_path)
    merged.to_netcdf(out_path, encoding=encoding)
    LOG.info("Done. Time steps written: %d", int(merged.dims.get("time", 0)))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert WRF-derived rain accumulation NetCDF(s) into a merged CF rain_time_dependent NetCDF."
    )
    ap.add_argument("--in", dest="in_paths", required=True, nargs="+",
                    help="One or more input NetCDFs (paths or globs). Example: wrf5_d02_20251221Z*.nc")
    ap.add_argument("--out", dest="out_path", required=True, help="Output merged NetCDF.")
    ap.add_argument("--rain-var", default=None,
                    help="Rain accumulation variable name (RAIN_DELTA or DELTA_RAIN). If omitted, auto-detect from first file.")
    ap.add_argument("--accum-hours", type=float, default=1.0,
                    help="Accumulation period in hours for the input variable (default: 1.0).")
    ap.add_argument("--no-sort", action="store_true", help="Do not sort time after concatenation.")
    ap.add_argument("--no-dedupe", action="store_true", help="Do not remove duplicate time values.")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = ap.parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")

    convert_many(
        in_paths=list(args.in_paths),
        out_path=args.out_path,
        rain_var=args.rain_var,
        accum_hours=float(args.accum_hours),
        sort_time=not args.no_sort,
        dedupe_time=not args.no_dedupe,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
