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
- hourly accumulated rain in: DELTA_RAIN (mm)
- reference domain (cdl/domain.cdl) for regridding

Output file contains:
- rain_rate(time, latitude, longitude) in mm h-1

Key feature (extended):
- You can pass multiple files (or a glob) and the script will **merge time steps**
  into a single output NetCDF, sorted by time and de-duplicated.

Conversion:
- If accumulation is in mm over Δt hours:
    rain_rate_mm_h = (accum_mm / Δt)
"""

from __future__ import annotations

import argparse
import glob
import logging
from pathlib import Path
from typing import Optional, Any, List

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


def _default_output_name() -> str:
    """Return the default filename used when --out points to a directory."""
    return "rain_rate_merged.nc"


def _resolve_output_path(out_path: str) -> Path:
    """Resolve the output path, handling directory targets gracefully."""
    target = Path(out_path)
    # Treat existing directories as output folders.
    if target.exists() and target.is_dir():
        return target / _default_output_name()
    # Treat paths with no suffix or a trailing separator as directories.
    if str(out_path).endswith(("/", "\\")) or target.suffix == "":
        target.mkdir(parents=True, exist_ok=True)
        return target / _default_output_name()
    # Ensure the parent directory exists for file targets.
    target.parent.mkdir(parents=True, exist_ok=True)
    return target


def _require_rioxarray() -> None:
    """Ensure rioxarray is available for reprojection/regridding."""
    try:
        import rioxarray as rxr  # noqa: F401
    except Exception as exc:  # pragma: no cover - handled via runtime error
        raise ImportError("rioxarray is required for regridding to the domain grid.") from exc


def _load_domain_reference(domain_path: str) -> xr.Dataset:
    """Load the reference domain NetCDF (compliant with cdl/domain.cdl)."""
    resolved = Path(domain_path)
    if not resolved.exists():
        raise FileNotFoundError(f"Domain NetCDF not found: {resolved}")
    ds = xr.open_dataset(resolved)
    if "latitude" not in ds.coords or "longitude" not in ds.coords:
        raise ValueError("Domain dataset must define latitude/longitude coordinates.")
    return ds


def _template_from_domain(domain: xr.Dataset) -> xr.DataArray:
    if domain.data_vars:
        template = next(iter(domain.data_vars.values()))
        return template.drop_vars([name for name in template.coords if name not in ("latitude", "longitude")])
    return xr.DataArray(
        np.empty((domain.dims["latitude"], domain.dims["longitude"])),
        coords={"latitude": domain["latitude"], "longitude": domain["longitude"]},
        dims=("latitude", "longitude"),
    )


def _extract_domain_crs(domain: xr.Dataset) -> Optional[str]:
    if "crs" not in domain:
        return None
    epsg_code = domain["crs"].attrs.get("epsg_code")
    if epsg_code:
        return str(epsg_code)
    grid_mapping_name = domain["crs"].attrs.get("grid_mapping_name")
    if grid_mapping_name:
        return str(grid_mapping_name)
    return None


def _prepare_for_regrid(da: xr.DataArray) -> xr.DataArray:
    """Set spatial dimension names for rioxarray operations."""
    try:
        return da.rio.set_spatial_dims(x_dim="longitude", y_dim="latitude", inplace=False)
    except Exception:
        return da


def _regrid_to_domain(da: xr.DataArray, domain: xr.Dataset) -> xr.DataArray:
    """Regrid a 2D DataArray onto the reference domain grid."""
    template = _prepare_for_regrid(_template_from_domain(domain))
    da = _prepare_for_regrid(da)
    domain_crs = _extract_domain_crs(domain)
    if domain_crs:
        template = template.rio.write_crs(domain_crs, inplace=False)
        if da.rio.crs is None:
            da = da.rio.write_crs(domain_crs, inplace=False)
    try:
        regridded = da.rio.reproject_match(template)
    except Exception:
        regridded = da.interp(
            latitude=domain["latitude"],
            longitude=domain["longitude"],
            method="nearest",
        )
    return regridded.assign_coords(latitude=domain["latitude"], longitude=domain["longitude"])


def _regrid_rain_rate_time_series(rain_rate: xr.DataArray, domain: xr.Dataset) -> xr.DataArray:
    """Apply regridding per time step to keep temporal dimension intact."""
    if "time" not in rain_rate.dims:
        return _regrid_to_domain(rain_rate, domain)

    regridded = []
    for idx in range(rain_rate.sizes["time"]):
        slice_rr = rain_rate.isel(time=idx)
        regridded_slice = _regrid_to_domain(slice_rr, domain)
        regridded_slice = regridded_slice.expand_dims(time=[rain_rate["time"].values[idx]])
        regridded.append(regridded_slice)
    return xr.concat(regridded, dim="time")


def _convert_one_dataset(
    ds: xr.Dataset,
    rain_var: str,
    accum_hours: float,
    time_name: str,
    lat_name: str,
    lon_name: str,
) -> xr.Dataset:
    """Convert one dataset to an intermediate dataset with rain_rate."""
    # Convert input accumulation into a floating rain field with missing values masked.
    rain_acc_mm = _mask_fill(ds[rain_var])

    for req in (time_name, lat_name, lon_name):
        if req not in rain_acc_mm.dims:
            raise ValueError(f"Rain variable dims {rain_acc_mm.dims} do not include required coord {req}.")

    # Convert mm over Δt hours into mm h-1 for the output forcing format.
    rain_rate = rain_acc_mm / float(accum_hours)

    # Attach the required CF-style metadata for rainfall forcing files.
    rain_rate.attrs = {
        "standard_name": "rainfall_rate",
        "long_name": "rainfall_rate",
        "units": "mm h-1",
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
    domain_path: str,
) -> None:
    """Convert multiple files and merge into a single time series output."""
    _require_rioxarray()

    paths = _expand_inputs(in_paths)
    if not paths:
        raise ValueError("No input files found. Check --in arguments/globs.")
    resolved_out_path = _resolve_output_path(out_path)

    domain = _load_domain_reference(domain_path)
    LOG.info("Loaded reference domain from %s", domain_path)

    LOG.info("Inputs (%d):", len(paths))
    for p in paths:
        LOG.info("  - %s", p)

    intermediate: List[xr.Dataset] = []

    chosen_rain_var: str = rain_var or "DELTA_RAIN"
    fill_out: float = -9999.0
    time_units: Optional[str] = None
    # Output fill value is fixed to match the expected forcing format.

    for idx, path in enumerate(paths):
        LOG.info("Opening input %d/%d: %s", idx + 1, len(paths), path)
        ds = xr.open_dataset(path)

        lat_name = _infer_coord_name(ds, ["latitude", "lat", "y"])
        lon_name = _infer_coord_name(ds, ["longitude", "lon", "x"])
        time_name = _infer_coord_name(ds, ["time"])
        if not lat_name or not lon_name or not time_name:
            raise ValueError(f"Could not infer coords in {path}. Found lat={lat_name} lon={lon_name} time={time_name}")

        if chosen_rain_var != "DELTA_RAIN":
            raise ValueError("This converter only supports the DELTA_RAIN variable.")

        if chosen_rain_var not in ds.variables:
            raise KeyError(f"Rain variable '{chosen_rain_var}' not found in {path}. Vars: {list(ds.variables)}")

        # Preserve input fill values for masking, but keep a fixed output fill value.
        if time_units is None:
            time_units = str(ds[time_name].attrs.get("units", "hours since 1900-01-01 00:00:0.0"))

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

        regridded_rr = _regrid_rain_rate_time_series(inter["rain_rate"], domain)
        regridded = xr.Dataset(
            data_vars={"rain_rate": regridded_rr},
            coords={"time": inter["time"], "latitude": domain["latitude"], "longitude": domain["longitude"]},
        )
        intermediate.append(regridded)

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

    # Match the expected global attributes for rain forcing inputs.
    merged.attrs = {
        "Conventions": "CF-1.10",
        "title": "Stations rainfall forcing",
        "institution": "Unknown",
        "source": "stations",
        "history": f"{np.datetime_as_string(np.datetime64('now'), unit='s')}: produced/ingested",
    }

    # Ensure coordinate metadata follows the CDL template expected by downstream tools.
    # NOTE: datetime64 time coordinates cannot keep a "units" attribute without
    # clashing with xarray's encoding step, so we move units into encoding below.
    time_units_to_use = time_units or "hours since 1900-01-01 00:00:0.0"
    merged["time"].attrs = {
        "description": "Time",
        "long_name": "time",
    }
    if not np.issubdtype(merged["time"].dtype, np.datetime64):
        merged["time"].attrs["units"] = time_units_to_use
    merged["latitude"].attrs = {
        "description": "Latitude",
        "long_name": "latitude",
        "units": "degrees_north",
    }
    merged["longitude"].attrs = {
        "description": "Longitude",
        "long_name": "longitude",
        "units": "degrees_east",
    }

    # Normalize rain rate attributes and add grid mapping reference.
    merged["rain_rate"].attrs = {
        "long_name": "rainfall_rate",
        "standard_name": "rainfall_rate",
        "units": "mm h-1",
        "grid_mapping": "crs",
    }

    # Add or normalize the grid mapping variable.
    if "crs" in domain:
        merged["crs"] = domain["crs"].astype("int64")
    else:
        merged["crs"] = xr.DataArray(np.int64(0))

    # Define encodings so xarray writes the expected dtypes and fill values.
    encoding = {
        "rain_rate": {"_FillValue": float(fill_out), "dtype": "float32", "zlib": True, "complevel": 4},
        "time": {"_FillValue": np.nan, "dtype": "float64"},
        "latitude": {"_FillValue": np.nan, "dtype": "float32"},
        "longitude": {"_FillValue": np.nan, "dtype": "float32"},
        "crs": {"dtype": "int64"},
    }
    if np.issubdtype(merged["time"].dtype, np.datetime64):
        encoding["time"]["units"] = time_units_to_use

    if Path(out_path).resolve() != resolved_out_path.resolve():
        LOG.info("Resolved output path from %s to %s", out_path, resolved_out_path)
    LOG.info("Writing output: %s", resolved_out_path)
    merged.to_netcdf(resolved_out_path, encoding=encoding)
    LOG.info("Done. Time steps written: %d", int(merged.dims.get("time", 0)))


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Convert WRF-derived rain accumulation NetCDF(s) into a merged CF rain_time_dependent NetCDF."
    )
    ap.add_argument("--in", dest="in_paths", required=True, nargs="+",
                    help="One or more input NetCDFs (paths or globs). Example: wrf5_d02_20251221Z*.nc")
    ap.add_argument("--out", dest="out_path", required=True, help="Output merged NetCDF.")
    ap.add_argument("--rain-var", default=None,
                    help="Rain accumulation variable name (only DELTA_RAIN is supported).")
    ap.add_argument("--accum-hours", type=float, default=1.0,
                    help="Accumulation period in hours for the input variable (default: 1.0).")
    ap.add_argument("--domain", required=True,
                    help="Domain NetCDF (matching cdl/domain.cdl) used as the target grid for regridding.")
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
        domain_path=args.domain,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
