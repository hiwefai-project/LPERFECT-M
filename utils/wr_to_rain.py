#!/usr/bin/env python3
"""Convert weather radar VMI GeoTIFFs into CF-compliant rainfall NetCDF."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
from typing import Optional

import numpy as np
import xarray as xr


TIME_UNITS = "hours since 1900-01-01 00:00:0.0"


def _require_rioxarray() -> None:
    if importlib.util.find_spec("rioxarray") is None:
        raise ImportError(
            "Reading GeoTIFF inputs requires rioxarray and rasterio. "
            "Install them or convert the input to NetCDF."
        )


def _parse_time(value: Optional[str]) -> datetime:
    if value is None:
        return datetime.now(timezone.utc)
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    dt = datetime.fromisoformat(text)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


def _time_to_hours(dt: datetime) -> int:
    base = datetime(1900, 1, 1, tzinfo=timezone.utc)
    hours = (dt - base).total_seconds() / 3600.0
    return int(round(hours))


def _normalize_dims(da: xr.DataArray) -> xr.DataArray:
    if "latitude" in da.dims and "longitude" in da.dims:
        return da.transpose("latitude", "longitude")
    rename_map = {}
    if "y" in da.dims:
        rename_map["y"] = "latitude"
    if "x" in da.dims:
        rename_map["x"] = "longitude"
    if rename_map:
        da = da.rename(rename_map)
    if da.dims != ("latitude", "longitude"):
        raise ValueError(f"Expected 2D raster with x/y or latitude/longitude dims, got {da.dims}")
    return da.transpose("latitude", "longitude")


def _crs_attrs(da: xr.DataArray, grid_mapping_name: Optional[str], epsg_code: Optional[str]) -> dict:
    attrs: dict[str, str] = {}
    crs = da.rio.crs if hasattr(da, "rio") else None
    if grid_mapping_name is None:
        if crs is not None:
            try:
                cf = crs.to_cf()
                grid_mapping_name = cf.get("grid_mapping_name")
            except Exception:
                grid_mapping_name = None
    if epsg_code is None and crs is not None:
        epsg = crs.to_epsg()
        if epsg is not None:
            epsg_code = f"EPSG:{epsg}"
    attrs["grid_mapping_name"] = grid_mapping_name or "unknown"
    attrs["epsg_code"] = epsg_code or "EPSG:0000"
    return attrs


def _pick_fill_value(da: xr.DataArray, default: float) -> float:
    for key in ("_FillValue", "missing_value"):
        if key in da.attrs:
            return float(da.attrs[key])
    if np.issubdtype(da.dtype, np.floating) and np.isnan(da).any():
        return np.nan
    return default


def convert_vmi_to_rain(
    input_path: str,
    output_path: str,
    timestamp: Optional[str],
    source_name: str,
    institution: str,
    source: str,
    grid_mapping_name: Optional[str],
    epsg_code: Optional[str],
    fill_value: float,
) -> None:
    _require_rioxarray()
    import rioxarray as rxr

    da = rxr.open_rasterio(input_path)
    if "band" in da.dims:
        if da.sizes["band"] < 1:
            raise ValueError(f"No bands found in raster {input_path}")
        da = da.isel(band=0, drop=True)
    da = da.rename("rain_rate")
    da = _normalize_dims(da)
    da = da.astype(np.float32)
    fill_value = _pick_fill_value(da, fill_value)

    dt = _parse_time(timestamp)
    time_value = _time_to_hours(dt)

    rain_rate = da.expand_dims(time=[time_value])
    rain_rate.attrs.update(
        {
            "long_name": "rainfall_rate",
            "standard_name": "rainfall_rate",
            "units": "mm h-1",
            "grid_mapping": "crs",
            "_FillValue": fill_value,
        }
    )

    ds = xr.Dataset(
        {
            "rain_rate": rain_rate,
            "crs": xr.DataArray(0, attrs=_crs_attrs(da, grid_mapping_name, epsg_code)),
        },
        coords={
            "time": ("time", np.array([time_value], dtype=np.int32)),
            "latitude": ("latitude", da["latitude"].values),
            "longitude": ("longitude", da["longitude"].values),
        },
        attrs={
            "Conventions": "CF-1.10",
            "title": f"{source_name} rainfall forcing",
            "institution": institution,
            "source": source,
            "history": f"{dt.isoformat()}: produced/ingested",
        },
    )

    ds["time"].attrs.update(
        {
            "description": "Time",
            "long_name": "time",
            "units": TIME_UNITS,
        }
    )
    ds["latitude"].attrs.update(
        {
            "description": "Latitude",
            "long_name": "latitude",
            "units": "degrees_north",
        }
    )
    ds["longitude"].attrs.update(
        {
            "description": "Longitude",
            "long_name": "longitude",
            "units": "degrees_east",
        }
    )

    encoding = {"rain_rate": {"_FillValue": fill_value}}
    ds.to_netcdf(output_path, encoding=encoding)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Convert a weather radar VMI GeoTIFF to a CF-compliant rainfall NetCDF."
    )
    parser.add_argument("--input", required=True, help="Path to input VMI GeoTIFF")
    parser.add_argument("--output", required=True, help="Path to output NetCDF")
    parser.add_argument(
        "--time",
        dest="timestamp",
        default=None,
        help="Timestamp for the raster (ISO-8601, defaults to now UTC)",
    )
    parser.add_argument("--source-name", default="Radar", help="Source name for global title")
    parser.add_argument("--institution", default="Unknown", help="Institution metadata")
    parser.add_argument("--source", default="radar", help="Source metadata")
    parser.add_argument(
        "--grid-mapping-name",
        default=None,
        help="Override grid_mapping_name (defaults from CRS if available)",
    )
    parser.add_argument(
        "--epsg",
        dest="epsg_code",
        default=None,
        help="Override EPSG code (e.g., EPSG:4326)",
    )
    parser.add_argument(
        "--fill-value",
        type=float,
        default=-9999.0,
        help="Fallback fill value when none is provided by the raster",
    )
    args = parser.parse_args()

    convert_vmi_to_rain(
        input_path=args.input,
        output_path=args.output,
        timestamp=args.timestamp,
        source_name=args.source_name,
        institution=args.institution,
        source=args.source,
        grid_mapping_name=args.grid_mapping_name,
        epsg_code=args.epsg_code,
        fill_value=args.fill_value,
    )


if __name__ == "__main__":
    main()
