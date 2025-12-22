#!/usr/bin/env python3
"""Convert weather radar VMI GeoTIFFs into CF-compliant rainfall NetCDF."""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.util
from pathlib import Path
from typing import Optional

import numpy as np
import xarray as xr


TIME_UNITS = "hours since 1900-01-01 00:00:0.0"


def dbz_to_rainrate(
    dbz: np.ndarray,
    a: float = 200.0,
    b: float = 1.6,
    dbz_min: Optional[float] = None,
    dbz_max: Optional[float] = None,
    nodata_value: Optional[float] = None,
) -> np.ndarray:
    """
    Convert reflectivity in dBZ to rain rate in mm/h using Z=a*R^b.

    Parameters
    ----------
    dbz : ndarray
        Reflectivity in dBZ.
    a, b : float
        Z–R parameters.
    dbz_min, dbz_max : float or None
        Optional clipping in dBZ (helps avoid absurd outliers).
    nodata_value : float or None
        If provided, treat that dbz value as nodata.

    Returns
    -------
    R : ndarray
        Rain rate in mm/h (float32), with NaN where nodata/invalid.
    """
    dbz = dbz.astype(np.float32, copy=False)

    # Build a validity mask
    valid = np.isfinite(dbz)
    if nodata_value is not None:
        valid &= (dbz != nodata_value)

    # Optional clipping
    if dbz_min is not None:
        valid &= (dbz >= dbz_min)
    if dbz_max is not None:
        valid &= (dbz <= dbz_max)

    # Prepare output with NaNs
    R = np.full(dbz.shape, np.nan, dtype=np.float32)
    if not np.any(valid):
        return R

    # Z (linear) from dBZ
    # dBZ = 10*log10(Z) -> Z = 10^(dBZ/10)
    Z = np.empty_like(dbz, dtype=np.float32)
    Z[valid] = np.power(10.0, dbz[valid] / 10.0, dtype=np.float32)

    # R from Z–R
    # R = (Z/a)^(1/b)
    R[valid] = np.power(Z[valid] / float(a), 1.0 / float(b), dtype=np.float32)

    # Guard against negative/inf (shouldn't happen, but be safe)
    R[~np.isfinite(R)] = np.nan
    R[R < 0] = np.nan

    return R


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


def _load_domain_reference(domain_path: str) -> xr.Dataset:
    resolved = Path(domain_path)
    if not resolved.exists():
        raise FileNotFoundError(f"Domain NetCDF not found: {resolved}")
    return xr.open_dataset(resolved)


def _template_from_domain(domain: xr.Dataset) -> xr.DataArray:
    if "latitude" not in domain.coords or "longitude" not in domain.coords:
        raise ValueError("Domain dataset must define latitude/longitude coordinates.")
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


def _regrid_to_domain(da: xr.DataArray, domain: xr.Dataset) -> xr.DataArray:
    template = _template_from_domain(domain)
    domain_crs = _extract_domain_crs(domain)
    if domain_crs:
        template = template.rio.write_crs(domain_crs, inplace=False)
        if da.rio.crs is None:
            da = da.rio.write_crs(domain_crs, inplace=False)
    try:
        return da.rio.reproject_match(template)
    except Exception:
        return da.interp(
            latitude=domain["latitude"],
            longitude=domain["longitude"],
            method="nearest",
        )


def convert_vmi_to_rain(
    input_path: str,
    output_path: str,
    timestamp: Optional[str],
    source_name: str,
    institution: str,
    source: str,
    grid_mapping_name: Optional[str],
    epsg_code: Optional[str],
    domain_path: str,
    fill_value: float,
    z_r_a: float,
    z_r_b: float,
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
    domain = _load_domain_reference(domain_path)
    da = _regrid_to_domain(da, domain)
    da = da.astype(np.float32)
    fill_value = _pick_fill_value(da, fill_value)
    nodata_value = fill_value if np.isfinite(fill_value) else None

    dt = _parse_time(timestamp)
    time_value = _time_to_hours(dt)

    rain_rate_data = dbz_to_rainrate(
        da.values,
        a=z_r_a,
        b=z_r_b,
        nodata_value=nodata_value,
    )
    rain_rate = xr.DataArray(
        rain_rate_data,
        dims=da.dims,
        coords=da.coords,
        name="rain_rate",
    ).expand_dims(time=[time_value])
    rain_rate.attrs.update(
        {
            "long_name": "rainfall_rate",
            "standard_name": "rainfall_rate",
            "units": "mm h-1",
            "grid_mapping": "crs",
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
        "--domain",
        required=True,
        help=(
            "Domain NetCDF (compliant with cdl/domain.cdl) to regrid onto; "
            "expects latitude/longitude coords."
        ),
    )
    parser.add_argument(
        "--fill-value",
        type=float,
        default=-9999.0,
        help="Fallback fill value when none is provided by the raster",
    )
    parser.add_argument(
        "--z-r-a",
        dest="z_r_a",
        type=float,
        default=200.0,
        help="Z-R parameter a in Z=a*R^b (default: 200.0)",
    )
    parser.add_argument(
        "--z-r-b",
        dest="z_r_b",
        type=float,
        default=1.6,
        help="Z-R parameter b in Z=a*R^b (default: 1.6)",
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
        domain_path=args.domain,
        fill_value=args.fill_value,
        z_r_a=args.z_r_a,
        z_r_b=args.z_r_b,
    )


if __name__ == "__main__":
    main()
