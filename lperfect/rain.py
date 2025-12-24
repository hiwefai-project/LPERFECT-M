# -*- coding: utf-8 -*-
"""Rain input and blending (rank0 read, broadcast)."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import dataclass for structured source definition.
from dataclasses import dataclass  # import dataclasses import dataclass

# Import typing primitives.
from typing import Any, Dict, List, Optional, Tuple  # import typing import Any, Dict, List, Optional, Tuple

# Import logging.
import logging  # import logging

# Import numpy.
import numpy as np  # import numpy as np

# Import xarray.
import xarray as xr  # import xarray as xr

# Import CF schema helpers.
from .cf_schema import (  # import .cf_schema
    RAIN_CRS_VAR,
    RAIN_GRID_MAPPING_ATTR,
    RAIN_LAT_VAR,
    RAIN_LON_VAR,
    RAIN_RATE_UNITS,
    RAIN_RATE_VAR,
    RAIN_TIME_UNITS,
    RAIN_TIME_VAR,
    hours_since_1900_to_datetime64,
    normalize_cf_time_units,
)

logger = logging.getLogger("lperfect.rain")  # set logger


@dataclass  # apply decorator
class RainSource:  # define class RainSource
    """Rainfall source configuration."""  # execute statement
    name: str  # execute statement
    kind: str  # execute statement
    weight: float  # execute statement
    mode: str  # execute statement
    path: Optional[str] = None  # execute statement
    var: Optional[str] = None  # execute statement
    time_var: str = RAIN_TIME_VAR  # execute statement
    lat_var: str = RAIN_LAT_VAR  # execute statement
    lon_var: str = RAIN_LON_VAR  # execute statement
    crs_var: str = RAIN_CRS_VAR  # execute statement
    time_units: str = RAIN_TIME_UNITS  # execute statement
    rate_units: str = RAIN_RATE_UNITS  # execute statement
    require_cf: bool = True  # execute statement
    require_time_dim: bool = True  # execute statement
    select: str = "nearest"  # execute statement
    value: Optional[float] = None  # execute statement


# Simple cache to avoid reopening NetCDF files each step.
_NC_CACHE: Dict[str, xr.Dataset] = {}  # execute statement
_NC_SCHEMA_CACHE: Dict[str, bool] = {}  # execute statement
# Track the last rain time index logged per source to avoid duplicate messages.
_RAIN_TIME_LOG: Dict[Tuple[str, str], int] = {}  # execute statement


def xr_open_cached(path: str) -> xr.Dataset:  # define function xr_open_cached
    """Open a NetCDF dataset with caching."""  # execute statement
    if path not in _NC_CACHE:  # check condition path not in _NC_CACHE:
        _NC_CACHE[path] = xr.open_dataset(path, decode_times=True)  # execute statement
    return _NC_CACHE[path]  # return _NC_CACHE[path]


def xr_close_cache() -> None:  # define function xr_close_cache
    """Close all cached datasets."""  # execute statement
    for ds in _NC_CACHE.values():  # loop over ds in _NC_CACHE.values():
        try:  # start exception handling
            ds.close()  # execute statement
        except Exception:  # handle exception Exception:
            pass  # no-op placeholder
    _NC_CACHE.clear()  # execute statement
    _NC_SCHEMA_CACHE.clear()  # execute statement


def build_rain_sources(cfg: Dict[str, Any]) -> List[RainSource]:  # define function build_rain_sources
    """Parse cfg['rain']['sources'] into a list of RainSource."""  # execute statement
    rain_cfg = cfg.get("rain", {})  # set rain_cfg
    schema_cfg = rain_cfg.get("schema", {})  # set schema_cfg
    sources_cfg = rain_cfg.get("sources", {})  # set sources_cfg
    out: List[RainSource] = []  # execute statement
    for name, sc in sources_cfg.items():  # loop over name, sc in sources_cfg.items():
        out.append(RainSource(  # execute statement
            name=name,  # set name
            kind=str(sc.get("kind", "netcdf")),  # set kind
            weight=float(sc.get("weight", 0.0)),  # set weight
            mode=str(sc.get("mode", "intensity_mmph")),  # set mode
            path=sc.get("path", None),  # set path
            var=sc.get("var", schema_cfg.get("rain_var", RAIN_RATE_VAR)),  # set var
            time_var=str(sc.get("time_var", schema_cfg.get("time_var", RAIN_TIME_VAR))),  # set time_var
            lat_var=str(sc.get("lat_var", schema_cfg.get("lat_var", RAIN_LAT_VAR))),  # set lat_var
            lon_var=str(sc.get("lon_var", schema_cfg.get("lon_var", RAIN_LON_VAR))),  # set lon_var
            crs_var=str(sc.get("crs_var", schema_cfg.get("crs_var", RAIN_CRS_VAR))),  # set crs_var
            time_units=str(sc.get("time_units", schema_cfg.get("time_units", RAIN_TIME_UNITS))),  # set time_units
            rate_units=str(sc.get("rate_units", schema_cfg.get("rate_units", RAIN_RATE_UNITS))),  # set rate_units
            require_cf=bool(sc.get("require_cf", schema_cfg.get("require_cf", True))),  # set require_cf
            require_time_dim=bool(sc.get("require_time_dim", schema_cfg.get("require_time_dim", True))),  # set require_time_dim
            select=str(sc.get("select", "nearest")),  # set select
            value=sc.get("value", None),  # set value
        ))  # execute statement
    return out  # return out


def rain_to_step_mm(field: np.ndarray, mode: str, dt_s: float) -> np.ndarray:  # define function rain_to_step_mm
    """Convert rainfall field to mm per model step."""  # execute statement
    f = np.asarray(field, dtype=np.float64)  # set f
    f = np.where(np.isfinite(f), f, 0.0)  # set f
    f = np.maximum(f, 0.0)  # set f
    if mode == "intensity_mmph":  # check condition mode == "intensity_mmph":
        return f * (dt_s / 3600.0)  # return f * (dt_s / 3600.0)
    if mode == "depth_mm_per_step":  # check condition mode == "depth_mm_per_step":
        return f  # return f
    raise ValueError(f"Unknown rain mode '{mode}'")  # raise ValueError(f"Unknown rain mode '{mode}'")


def pick_time_index(time_vals: np.ndarray, target: np.datetime64) -> int:  # define function pick_time_index
    """Pick nearest time index for a datetime64 time axis."""  # execute statement
    tv = np.asarray(time_vals)  # set tv
    if tv.dtype.kind != "M":  # check condition tv.dtype.kind != "M":
        raise ValueError("Time axis is not datetime64; verify CF time decoding.")  # raise ValueError("Time axis is not datetime64; verify CF time decoding.")
    return int(np.argmin(np.abs(tv - target)))  # return int(np.argmin(np.abs(tv - target)))


def _decode_cf_time_axis(ds: xr.Dataset, src: RainSource) -> np.ndarray:  # define function _decode_cf_time_axis
    """Return datetime64 time axis from a CF time coordinate."""  # execute statement
    if src.time_var not in ds:  # check condition src.time_var not in ds:
        raise ValueError(f"Rain dataset missing time variable '{src.time_var}'")  # raise ValueError(f"Rain dataset missing time variable '{src.time_var}'")
    time_da = ds[src.time_var]  # set time_da
    time_vals = np.asarray(time_da.values)  # set time_vals
    if time_vals.dtype.kind == "M":  # check condition time_vals.dtype.kind == "M":
        return time_vals  # return time_vals
    if src.require_cf and "units" not in time_da.attrs:  # check condition src.require_cf and "units" not in time_da.attrs:
        raise ValueError(f"Rain time variable '{src.time_var}' missing CF units attribute")  # raise ValueError(f"Rain time variable '{src.time_var}' missing CF units attribute")
    units = normalize_cf_time_units(str(time_da.attrs.get("units", src.time_units)))  # set units
    expected = normalize_cf_time_units(src.time_units)  # set expected
    if src.require_cf and units != expected:  # check condition src.require_cf and units != expected:
        raise ValueError(f"Rain time units '{units}' do not match expected '{expected}'")  # raise ValueError(f"Rain time units '{units}' do not match expected '{expected}'")
    return hours_since_1900_to_datetime64(time_vals)  # return hours_since_1900_to_datetime64(time_vals)


def _validate_rain_dataset(ds: xr.Dataset, src: RainSource) -> None:  # define function _validate_rain_dataset
    """Validate a rain dataset against the CF time-dependent schema."""  # execute statement
    if not src.require_cf:  # check condition not src.require_cf:
        return  # return None
    if src.var not in ds:  # check condition src.var not in ds:
        raise ValueError(f"Rain variable '{src.var}' not found in {src.path}")  # raise ValueError(f"Rain variable '{src.var}' not found in {src.path}")
    da = ds[src.var]  # set da

    if src.time_var not in ds:  # check condition src.time_var not in ds:
        raise ValueError(f"Rain dataset missing time coordinate '{src.time_var}'")  # raise ValueError(f"Rain dataset missing time coordinate '{src.time_var}'")
    if src.lat_var not in ds:  # check condition src.lat_var not in ds:
        raise ValueError(f"Rain dataset missing latitude coordinate '{src.lat_var}'")  # raise ValueError(f"Rain dataset missing latitude coordinate '{src.lat_var}'")
    if src.lon_var not in ds:  # check condition src.lon_var not in ds:
        raise ValueError(f"Rain dataset missing longitude coordinate '{src.lon_var}'")  # raise ValueError(f"Rain dataset missing longitude coordinate '{src.lon_var}'")

    if src.require_time_dim and src.time_var not in da.dims:  # check condition src.require_time_dim and src.time_var not in da.dims:
        raise ValueError(f"Rain variable '{src.var}' must include time dimension '{src.time_var}'")  # raise ValueError(f"Rain variable '{src.var}' must include time dimension '{src.time_var}'")
    if src.lat_var not in da.dims or src.lon_var not in da.dims:  # check condition src.lat_var not in da.dims or src.lon_var not in da.dims:
        raise ValueError(f"Rain variable '{src.var}' must include dims '{src.lat_var}', '{src.lon_var}'")  # raise ValueError(f"Rain variable '{src.var}' must include dims '{src.lat_var}', '{src.lon_var}'")

    lat_units = str(ds[src.lat_var].attrs.get("units", "")).lower().strip()  # set lat_units
    lon_units = str(ds[src.lon_var].attrs.get("units", "")).lower().strip()  # set lon_units
    if src.require_cf and (not lat_units or "degrees_north" not in lat_units):  # check condition src.require_cf and (not lat_units or "degrees_north" not in lat_units):
        raise ValueError(f"Latitude '{src.lat_var}' missing CF units 'degrees_north'")  # raise ValueError(f"Latitude '{src.lat_var}' missing CF units 'degrees_north'")
    if src.require_cf and (not lon_units or "degrees_east" not in lon_units):  # check condition src.require_cf and (not lon_units or "degrees_east" not in lon_units):
        raise ValueError(f"Longitude '{src.lon_var}' missing CF units 'degrees_east'")  # raise ValueError(f"Longitude '{src.lon_var}' missing CF units 'degrees_east'")

    grid_mapping = str(da.attrs.get(RAIN_GRID_MAPPING_ATTR, ""))  # set grid_mapping
    if src.crs_var and (src.crs_var not in ds or grid_mapping != src.crs_var):  # check condition src.crs_var and (src.crs_var not in ds or grid_mapping != src.crs_var):
        raise ValueError(f"Rain variable '{src.var}' must reference grid mapping '{src.crs_var}'")  # raise ValueError(f"Rain variable '{src.var}' must reference grid mapping '{src.crs_var}'")

    if src.require_cf and "units" not in da.attrs:  # check condition src.require_cf and "units" not in da.attrs:
        raise ValueError(f"Rain variable '{src.var}' missing CF units attribute")  # raise ValueError(f"Rain variable '{src.var}' missing CF units attribute")
    units = str(da.attrs.get("units", "")).strip()  # set units
    if units and units != src.rate_units:  # check condition units and units != src.rate_units:
        raise ValueError(f"Rain units '{units}' do not match expected '{src.rate_units}'")  # raise ValueError(f"Rain units '{units}' do not match expected '{src.rate_units}'")


def _format_time_value(val: Any) -> str:  # define function _format_time_value
    """Format a time coordinate value for logging."""  # execute statement
    if isinstance(val, np.datetime64):  # check condition isinstance(val, np.datetime64):
        return np.datetime_as_string(val, unit="s")  # return np.datetime_as_string(val, unit="s")
    return str(val)  # return str(val)


def _log_rain_time_usage(src: RainSource, time_vals: Optional[np.ndarray], idx: int) -> None:  # define function _log_rain_time_usage
    """Log when a rain source advances to a new time index."""  # execute statement
    key = (src.name, src.path or src.var or src.kind)  # set key
    last_idx = _RAIN_TIME_LOG.get(key, None)  # set last_idx
    if last_idx == idx:  # check condition last_idx == idx:
        return  # return None

    time_label = None  # set time_label
    if time_vals is not None:  # check condition time_vals is not None:
        arr = np.asarray(time_vals)  # set arr
        if arr.size > idx >= 0:  # check condition arr.size > idx >= 0:
            time_label = _format_time_value(arr[idx])  # set time_label

    logger.info(  # execute statement
        "Rain source '%s': using rain rate at %s (index=%d)",  # set format string
        src.name,  # execute statement
        time_label if time_label is not None else f"index {idx}",  # execute statement
        idx,  # execute statement
    )  # execute statement
    _RAIN_TIME_LOG[key] = idx  # execute statement


def blended_rain_step_mm_rank0(  # define function blended_rain_step_mm_rank0
    sources: List[RainSource],  # execute statement
    shape: Tuple[int, int],  # execute statement
    dt_s: float,  # execute statement
    step_idx: int,  # execute statement
    sim_time: Optional[np.datetime64],  # execute statement
) -> np.ndarray:  # execute statement
    """Compute blended rainfall (mm/step) on rank0."""  # execute statement
    H, W = shape  # set H, W
    total = np.zeros((H, W), dtype=np.float64)  # set total

    for src in sources:  # loop over src in sources:
        if src.weight == 0.0:  # check condition src.weight == 0.0:
            continue  # continue loop

        if src.kind == "scalar":  # check condition src.kind == "scalar":
            if src.value is None:  # check condition src.value is None:
                raise ValueError(f"Rain source '{src.name}' scalar requires 'value'")  # raise ValueError(f"Rain source '{src.name}' scalar requires 'value'")
            field = np.full((H, W), float(src.value), dtype=np.float64)  # set field

        elif src.kind == "netcdf":  # check alternate condition src.kind == "netcdf":
            if not src.path or not src.var:  # check condition not src.path or not src.var:
                raise ValueError(f"Rain source '{src.name}' netcdf requires 'path' and 'var'")  # raise ValueError(f"Rain source '{src.name}' netcdf requires 'path' and 'var'")
            ds = xr_open_cached(src.path)  # set ds
            if not _NC_SCHEMA_CACHE.get(src.path, False):  # check condition not _NC_SCHEMA_CACHE.get(src.path, False):
                _validate_rain_dataset(ds, src)  # execute statement
                _NC_SCHEMA_CACHE[src.path] = True  # execute statement
            da = ds[src.var]  # set da

            if da.ndim == 2:  # check condition da.ndim == 2:
                field = np.asarray(da.values)  # set field
            elif da.ndim == 3:  # check alternate condition da.ndim == 3:
                tdim = src.time_var if src.time_var in da.dims else da.dims[0]  # set tdim
                time_vals: Optional[np.ndarray] = None  # set time_vals
                if src.select == "step" or sim_time is None:  # check condition src.select == "step" or sim_time is None:
                    it = min(step_idx, da.sizes[tdim] - 1)  # set it
                    try:  # start exception handling
                        time_vals = _decode_cf_time_axis(ds, src)  # set time_vals
                    except Exception:  # handle exception Exception:
                        time_vals = None  # set time_vals
                else:  # fallback branch
                    time_vals = _decode_cf_time_axis(ds, src)  # set time_vals
                    it = pick_time_index(time_vals, sim_time)  # set it
                _log_rain_time_usage(src, time_vals, it)  # execute statement
                field = np.asarray(da.isel({tdim: it}).values)  # set field
            else:  # fallback branch
                raise ValueError("Rain var must be 2D or 3D (time,y,x)")  # raise ValueError("Rain var must be 2D or 3D (time,y,x)")

            if field.shape != (H, W):  # check condition field.shape != (H, W):
                raise ValueError(f"Rain shape {field.shape} != domain shape {(H, W)}")  # raise ValueError(f"Rain shape {field.shape} != domain shape {(H, W)}")

        else:  # fallback branch
            raise ValueError(f"Unknown rain kind '{src.kind}' for '{src.name}'")  # raise ValueError(f"Unknown rain kind '{src.kind}' for '{src.name}'")

        total += src.weight * rain_to_step_mm(field, src.mode, dt_s)  # execute statement

    return total  # return total
