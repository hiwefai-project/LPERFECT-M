# -*- coding: utf-8 -*-
"""NetCDF I/O for outputs and restart state (rank0)."""  # execute statement

# Import JSON for embedding config as provenance attribute.
import json  # import json

# Import datetime utilities for time encoding.
from datetime import datetime, timezone  # import datetime, timezone

# Import typing primitives.
from typing import Any, Dict  # import typing import Any, Dict

# Import numpy.
import numpy as np  # import numpy as np

# Import xarray.
import xarray as xr  # import xarray as xr

# Import local time helper.
from .time_utils import utc_now_iso  # import .time_utils import utc_now_iso

# Import local Domain and Particles.
from .domain import Domain  # import .domain import Domain
from .particles import Particles  # import .particles import Particles

TIME_UNITS = "hours since 1900-01-01 00:00:0.0"  # execute statement


def _parse_iso_datetime(value: str | None) -> datetime:  # define function _parse_iso_datetime
    """Parse ISO-8601 string into UTC datetime, fallback to now."""  # execute statement
    if not value:  # check condition not value:
        return datetime.now(timezone.utc)  # return datetime.now(timezone.utc)
    text = value.strip()  # set text
    if text.endswith("Z"):  # check condition text.endswith("Z"):
        text = text[:-1] + "+00:00"  # set text
    dt = datetime.fromisoformat(text)  # set dt
    if dt.tzinfo is None:  # check condition dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)  # set dt
    return dt.astimezone(timezone.utc)  # return dt.astimezone(timezone.utc)


def _time_to_hours_since_1900(dt: datetime) -> int:  # define function _time_to_hours_since_1900
    """Convert datetime to integer hours since 1900-01-01 UTC."""  # execute statement
    base = datetime(1900, 1, 1, tzinfo=timezone.utc)  # set base
    hours = (dt - base).total_seconds() / 3600.0  # set hours
    return int(round(hours))  # return int(round(hours))


def _coord_attrs(name: str) -> Dict[str, str]:  # define function _coord_attrs
    """Return CF-style coordinate attributes for lat/lon axes."""  # execute statement
    lname = name.lower()  # set lname
    if "lat" in lname:  # check condition "lat" in lname:
        return {"description": "Latitude", "long_name": "latitude", "units": "degrees_north"}  # return {"description": "Latitude", "long_name": "latitude", "units": "degrees_north"}
    if "lon" in lname:  # check condition "lon" in lname:
        return {"description": "Longitude", "long_name": "longitude", "units": "degrees_east"}  # return {"description": "Longitude", "long_name": "longitude", "units": "degrees_east"}
    return {}  # return {}


def write_results_netcdf_rank0(out_path: str, cfg: Dict[str, Any], dom: Domain,  # define function write_results_netcdf_rank0
                              flood_depth_m: np.ndarray, risk_index: np.ndarray) -> None:  # execute statement
    """Write final results in CF-friendly NetCDF."""  # execute statement
    out_cfg = cfg.get("output", {})  # set out_cfg

    ds = xr.Dataset()  # set ds
    time_value = _time_to_hours_since_1900(_parse_iso_datetime(cfg.get("model", {}).get("start_time")))  # execute statement
    fill_value = float(out_cfg.get("fill_value", -9999.0))  # set fill_value

    ds = ds.assign_coords({  # set ds
        "time": xr.DataArray(np.array([time_value], dtype=np.int32), dims=("time",), attrs={"description": "Time", "long_name": "time", "units": TIME_UNITS}),  # execute statement
        dom.x_name: xr.DataArray(dom.x_vals, dims=(dom.x_name,), attrs=_coord_attrs(dom.x_name)),  # execute statement
        dom.y_name: xr.DataArray(dom.y_vals, dims=(dom.y_name,), attrs=_coord_attrs(dom.y_name)),  # execute statement
    })  # execute statement

    ds["flood_depth"] = xr.DataArray(  # execute statement
        flood_depth_m.astype(np.float32)[None, ...],  # execute statement
        dims=("time", dom.y_name, dom.x_name),  # set dims
        attrs={  # set attrs
            "standard_name": "water_depth",  # execute statement
            "long_name": "flooded_water_depth",  # execute statement
            "units": "m",  # execute statement
            "_FillValue": fill_value,  # execute statement
        },  # execute statement
    )  # execute statement

    ds["risk_index"] = xr.DataArray(  # execute statement
        risk_index.astype(np.float32)[None, ...],  # execute statement
        dims=("time", dom.y_name, dom.x_name),  # set dims
        attrs={  # set attrs
            "long_name": "hydrogeological_risk_index",  # execute statement
            "units": "1",  # execute statement
            "_FillValue": fill_value,  # execute statement
        },  # execute statement
    )  # execute statement

    if dom.grid_mapping_name and dom.grid_mapping_attrs:  # check condition dom.grid_mapping_name and dom.grid_mapping_attrs:
        gm = dom.grid_mapping_name  # set gm
        ds[gm] = xr.DataArray(0, attrs=dom.grid_mapping_attrs)  # execute statement
        ds["flood_depth"].attrs["grid_mapping"] = gm  # execute statement
        ds["risk_index"].attrs["grid_mapping"] = gm  # execute statement

    ds.attrs["title"] = out_cfg.get("title", "LPERFECT flood depth + hydrogeological risk index")  # execute statement
    ds.attrs["institution"] = out_cfg.get("institution", "")  # execute statement
    ds.attrs["source"] = "LPERFECT"  # execute statement
    ds.attrs["history"] = f"{utc_now_iso()}: results written by LPERFECT"  # execute statement
    ds.attrs["Conventions"] = out_cfg.get("Conventions", "CF-1.10")  # execute statement
    ds.attrs["lperfect_config_json"] = json.dumps(cfg, separators=(",", ":"), sort_keys=True)  # execute statement

    ds.to_netcdf(out_path, encoding={"flood_depth": {"_FillValue": fill_value}, "risk_index": {"_FillValue": fill_value}})  # execute statement


def save_restart_netcdf_rank0(out_path: str, cfg: Dict[str, Any], dom: Domain,  # define function save_restart_netcdf_rank0
                             elapsed_s: float, cum_rain_vol_m3: float, cum_runoff_vol_m3: float, cum_outflow_vol_m3: float,  # execute statement
                             P_cum_mm_full: np.ndarray, Q_cum_mm_full: np.ndarray, particles_all: Particles) -> None:  # execute statement
    """Save restart state to NetCDF."""  # execute statement
    ds = xr.Dataset()  # set ds
    fill_value = float(cfg.get("restart", {}).get("fill_value", -9999.0))  # set fill_value

    ds = ds.assign_coords({  # set ds
        dom.x_name: xr.DataArray(dom.x_vals, dims=(dom.x_name,), attrs=_coord_attrs(dom.x_name)),  # execute statement
        dom.y_name: xr.DataArray(dom.y_vals, dims=(dom.y_name,), attrs=_coord_attrs(dom.y_name)),  # execute statement
        "particle": xr.DataArray(np.arange(particles_all.r.size, dtype=np.int64), dims=("particle",)),  # execute statement
    })  # execute statement

    ds["P_cum_mm"] = xr.DataArray(  # execute statement
        P_cum_mm_full.astype(np.float64),  # execute statement
        dims=(dom.y_name, dom.x_name),  # set dims
        attrs={"long_name": "cumulative_precipitation", "units": "mm", "_FillValue": fill_value},  # execute statement
    )  # execute statement
    ds["Q_cum_mm"] = xr.DataArray(  # execute statement
        Q_cum_mm_full.astype(np.float64),  # execute statement
        dims=(dom.y_name, dom.x_name),  # set dims
        attrs={"long_name": "cumulative_runoff_depth", "units": "mm", "_FillValue": fill_value},  # execute statement
    )  # execute statement

    ds["particle_r"] = xr.DataArray(particles_all.r.astype(np.int32), dims=("particle",), attrs={"units": "1"})  # execute statement
    ds["particle_c"] = xr.DataArray(particles_all.c.astype(np.int32), dims=("particle",), attrs={"units": "1"})  # execute statement
    ds["particle_vol"] = xr.DataArray(particles_all.vol.astype(np.float64), dims=("particle",), attrs={"units": "m3"})  # execute statement
    ds["particle_tau"] = xr.DataArray(particles_all.tau.astype(np.float64), dims=("particle",), attrs={"units": "s"})  # execute statement

    ds["elapsed_s"] = xr.DataArray(np.array(elapsed_s, dtype=np.float64), attrs={"units": "s"})  # execute statement
    ds["cum_rain_vol_m3"] = xr.DataArray(np.array(cum_rain_vol_m3, dtype=np.float64), attrs={"units": "m3"})  # execute statement
    ds["cum_runoff_vol_m3"] = xr.DataArray(np.array(cum_runoff_vol_m3, dtype=np.float64), attrs={"units": "m3"})  # execute statement
    ds["cum_outflow_vol_m3"] = xr.DataArray(np.array(cum_outflow_vol_m3, dtype=np.float64), attrs={"units": "m3"})  # execute statement

    if dom.grid_mapping_name and dom.grid_mapping_attrs:  # check condition dom.grid_mapping_name and dom.grid_mapping_attrs:
        gm = dom.grid_mapping_name  # set gm
        ds[gm] = xr.DataArray(0, attrs=dom.grid_mapping_attrs)  # execute statement
        ds["P_cum_mm"].attrs["grid_mapping"] = gm  # execute statement
        ds["Q_cum_mm"].attrs["grid_mapping"] = gm  # execute statement

    ds.attrs["title"] = "LPERFECT restart"  # execute statement
    ds.attrs["source"] = "LPERFECT"  # execute statement
    ds.attrs["history"] = f"{utc_now_iso()}: restart written by LPERFECT"  # execute statement
    ds.attrs["Conventions"] = cfg.get("output", {}).get("Conventions", "CF-1.10")  # execute statement
    ds.attrs["lperfect_config_json"] = json.dumps(cfg, separators=(",", ":"), sort_keys=True)  # execute statement

    ds.to_netcdf(out_path, encoding={"P_cum_mm": {"_FillValue": fill_value}, "Q_cum_mm": {"_FillValue": fill_value}})  # execute statement


def load_restart_netcdf_rank0(path: str) -> Dict[str, Any]:  # define function load_restart_netcdf_rank0
    """Load restart NetCDF and return state dict."""  # execute statement
    ds = xr.open_dataset(path)  # set ds

    out = {  # set out
        "P_cum_mm": np.asarray(ds["P_cum_mm"].values).astype(np.float64),  # execute statement
        "Q_cum_mm": np.asarray(ds["Q_cum_mm"].values).astype(np.float64),  # execute statement
        "r": np.asarray(ds["particle_r"].values).astype(np.int32),  # execute statement
        "c": np.asarray(ds["particle_c"].values).astype(np.int32),  # execute statement
        "vol": np.asarray(ds["particle_vol"].values).astype(np.float64),  # execute statement
        "tau": np.asarray(ds["particle_tau"].values).astype(np.float64),  # execute statement
        "elapsed_s": float(ds["elapsed_s"].values),  # execute statement
        "cum_rain_vol_m3": float(ds["cum_rain_vol_m3"].values),  # execute statement
        "cum_runoff_vol_m3": float(ds["cum_runoff_vol_m3"].values),  # execute statement
        "cum_outflow_vol_m3": float(ds["cum_outflow_vol_m3"].values),  # execute statement
    }  # execute statement

    ds.close()  # execute statement
    return out  # return out
