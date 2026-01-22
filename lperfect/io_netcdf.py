# -*- coding: utf-8 -*-
"""NetCDF I/O for outputs and restart state (rank0)."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import JSON for embedding config as provenance attribute.
import json  # import json

# Import typing primitives.
from typing import Any, Dict  # import typing import Any, Dict

# Import numpy.
import numpy as np  # import numpy as np

# Import xarray.
import xarray as xr  # import xarray as xr

# Import local time helper.
from .time_utils import utc_now_iso  # import .time_utils import utc_now_iso

# Import CF schema constants.
from .cf_schema import CF_CONVENTIONS, RAIN_TIME_UNITS  # import .cf_schema import CF_CONVENTIONS, RAIN_TIME_UNITS

# Import local Domain and Particles.
from .domain import Domain  # import .domain import Domain
from .particles import Particles  # import .particles import Particles

TIME_UNITS = RAIN_TIME_UNITS  # execute statement
OUTPUT_VARIABLES = (  # execute statement
    "flood_depth",
    "risk_index",
    "inundation_mask",
    "flood_depth_max",
    "inundation_mask_max",
)


def normalize_output_variables(cfg: Dict[str, Any]) -> set[str]:  # define function normalize_output_variables
    """Normalize output.variables into a validated set."""  # execute statement
    out_cfg = cfg.get("output", {})  # set out_cfg
    selected = out_cfg.get("variables", None)  # set selected
    if selected is None:  # check condition selected is None
        return set(OUTPUT_VARIABLES)  # return set(OUTPUT_VARIABLES)
    if isinstance(selected, str):  # check condition isinstance(selected, str)
        if selected.strip().lower() == "all":  # check condition selected is "all"
            return set(OUTPUT_VARIABLES)  # return set(OUTPUT_VARIABLES)
        raise ValueError("output.variables must be a list or 'all'.")  # raise ValueError
    if not isinstance(selected, list):  # check condition not list
        raise ValueError("output.variables must be a list of variable names.")  # raise ValueError
    if not selected:  # check condition not selected
        raise ValueError("output.variables must include at least one variable.")  # raise ValueError
    normalized: list[str] = []  # set normalized
    for item in selected:  # loop over selected
        if not isinstance(item, str):  # check condition not string
            raise ValueError("output.variables entries must be strings.")  # raise ValueError
        name = item.strip()  # set name
        if not name:  # check condition not name
            raise ValueError("output.variables entries must be non-empty strings.")  # raise ValueError
        if name.lower() == "all":  # check condition name is "all"
            return set(OUTPUT_VARIABLES)  # return set(OUTPUT_VARIABLES)
        normalized.append(name)  # execute statement
    requested = {name for name in normalized}  # set requested
    unknown = sorted(requested.difference(OUTPUT_VARIABLES))  # set unknown
    if unknown:  # check condition unknown
        raise ValueError(f"output.variables contains unknown entries: {', '.join(unknown)}")  # raise ValueError
    return set(normalized)  # return set(normalized)


def _coord_attrs(name: str) -> Dict[str, str]:  # define function _coord_attrs
    """Return CF-style coordinate attributes for lat/lon axes."""  # execute statement
    lname = name.lower()  # set lname
    if "lat" in lname:  # check condition "lat" in lname:
        return {"description": "Latitude", "long_name": "latitude", "units": "degrees_north"}  # return {"description": "Latitude", "long_name": "latitude", "units": "degrees_north"}
    if "lon" in lname:  # check condition "lon" in lname:
        return {"description": "Longitude", "long_name": "longitude", "units": "degrees_east"}  # return {"description": "Longitude", "long_name": "longitude", "units": "degrees_east"}
    return {}  # return {}


def write_results_netcdf_rank0(
    out_path: str,
    cfg: Dict[str, Any],
    dom: Domain,
    flood_depth_m: np.ndarray,
    risk_index: np.ndarray,
    inundation_mask: np.ndarray,
    flood_depth_max: np.ndarray,
    inundation_mask_max: np.ndarray,
    time_hours: float,
    mode: str = "w",
) -> None:  # execute statement
    """Write results in CF-friendly NetCDF (append-safe on time axis)."""  # execute statement
    out_cfg = cfg.get("output", {})  # set out_cfg
    inundation_threshold_m = float(out_cfg.get("inundation_threshold_m", 0.01))  # set inundation_threshold_m
    output_vars = normalize_output_variables(cfg)  # set output_vars

    ds = xr.Dataset()  # set ds
    fill_value = float(out_cfg.get("fill_value", -9999.0))  # set fill_value

    ds = ds.assign_coords(  # set ds
        {  # execute statement
            "time": xr.DataArray(np.array([time_hours], dtype=np.float64), dims=("time",), attrs={"description": "Time", "long_name": "time", "units": TIME_UNITS}),  # execute statement
            dom.x_name: xr.DataArray(dom.x_vals, dims=(dom.x_name,), attrs=_coord_attrs(dom.x_name)),  # execute statement
            dom.y_name: xr.DataArray(dom.y_vals, dims=(dom.y_name,), attrs=_coord_attrs(dom.y_name)),  # execute statement
        }  # execute statement
    )  # execute statement

    if "flood_depth" in output_vars:  # check condition "flood_depth" in output_vars
        ds["flood_depth"] = xr.DataArray(  # execute statement
            flood_depth_m.astype(np.float32, copy=False)[None, ...],  # execute statement
            dims=("time", dom.y_name, dom.x_name),  # set dims
            attrs={  # set attrs
                "standard_name": "water_depth",  # execute statement
                "long_name": "Flood water depth",  # execute statement
                "units": "m",  # execute statement
            },  # execute statement
        )  # execute statement

    if "risk_index" in output_vars:  # check condition "risk_index" in output_vars
        ds["risk_index"] = xr.DataArray(  # execute statement
            risk_index.astype(np.float32, copy=False)[None, ...],  # execute statement
            dims=("time", dom.y_name, dom.x_name),  # set dims
            attrs={  # set attrs
                "long_name": "Hydrogeological risk index",  # execute statement
                "units": "1",  # execute statement
            },  # execute statement
        )  # execute statement

    if "inundation_mask" in output_vars:  # check condition "inundation_mask" in output_vars
        ds["inundation_mask"] = xr.DataArray(  # execute statement
            inundation_mask.astype(np.int8, copy=False)[None, ...],  # execute statement
            dims=("time", dom.y_name, dom.x_name),  # set dims
            attrs={  # set attrs
                "long_name": "Inundation mask (1=inundated, 0=dry)",  # execute statement
                "units": "1",  # execute statement
                "flag_values": np.array([0, 1], dtype=np.int8),  # execute statement
                "flag_meanings": "dry inundated",  # execute statement
                "threshold_depth_m": inundation_threshold_m,  # execute statement
            },  # execute statement
        )  # execute statement

    if "flood_depth_max" in output_vars:  # check condition "flood_depth_max" in output_vars
        ds["flood_depth_max"] = xr.DataArray(  # execute statement
            flood_depth_max.astype(np.float32, copy=False),  # execute statement
            dims=(dom.y_name, dom.x_name),  # set dims
            attrs={  # set attrs
                "long_name": "Maximum flood water depth over simulation",  # execute statement
                "units": "m",  # execute statement
            },  # execute statement
        )  # execute statement

    if "inundation_mask_max" in output_vars:  # check condition "inundation_mask_max" in output_vars
        ds["inundation_mask_max"] = xr.DataArray(  # execute statement
            inundation_mask_max.astype(np.int8, copy=False),  # execute statement
            dims=(dom.y_name, dom.x_name),  # set dims
            attrs={  # set attrs
                "long_name": "Ever inundated during simulation",  # execute statement
                "units": "1",  # execute statement
                "flag_values": np.array([0, 1], dtype=np.int8),  # execute statement
                "flag_meanings": "never_inundated inundated",  # execute statement
                "threshold_depth_m": inundation_threshold_m,  # execute statement
            },  # execute statement
        )  # execute statement

    if dom.grid_mapping_name and dom.grid_mapping_attrs:  # check condition dom.grid_mapping_name and dom.grid_mapping_attrs:
        gm = dom.grid_mapping_name  # set gm
        ds[gm] = xr.DataArray(0, attrs=dom.grid_mapping_attrs)  # execute statement
        if "flood_depth" in ds:  # check condition "flood_depth" in ds
            ds["flood_depth"].attrs["grid_mapping"] = gm  # execute statement
        if "risk_index" in ds:  # check condition "risk_index" in ds
            ds["risk_index"].attrs["grid_mapping"] = gm  # execute statement
        if "inundation_mask" in ds:  # check condition "inundation_mask" in ds
            ds["inundation_mask"].attrs["grid_mapping"] = gm  # execute statement
        if "flood_depth_max" in ds:  # check condition "flood_depth_max" in ds
            ds["flood_depth_max"].attrs["grid_mapping"] = gm  # execute statement
        if "inundation_mask_max" in ds:  # check condition "inundation_mask_max" in ds
            ds["inundation_mask_max"].attrs["grid_mapping"] = gm  # execute statement

    ds.attrs["title"] = out_cfg.get("title", "LPERFECT flood depth + hydrogeological risk index")  # execute statement
    ds.attrs["institution"] = out_cfg.get("institution", "")  # execute statement
    ds.attrs["source"] = "LPERFECT"  # execute statement
    ds.attrs["history"] = f"{utc_now_iso()}: results written by LPERFECT"  # execute statement
    ds.attrs["Conventions"] = out_cfg.get("Conventions", CF_CONVENTIONS)  # execute statement
    ds.attrs["lperfect_config_json"] = json.dumps(cfg, separators=(",", ":"), sort_keys=True)  # execute statement
    ds.attrs["inundation_threshold_m"] = inundation_threshold_m  # execute statement

    encoding: Dict[str, Dict[str, Any]] = {}  # set encoding
    for name in ("flood_depth", "risk_index", "flood_depth_max"):  # loop over float outputs
        if name in ds:  # check condition name in ds
            encoding[name] = {"_FillValue": fill_value}  # execute statement
    to_netcdf_kwargs: Dict[str, Any] = {"encoding": encoding, "mode": mode}  # set to_netcdf_kwargs
    if mode == "a":  # check condition mode == "a":
        to_netcdf_kwargs["unlimited_dims"] = ("time",)  # execute statement
    ds.to_netcdf(out_path, **to_netcdf_kwargs)  # execute statement


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
        attrs={"long_name": "cumulative_precipitation", "units": "mm"},  # execute statement
    )  # execute statement
    ds["Q_cum_mm"] = xr.DataArray(  # execute statement
        Q_cum_mm_full.astype(np.float64),  # execute statement
        dims=(dom.y_name, dom.x_name),  # set dims
        attrs={"long_name": "cumulative_runoff_depth", "units": "mm"},  # execute statement
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
    ds.attrs["Conventions"] = cfg.get("output", {}).get("Conventions", CF_CONVENTIONS)  # execute statement
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
