# -*- coding: utf-8 -*-
"""Domain reading and broadcasting (NetCDF-only)."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import typing primitives.
from typing import Any, Dict, Optional  # import typing import Any, Dict, Optional

# Import dataclass for structured domain object.
from dataclasses import dataclass  # import dataclasses import dataclass

# Import numpy for arrays.
import numpy as np  # import numpy as np

# Import xarray for NetCDF reading.
import xarray as xr  # import xarray as xr


@dataclass  # apply decorator
class Domain:  # define class Domain
    """Domain data required by the model."""  # execute statement
    dem: np.ndarray  # execute statement
    d8: np.ndarray  # execute statement
    cn: np.ndarray  # execute statement
    channel_mask: Optional[np.ndarray]  # execute statement
    active_mask: np.ndarray  # execute statement
    x_name: str  # execute statement
    y_name: str  # execute statement
    x_vals: np.ndarray  # execute statement
    y_vals: np.ndarray  # execute statement
    grid_mapping_name: Optional[str]  # execute statement
    grid_mapping_attrs: Dict[str, Any]  # execute statement
    cell_area_m2: float | np.ndarray  # execute statement


def cell_area_m2_from_domain(ds: xr.Dataset, x_name: str, y_name: str) -> float | np.ndarray:  # define function cell_area_m2_from_domain
    """Estimate cell area from coordinate spacing and units.

    If x/y units look like meters -> constant area dx*dy.
    If x/y are degrees -> compute per-row geodesic area (WGS84) when possible.
    """
    # Extract coordinate vectors.
    x = np.asarray(ds[x_name].values)  # set x
    y = np.asarray(ds[y_name].values)  # set y
    # Estimate spacings robustly via median.
    dx = float(np.median(np.abs(np.diff(x))))  # set dx
    dy = float(np.median(np.abs(np.diff(y))))  # set dy
    # Read units (if any).
    xu = str(ds[x_name].attrs.get("units", "")).lower()  # set xu
    yu = str(ds[y_name].attrs.get("units", "")).lower()  # set yu

    # Projected meters case -> constant area.
    if ("m" in xu) and ("m" in yu):  # check condition ("m" in xu) and ("m" in yu):
        return float(dx * dy)  # return float(dx * dy)

    # Geographic degrees case -> compute per-row areas.
    H = int(ds.sizes[y_name])  # set H
    W = int(ds.sizes[x_name])  # set W

    # Try pyproj Geod for accurate ellipsoidal area.
    try:  # start exception handling
        from pyproj import Geod  # import pyproj import Geod
        geod = Geod(ellps="WGS84")  # set geod
        # Build one cell polygon width using dx.
        lon_left = float(x.min())  # set lon_left
        lon_right = float(lon_left + (dx if x[1] > x[0] else -dx))  # set lon_right
        # Allocate per-row areas.
        areas_row = np.zeros(H, dtype=np.float64)  # set areas_row
        # Loop rows (y dimension).
        for i in range(H):  # loop over i in range(H):
            lat_top = float(y[i])  # set lat_top
            lat_bot = float(lat_top - dy if y[1] < y[0] else lat_top + dy)  # set lat_bot
            lons = [lon_left, lon_right, lon_right, lon_left]  # set lons
            lats = [lat_top, lat_top, lat_bot, lat_bot]  # set lats
            poly_area, _ = geod.polygon_area_perimeter(lons, lats)  # set poly_area, _
            areas_row[i] = abs(poly_area)  # execute statement
        # Expand to full grid.
        return np.repeat(areas_row[:, None], W, axis=1)  # return np.repeat(areas_row[:, None], W, axis=1)
    except Exception:  # handle exception Exception:
        # Fallback: simple spherical approximation.
        R = 6371000.0  # set R
        dlon = np.deg2rad(abs(dx))  # set dlon
        areas_row = np.zeros(H, dtype=np.float64)  # set areas_row
        for i in range(H):  # loop over i in range(H):
            lat_top = np.deg2rad(float(y[i]))  # set lat_top
            lat_bot = np.deg2rad(float(y[i] - dy if y[1] < y[0] else y[i] + dy))  # set lat_bot
            areas_row[i] = abs((R * R) * dlon * (np.sin(lat_bot) - np.sin(lat_top)))  # execute statement
        return np.repeat(areas_row[:, None], W, axis=1)  # return np.repeat(areas_row[:, None], W, axis=1)


def read_domain_netcdf_rank0(cfg: Dict[str, Any]) -> Domain:  # define function read_domain_netcdf_rank0
    """Read domain NetCDF on rank 0."""  # execute statement
    # Extract domain configuration.
    dom_cfg = cfg["domain"]  # set dom_cfg
    path = dom_cfg["domain_nc"]  # set path
    varmap = dom_cfg.get("varmap", {})  # set varmap

    # Open dataset.
    ds = xr.open_dataset(path)  # set ds

    # Resolve variable names.
    dem_name = varmap.get("dem", "dem")  # set dem_name
    d8_name = varmap.get("d8", "d8")  # set d8_name
    cn_name = varmap.get("cn", "cn")  # set cn_name
    ch_name = varmap.get("channel_mask", "channel_mask")  # set ch_name

    # Read arrays.
    dem = np.asarray(ds[dem_name].values).astype(np.float32)  # set dem
    active_mask = np.isfinite(dem)  # set active_mask

    d8_da = ds[d8_name]  # set d8_da
    d8_raw = np.asarray(d8_da.values)  # set d8_raw
    d8_fill = d8_da.attrs.get("_FillValue", None)  # set d8_fill
    d8_mask = ~np.isfinite(d8_raw)  # set d8_mask
    if d8_fill is not None:  # check condition d8_fill is not None:
        d8_mask |= d8_raw == d8_fill  # execute statement
    d8_clean = np.where(d8_mask, 0, d8_raw)  # set d8_clean
    d8 = np.where(active_mask, d8_clean, 0).astype(np.int32)  # set d8
    cn = np.asarray(ds[cn_name].values).astype(np.float32)  # set cn

    # Optional channel mask.
    channel_mask = None  # set channel_mask
    if ch_name in ds:  # check condition ch_name in ds:
        channel_mask = (np.asarray(ds[ch_name].values) > 0)  # set channel_mask

    # Resolve coordinate names.
    x_name = varmap.get("x", "x")  # set x_name
    y_name = varmap.get("y", "y")  # set y_name

    # Read coordinates from coords or variables.
    x_vals = np.asarray(ds.coords[x_name].values if x_name in ds.coords else ds[x_name].values)  # set x_vals
    y_vals = np.asarray(ds.coords[y_name].values if y_name in ds.coords else ds[y_name].values)  # set y_vals

    # Estimate cell area.
    cell_area_m2 = cell_area_m2_from_domain(ds, x_name=x_name, y_name=y_name)  # set cell_area_m2
    if isinstance(cell_area_m2, np.ndarray):
        cell_area_m2 = cell_area_m2.astype(np.float32)

    # Preserve CF grid mapping if present.
    gm_name = ds[dem_name].attrs.get("grid_mapping", None)  # set gm_name
    gm_attrs: Dict[str, Any] = {}  # execute statement
    if gm_name and gm_name in ds:  # check condition gm_name and gm_name in ds:
        gm_attrs = dict(ds[gm_name].attrs)  # set gm_attrs

    # Close dataset.
    ds.close()  # execute statement

    # Clean CN outside active cells.
    cn = np.where(active_mask & np.isfinite(cn), cn, 0.0)  # set cn

    # Apply active mask to channel mask if present.
    if channel_mask is not None:  # check condition channel_mask is not None:
        channel_mask = channel_mask & active_mask  # set channel_mask

    # Return domain object.
    return Domain(  # return Domain(
        dem=dem,  # set dem
        d8=d8,  # set d8
        cn=cn,  # set cn
        channel_mask=channel_mask,  # set channel_mask
        active_mask=active_mask,  # set active_mask
        x_name=x_name,  # set x_name
        y_name=y_name,  # set y_name
        x_vals=x_vals,  # set x_vals
        y_vals=y_vals,  # set y_vals
        grid_mapping_name=gm_name,  # set grid_mapping_name
        grid_mapping_attrs=gm_attrs,  # set grid_mapping_attrs
        cell_area_m2=cell_area_m2,  # set cell_area_m2
    )  # execute statement


def bcast_domain(comm, dom0: Optional[Domain]) -> Domain:  # define function bcast_domain
    """Broadcast Domain from rank0 to all ranks."""  # execute statement
    # Rank id.
    rank = comm.Get_rank()  # set rank

    # Prepare metadata dict on root.
    if rank == 0:  # check condition rank == 0:
        meta = {  # set meta
            "shape": dom0.dem.shape,  # execute statement
            "x_name": dom0.x_name,  # execute statement
            "y_name": dom0.y_name,  # execute statement
            "x_vals": dom0.x_vals,  # execute statement
            "y_vals": dom0.y_vals,  # execute statement
            "grid_mapping_name": dom0.grid_mapping_name,  # execute statement
            "grid_mapping_attrs": dom0.grid_mapping_attrs,  # execute statement
            "has_channel_mask": dom0.channel_mask is not None,  # execute statement
            "cell_area_is_scalar": bool(np.isscalar(dom0.cell_area_m2)),  # execute statement
            "cell_area_scalar": float(dom0.cell_area_m2) if np.isscalar(dom0.cell_area_m2) else None,  # execute statement
        }  # execute statement
    else:  # fallback branch
        meta = None  # set meta

    # Broadcast metadata (pickle-based).
    meta = comm.bcast(meta, root=0)  # set meta

    # Extract shape.
    H, W = meta["shape"]  # set H, W

    # Allocate or reuse arrays.
    if rank != 0:  # check condition rank != 0:
        dem = np.empty((H, W), dtype=np.float64)  # set dem
        d8 = np.empty((H, W), dtype=np.int32)  # set d8
        cn = np.empty((H, W), dtype=np.float64)  # set cn
        active = np.empty((H, W), dtype=np.bool_)  # set active
    else:  # fallback branch
        dem = dom0.dem.astype(np.float64)  # set dem
        d8 = dom0.d8.astype(np.int32)  # set d8
        cn = dom0.cn.astype(np.float64)  # set cn
        active = dom0.active_mask.astype(np.bool_)  # set active

    # Broadcast arrays.
    comm.Bcast(dem, root=0)  # execute statement
    comm.Bcast(d8, root=0)  # execute statement
    comm.Bcast(cn, root=0)  # execute statement
    comm.Bcast(active, root=0)  # execute statement

    # Optional channel mask.
    channel_mask = None  # set channel_mask
    if meta["has_channel_mask"]:  # check condition meta["has_channel_mask"]:
        if rank != 0:  # check condition rank != 0:
            cm = np.empty((H, W), dtype=np.bool_)  # set cm
        else:  # fallback branch
            cm = dom0.channel_mask.astype(np.bool_)  # set cm
        comm.Bcast(cm, root=0)  # execute statement
        channel_mask = cm  # set channel_mask

    # Cell area.
    if meta["cell_area_is_scalar"]:  # check condition meta["cell_area_is_scalar"]:
        cell_area_m2: float | np.ndarray = float(meta["cell_area_scalar"])  # execute statement
    else:  # fallback branch
        if rank != 0:  # check condition rank != 0:
            ca = np.empty((H, W), dtype=np.float64)  # set ca
        else:  # fallback branch
            ca = dom0.cell_area_m2.astype(np.float64)  # set ca
        comm.Bcast(ca, root=0)  # execute statement
        cell_area_m2 = ca  # set cell_area_m2

    # Return domain.
    return Domain(  # return Domain(
        dem=dem,  # set dem
        d8=d8,  # set d8
        cn=cn,  # set cn
        channel_mask=channel_mask,  # set channel_mask
        active_mask=active,  # set active_mask
        x_name=meta["x_name"],  # set x_name
        y_name=meta["y_name"],  # set y_name
        x_vals=np.asarray(meta["x_vals"]),  # set x_vals
        y_vals=np.asarray(meta["y_vals"]),  # set y_vals
        grid_mapping_name=meta["grid_mapping_name"],  # set grid_mapping_name
        grid_mapping_attrs=dict(meta["grid_mapping_attrs"]),  # set grid_mapping_attrs
        cell_area_m2=cell_area_m2,  # set cell_area_m2
    )  # execute statement
