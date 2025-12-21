#!/usr/bin/env python3
"""Create a domain NetCDF from DEM/CN/D8/river mask inputs."""
from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable, Optional

import numpy as np
import xarray as xr


@dataclass
class GridSpec:
    x_name: str
    y_name: str
    x: np.ndarray
    y: np.ndarray


D8_CODES = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.int32)
D8_DIRS = [
    (0, 1),    # E
    (1, 1),    # SE
    (1, 0),    # S
    (1, -1),   # SW
    (0, -1),   # W
    (-1, -1),  # NW
    (-1, 0),   # N
    (-1, 1),   # NE
]


def _normalize_da(da: xr.DataArray, x_name: str, y_name: str) -> xr.DataArray:
    if len(da.dims) != 2:
        raise ValueError(f"Expected 2D variable; got dims={da.dims}")
    if da.dims != (y_name, x_name):
        da = da.rename({da.dims[0]: y_name, da.dims[1]: x_name})
    return da


def _slice_for_bbox(coord: np.ndarray, minv: float, maxv: float) -> slice:
    ascending = coord[0] < coord[-1]
    if ascending:
        return slice(minv, maxv)
    return slice(maxv, minv)


def _build_axis(minv: float, maxv: float, step: float, ascending: bool) -> np.ndarray:
    if step <= 0:
        raise ValueError("Resolution must be positive")
    step = float(step)
    if not ascending:
        step = -step
    start = minv if ascending else maxv
    end = maxv if ascending else minv
    return np.arange(start, end + step * 0.5, step)


def _grid_from_da(da: xr.DataArray, x_name: str, y_name: str) -> GridSpec:
    x = np.asarray(da[x_name].values)
    y = np.asarray(da[y_name].values)
    return GridSpec(x_name=x_name, y_name=y_name, x=x, y=y)


def _apply_bbox(da: xr.DataArray, grid: GridSpec, bbox: Optional[Iterable[float]]) -> xr.DataArray:
    if bbox is None:
        return da
    xmin, ymin, xmax, ymax = bbox
    return da.sel(
        {
            grid.x_name: _slice_for_bbox(grid.x, xmin, xmax),
            grid.y_name: _slice_for_bbox(grid.y, ymin, ymax),
        }
    )


def _apply_resolution(da: xr.DataArray, grid: GridSpec, res: Optional[Iterable[float]], method: str) -> xr.DataArray:
    if res is None:
        return da
    res_vals = list(res)
    if len(res_vals) == 1:
        dx = dy = res_vals[0]
    elif len(res_vals) == 2:
        dx, dy = res_vals
    else:
        raise ValueError("Resolution expects 1 or 2 values")
    xmin, xmax = float(np.min(grid.x)), float(np.max(grid.x))
    ymin, ymax = float(np.min(grid.y)), float(np.max(grid.y))
    x_asc = grid.x[0] < grid.x[-1]
    y_asc = grid.y[0] < grid.y[-1]
    new_x = _build_axis(xmin, xmax, dx, x_asc)
    new_y = _build_axis(ymin, ymax, dy, y_asc)
    return da.interp({grid.x_name: new_x, grid.y_name: new_y}, method=method)


def _regrid_to_target(da: xr.DataArray, target: GridSpec, method: str) -> xr.DataArray:
    return da.interp({target.x_name: target.x, target.y_name: target.y}, method=method)


def _load_var(path: str, var_name: str, x_name: str, y_name: str) -> tuple[xr.DataArray, xr.Dataset]:
    ds = xr.open_dataset(path)
    if var_name not in ds:
        raise KeyError(f"Variable '{var_name}' not found in {path}")
    da = _normalize_da(ds[var_name], x_name=x_name, y_name=y_name)
    return da, ds


def _compute_d8(dem: np.ndarray, dx: float, dy: float) -> np.ndarray:
    dem = dem.astype(np.float64)
    nrows, ncols = dem.shape
    pad = np.pad(dem, 1, mode="constant", constant_values=np.nan)
    slopes = []
    dist_diag = float(np.hypot(dx, dy))
    dist_card = float(np.mean([abs(dx), abs(dy)]))
    distances = [dist_card, dist_diag, dist_card, dist_diag, dist_card, dist_diag, dist_card, dist_diag]
    for (dr, dc), dist in zip(D8_DIRS, distances):
        neigh = pad[1 + dr:1 + dr + nrows, 1 + dc:1 + dc + ncols]
        slope = (dem - neigh) / dist
        slope = np.where(np.isfinite(dem) & np.isfinite(neigh), slope, -np.inf)
        slopes.append(slope)
    slope_stack = np.stack(slopes, axis=0)
    max_idx = np.argmax(slope_stack, axis=0)
    max_slope = np.max(slope_stack, axis=0)
    d8 = np.zeros_like(dem, dtype=np.int32)
    valid = max_slope > 0
    for idx, code in enumerate(D8_CODES):
        mask = valid & (max_idx == idx)
        d8[mask] = code
    return d8


def build_domain(
    dem_path: str,
    output_path: str,
    cn_path: Optional[str],
    d8_path: Optional[str],
    mask_path: Optional[str],
    dem_var: str,
    cn_var: str,
    d8_var: str,
    mask_var: str,
    x_name: str,
    y_name: str,
    bbox: Optional[Iterable[float]],
    resolution: Optional[Iterable[float]],
) -> None:
    dem_da, dem_ds = _load_var(dem_path, dem_var, x_name=x_name, y_name=y_name)

    dem_da = _apply_bbox(dem_da, _grid_from_da(dem_da, x_name, y_name), bbox)
    dem_da = _apply_resolution(dem_da, _grid_from_da(dem_da, x_name, y_name), resolution, method="linear")

    grid = _grid_from_da(dem_da, x_name, y_name)
    dx = float(np.median(np.abs(np.diff(grid.x))))
    dy = float(np.median(np.abs(np.diff(grid.y))))

    if d8_path:
        d8_da, d8_ds = _load_var(d8_path, d8_var, x_name=x_name, y_name=y_name)
        d8_da = _apply_bbox(d8_da, _grid_from_da(d8_da, x_name, y_name), bbox)
        d8_da = _regrid_to_target(d8_da, grid, method="nearest")
        d8 = np.asarray(d8_da.values).astype(np.int32)
        d8_ds.close()
    else:
        d8 = _compute_d8(np.asarray(dem_da.values), dx=dx, dy=dy)

    if cn_path:
        cn_da, cn_ds = _load_var(cn_path, cn_var, x_name=x_name, y_name=y_name)
        cn_da = _apply_bbox(cn_da, _grid_from_da(cn_da, x_name, y_name), bbox)
        cn_da = _regrid_to_target(cn_da, grid, method="nearest")
        cn = np.asarray(cn_da.values).astype(np.float64)
        cn_ds.close()
    else:
        cn = np.zeros_like(dem_da.values, dtype=np.float64)

    channel_mask = None
    if mask_path:
        mask_da, mask_ds = _load_var(mask_path, mask_var, x_name=x_name, y_name=y_name)
        mask_da = _apply_bbox(mask_da, _grid_from_da(mask_da, x_name, y_name), bbox)
        mask_da = _regrid_to_target(mask_da, grid, method="nearest")
        channel_mask = (np.asarray(mask_da.values) > 0).astype(np.int8)
        mask_ds.close()

    ds_out = xr.Dataset()
    ds_out = ds_out.assign_coords(
        {
            x_name: xr.DataArray(grid.x, dims=(x_name,)),
            y_name: xr.DataArray(grid.y, dims=(y_name,)),
        }
    )

    dem_attrs = dict(dem_da.attrs)
    ds_out[dem_var] = xr.DataArray(
        np.asarray(dem_da.values, dtype=np.float64),
        dims=(y_name, x_name),
        attrs={
            "standard_name": dem_attrs.get("standard_name", "surface_altitude"),
            "long_name": dem_attrs.get("long_name", "digital_elevation_model"),
            "units": dem_attrs.get("units", "m"),
        },
    )

    ds_out[d8_var] = xr.DataArray(
        d8.astype(np.int32),
        dims=(y_name, x_name),
        attrs={
            "long_name": "D8_flow_direction",
            "flag_values": "1 2 4 8 16 32 64 128",
            "flag_meanings": "E SE S SW W NW N NE",
            "comment": "ESRI D8 encoding (see LPERFECT model.encoding).",
        },
    )

    ds_out[cn_var] = xr.DataArray(
        cn.astype(np.float64),
        dims=(y_name, x_name),
        attrs={
            "long_name": "SCS_curve_number",
            "units": "1",
            "valid_min": 0.0,
            "valid_max": 100.0,
        },
    )

    if channel_mask is not None:
        ds_out[mask_var] = xr.DataArray(
            channel_mask,
            dims=(y_name, x_name),
            attrs={
                "long_name": "channel_mask",
                "units": "1",
                "flag_values": "0 1",
                "flag_meanings": "no_channel channel",
            },
        )

    grid_mapping = dem_attrs.get("grid_mapping")
    if grid_mapping and grid_mapping in dem_ds:
        ds_out[grid_mapping] = xr.DataArray(0, attrs=dict(dem_ds[grid_mapping].attrs))
        ds_out[dem_var].attrs["grid_mapping"] = grid_mapping

    ds_out.attrs.update(
        {
            "title": "LPERFECT simulation domain",
            "source": "Prepared for LPERFECT (DEM + D8 + CN)",
            "Conventions": "CF-1.10",
        }
    )

    ds_out.to_netcdf(output_path)
    dem_ds.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dem", required=True, help="Path to DEM NetCDF")
    parser.add_argument("--cn", help="Path to CN NetCDF (optional)")
    parser.add_argument("--d8", help="Path to D8 NetCDF (optional)")
    parser.add_argument("--mask", help="Path to channel/river mask NetCDF (optional)")
    parser.add_argument("--dem-var", default="dem", help="DEM variable name")
    parser.add_argument("--cn-var", default="cn", help="CN variable name")
    parser.add_argument("--d8-var", default="d8", help="D8 variable name")
    parser.add_argument("--mask-var", default="channel_mask", help="Mask variable name")
    parser.add_argument("--x-name", default="x", help="X/longitude coordinate name")
    parser.add_argument("--y-name", default="y", help="Y/latitude coordinate name")
    parser.add_argument("--bbox", nargs=4, type=float, metavar=("MINX", "MINY", "MAXX", "MAXY"))
    parser.add_argument(
        "--resolution",
        nargs="+",
        type=float,
        metavar=("DX", "DY"),
        help="Target resolution (dx [dy]) in degrees or meters",
    )
    parser.add_argument("--output", required=True, help="Output NetCDF path")
    args = parser.parse_args()

    build_domain(
        dem_path=args.dem,
        output_path=args.output,
        cn_path=args.cn,
        d8_path=args.d8,
        mask_path=args.mask,
        dem_var=args.dem_var,
        cn_var=args.cn_var,
        d8_var=args.d8_var,
        mask_var=args.mask_var,
        x_name=args.x_name,
        y_name=args.y_name,
        bbox=args.bbox,
        resolution=args.resolution,
    )


if __name__ == "__main__":
    main()
