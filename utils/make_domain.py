#!/usr/bin/env python3
"""Create a domain NetCDF from DEM/CN/D8/river mask inputs."""
from __future__ import annotations

import argparse  # Parse command-line arguments for the CLI.
from dataclasses import dataclass  # Provide a simple data container for grid metadata.
from datetime import datetime, timezone  # Timestamp domain creation for metadata.
import importlib.util  # Detect optional dependencies at runtime.
from typing import Iterable, Optional  # Define type hints for optional iterables.

import numpy as np  # Numerical arrays and math utilities.
from tqdm import tqdm  # Progress bar for long-running workflows.
import xarray as xr  # Dataset and DataArray abstractions for NetCDF.


@dataclass
class GridSpec:
    longitude_name: str  # Name of the x-coordinate dimension.
    latitude_name: str  # Name of the y-coordinate dimension.
    x: np.ndarray  # 1D coordinate array for x.
    y: np.ndarray  # 1D coordinate array for y.


D8_CODES = np.array([1, 2, 4, 8, 16, 32, 64, 128], dtype=np.int32)  # ESRI D8 codes.
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


def _normalize_da(da: xr.DataArray, longitude_name: str, latitude_name: str) -> xr.DataArray:
    # Ensure the variable is 2D so downstream math behaves as expected.
    if len(da.dims) != 2:
        raise ValueError(f"Expected 2D variable; got dims={da.dims}")
    # Rename dimensions to the expected x/y names if needed.
    if da.dims != (latitude_name, longitude_name):
        da = da.rename({da.dims[0]: latitude_name, da.dims[1]: longitude_name})
    # Return the normalized DataArray to the caller.
    return da


def _slice_for_bbox(coord: np.ndarray, minv: float, maxv: float) -> slice:
    # Detect whether the coordinate array increases or decreases.
    ascending = coord[0] < coord[-1]
    # Build a slice that respects axis direction.
    if ascending:
        return slice(minv, maxv)
    return slice(maxv, minv)


def _build_axis(minv: float, maxv: float, step: float, ascending: bool) -> np.ndarray:
    # Guard against invalid resolution values early.
    if step <= 0:
        raise ValueError("Resolution must be positive")
    step = float(step)  # Normalize step to a float for numpy arange.
    # Flip the sign if the axis should be descending.
    if not ascending:
        step = -step
    start = minv if ascending else maxv  # Pick the correct start bound.
    end = maxv if ascending else minv  # Pick the correct end bound.
    # Include the last cell by extending the end slightly.
    return np.arange(start, end + step * 0.5, step)


def _grid_from_da(da: xr.DataArray, longitude_name: str, latitude_name: str) -> GridSpec:
    # Extract coordinate arrays into numpy for downstream operations.
    x = np.asarray(da[longitude_name].values)
    y = np.asarray(da[latitude_name].values)
    # Return a GridSpec for consistent grid metadata handling.
    return GridSpec(longitude_name=longitude_name, latitude_name=latitude_name, x=x, y=y)


def _apply_bbox(da: xr.DataArray, grid: GridSpec, bbox: Optional[Iterable[float]]) -> xr.DataArray:
    # Short-circuit if no bounding box was requested.
    if bbox is None:
        return da
    xmin, ymin, xmax, ymax = bbox  # Unpack the bounding box values.
    # Select the subset of the domain inside the bounding box.
    return da.sel(
        {
            grid.longitude_name: _slice_for_bbox(grid.x, xmin, xmax),
            grid.latitude_name: _slice_for_bbox(grid.y, ymin, ymax),
        }
    )


def _apply_resolution(da: xr.DataArray, grid: GridSpec, res: Optional[Iterable[float]], method: str) -> xr.DataArray:
    # Skip resampling if no target resolution was provided.
    if res is None:
        return da
    res_vals = list(res)  # Convert iterable to list for length checks.
    # Allow a single value to apply to both x and y axes.
    if len(res_vals) == 1:
        dx = dy = res_vals[0]
    # Allow two values for explicit dx/dy control.
    elif len(res_vals) == 2:
        dx, dy = res_vals
    else:
        raise ValueError("Resolution expects 1 or 2 values")
    # Determine bounds of the original grid.
    xmin, xmax = float(np.min(grid.x)), float(np.max(grid.x))
    ymin, ymax = float(np.min(grid.y)), float(np.max(grid.y))
    x_asc = grid.x[0] < grid.x[-1]  # Check axis orientation for x.
    y_asc = grid.y[0] < grid.y[-1]  # Check axis orientation for y.
    # Build the new coordinate axes with the desired spacing.
    new_x = _build_axis(xmin, xmax, dx, x_asc)
    new_y = _build_axis(ymin, ymax, dy, y_asc)
    # Interpolate the data onto the new grid.
    return _interp_da(da, {grid.longitude_name: new_x, grid.latitude_name: new_y}, method=method)


def _regrid_to_target(da: xr.DataArray, target: GridSpec, method: str) -> xr.DataArray:
    # Interpolate onto the target grid coordinates.
    return _interp_da(da, {target.longitude_name: target.x, target.latitude_name: target.y}, method=method)


def _interp_da(da: xr.DataArray, coords: dict[str, np.ndarray], method: str) -> xr.DataArray:
    # Use xarray's nearest-neighbor selection if requested.
    if method == "nearest":
        return da.sel(coords, method="nearest")
    # For linear interpolation we need scipy available.
    if importlib.util.find_spec("scipy") is None:
        raise ImportError(
            "scipy is required for linear interpolation. "
            "Install scipy or avoid --resolution/linear interpolation."
        )
    # Delegate to xarray's interpolation utility.
    return da.interp(coords, method=method)


def _pick_fill_value(da: xr.DataArray, default: float | int) -> float | int:
    for key in ("_FillValue", "missing_value"):
        if key in da.attrs:
            return da.attrs[key]
    return default


def _coord_attrs(name: str, src_attrs: dict) -> dict:
    attrs = dict(src_attrs)
    lname = name.lower()
    if "lat" in lname:
        attrs.setdefault("description", "Latitude")
        attrs.setdefault("units", "degrees_north")
        attrs.setdefault("long_name", "latitude")
    if "lon" in lname:
        attrs.setdefault("description", "Longitude")
        attrs.setdefault("units", "degrees_east")
        attrs.setdefault("long_name", "longitude")
    return attrs


def _load_var(path: str, var_name: str, longitude_name: str, latitude_name: str) -> tuple[xr.DataArray, xr.Dataset]:
    raster_exts = (".tif", ".tiff", ".geotiff", ".gtiff")  # Supported raster file extensions.
    # Branch based on file type to load raster or NetCDF inputs.
    if path.lower().endswith(raster_exts):
        try:
            import rioxarray as rxr  # Lazy import to keep raster optional.
        except ImportError as exc:
            raise ImportError(
                "Reading GeoTIFF inputs requires rioxarray and rasterio. "
                "Install them or convert the input to NetCDF."
            ) from exc
        da = rxr.open_rasterio(path)  # Load raster data as a DataArray.
        # Ensure a single-band raster is selected for use.
        if "band" in da.dims:
            if da.sizes["band"] < 1:
                raise ValueError(f"No bands found in raster {path}")
            da = da.isel(band=0, drop=True)
        da = da.rename(var_name)  # Rename the DataArray to the desired variable name.
        ds = da.to_dataset()  # Convert to a Dataset for consistent handling.
    else:
        ds = xr.open_dataset(path)  # Load NetCDF (or similar) via xarray.
    # Validate that the requested variable is present.
    if var_name not in ds:
        raise KeyError(f"Variable '{var_name}' not found in {path}")
    # Normalize the variable to a consistent 2D layout.
    da = _normalize_da(ds[var_name], longitude_name=longitude_name, latitude_name=latitude_name)
    return da, ds  # Return both DataArray and its parent Dataset.


def _compute_d8(dem: np.ndarray, dx: float, dy: float) -> np.ndarray:
    dem = dem.astype(np.float64)  # Work in floating-point for slope math.
    nrows, ncols = dem.shape  # Capture the grid shape.
    # Pad edges with NaNs so neighbor lookups stay in bounds.
    pad = np.pad(dem, 1, mode="constant", constant_values=np.nan)
    slopes = []  # Collect slope arrays for each D8 direction.
    dist_diag = float(np.hypot(dx, dy))  # Diagonal distance between cells.
    dist_card = float(np.mean([abs(dx), abs(dy)]))  # Cardinal distance between cells.
    distances = [dist_card, dist_diag, dist_card, dist_diag, dist_card, dist_diag, dist_card, dist_diag]
    # Compute slopes to each neighbor direction.
    for (dr, dc), dist in zip(D8_DIRS, distances):
        neigh = pad[1 + dr:1 + dr + nrows, 1 + dc:1 + dc + ncols]  # Neighbor elevation.
        slope = (dem - neigh) / dist  # Positive slope means flow toward neighbor.
        slope = np.where(np.isfinite(dem) & np.isfinite(neigh), slope, -np.inf)  # Mask invalid cells.
        slopes.append(slope)  # Store slope surface for this direction.
    slope_stack = np.stack(slopes, axis=0)  # Shape (8, nrows, ncols).
    max_idx = np.argmax(slope_stack, axis=0)  # Direction index of steepest descent.
    max_slope = np.max(slope_stack, axis=0)  # Steepest slope value per cell.
    d8 = np.zeros_like(dem, dtype=np.int32)  # Initialize output with zeros.
    valid = max_slope > 0  # Only assign flow directions where descent exists.
    for idx, code in enumerate(D8_CODES):
        mask = valid & (max_idx == idx)  # Select cells matching this direction.
        d8[mask] = code  # Assign the D8 code.
    return d8  # Return the completed D8 grid.


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
    longitude_name: str,
    latitude_name: str,
    bbox: Optional[Iterable[float]],
    resolution: Optional[Iterable[float]],
) -> None:
    # Count progress steps so the bar reflects optional inputs.
    total_steps = 4  # DEM load, bbox, resolution, grid prep.
    total_steps += 1  # D8 step (load or compute).
    total_steps += 1  # CN step (load or default).
    if mask_path:
        total_steps += 1  # Optional mask step.
    total_steps += 3  # Dataset assembly, attributes, write to disk.

    with tqdm(total=total_steps, desc="Building domain", unit="step") as progress:
        progress.set_description("Loading DEM")
        dem_da, dem_ds = _load_var(dem_path, dem_var, longitude_name=longitude_name, latitude_name=latitude_name)
        progress.update(1)

        progress.set_description("Applying bounding box")
        dem_da = _apply_bbox(dem_da, _grid_from_da(dem_da, longitude_name, latitude_name), bbox)
        progress.update(1)

        progress.set_description("Applying resolution")
        dem_da = _apply_resolution(dem_da, _grid_from_da(dem_da, longitude_name, latitude_name), resolution, method="linear")
        progress.update(1)

        progress.set_description("Preparing grid")
        grid = _grid_from_da(dem_da, longitude_name, latitude_name)
        dx = float(np.median(np.abs(np.diff(grid.x))))
        dy = float(np.median(np.abs(np.diff(grid.y))))
        progress.update(1)

        progress.set_description("Preparing D8")
        d8_da = None
        if d8_path:
            d8_da, d8_ds = _load_var(d8_path, d8_var, longitude_name=longitude_name, latitude_name=latitude_name)
            d8_da = _apply_bbox(d8_da, _grid_from_da(d8_da, longitude_name, latitude_name), bbox)
            d8_da = _regrid_to_target(d8_da, grid, method="nearest")
            d8 = np.asarray(d8_da.values).astype(np.int32)
            d8_ds.close()
        else:
            d8 = _compute_d8(np.asarray(dem_da.values), dx=dx, dy=dy)
        progress.update(1)

        progress.set_description("Preparing CN")
        cn_da = None
        if cn_path:
            cn_da, cn_ds = _load_var(cn_path, cn_var, longitude_name=longitude_name, latitude_name=latitude_name)
            cn_da = _apply_bbox(cn_da, _grid_from_da(cn_da, longitude_name, latitude_name), bbox)
            cn_da = _regrid_to_target(cn_da, grid, method="nearest")
            cn = np.asarray(cn_da.values).astype(np.float64)
            cn_ds.close()
        else:
            cn = np.zeros_like(dem_da.values, dtype=np.float64)
        progress.update(1)

        channel_mask = None  # Default to no channel mask if none provided.
        mask_da = None
        if mask_path:
            progress.set_description("Preparing channel mask")
            mask_da, mask_ds = _load_var(mask_path, mask_var, longitude_name=longitude_name, latitude_name=latitude_name)
            mask_da = _apply_bbox(mask_da, _grid_from_da(mask_da, longitude_name, latitude_name), bbox)
            mask_da = _regrid_to_target(mask_da, grid, method="nearest")
            channel_mask = (np.asarray(mask_da.values) > 0).astype(np.int8)
            mask_ds.close()
            progress.update(1)

        progress.set_description("Assembling dataset")
        x_coord = xr.DataArray(
            grid.x,
            dims=(longitude_name,),
            attrs=_coord_attrs(longitude_name, dem_da[longitude_name].attrs),
        )
        y_coord = xr.DataArray(
            grid.y,
            dims=(latitude_name,),
            attrs=_coord_attrs(latitude_name, dem_da[latitude_name].attrs),
        )
        fill_dem = _pick_fill_value(dem_da, np.nan)
        fill_d8 = None
        fill_cn = None
        fill_mask = None
        ds_out = xr.Dataset()  # Start with an empty dataset.
        ds_out = ds_out.assign_coords(
            {
                longitude_name: x_coord,
                latitude_name: y_coord,
            }
        )
        fill_d8 = _pick_fill_value(d8_da if d8_path else dem_da, -9999)
        fill_cn = _pick_fill_value(cn_da if cn_path else dem_da, -9999.0)
        if mask_da is not None:
            fill_mask = _pick_fill_value(mask_da, 0)
        progress.update(1)

        progress.set_description("Writing variables")
        dem_attrs = dict(dem_da.attrs)  # Copy DEM attributes for reuse.
        ds_out[dem_var] = xr.DataArray(
            np.asarray(dem_da.values, dtype=np.float64),
            dims=(latitude_name, longitude_name),
            attrs={
                "standard_name": dem_attrs.get("standard_name", "surface_altitude"),
                "long_name": dem_attrs.get("long_name", "digital_elevation_model"),
                "units": dem_attrs.get("units", "m"),
                "_FillValue": fill_dem,
            },
        )

        ds_out[d8_var] = xr.DataArray(
            d8.astype(np.int32),
            dims=(latitude_name, longitude_name),
            attrs={
                "long_name": "D8_flow_direction",
                "flag_values": D8_CODES.tolist(),
                "flag_meanings": "E SE S SW W NW N NE",
                "comment": "ESRI D8 encoding (see LPERFECT model.encoding).",
                "_FillValue": int(fill_d8),
            },
        )

        ds_out[cn_var] = xr.DataArray(
            cn.astype(np.float64),
            dims=(latitude_name, longitude_name),
            attrs={
                "long_name": "SCS_curve_number",
                "units": "1",
                "valid_min": 0.0,
                "valid_max": 100.0,
                "_FillValue": float(fill_cn),
            },
        )

        if channel_mask is not None:
            ds_out[mask_var] = xr.DataArray(
                channel_mask,
                dims=(latitude_name, longitude_name),
                attrs={
                    "long_name": "channel_mask",
                    "units": "1",
                    "flag_values": "0 1",
                    "flag_meanings": "no_channel channel",
                    "_FillValue": int(fill_mask) if fill_mask is not None else 0,
                },
            )
        progress.update(1)

        progress.set_description("Finalizing attributes")
        grid_mapping = dem_attrs.get("grid_mapping")
        if grid_mapping and grid_mapping in dem_ds:
            ds_out[grid_mapping] = xr.DataArray(0, attrs=dict(dem_ds[grid_mapping].attrs))
            ds_out[dem_var].attrs["grid_mapping"] = grid_mapping
            ds_out[d8_var].attrs["grid_mapping"] = grid_mapping
            ds_out[cn_var].attrs["grid_mapping"] = grid_mapping
            if channel_mask is not None:
                ds_out[mask_var].attrs["grid_mapping"] = grid_mapping

        ds_out.attrs.update(
            {
                "title": "LPERFECT simulation domain",
                "source": "Prepared for LPERFECT (DEM + D8 + CN)",
                "institution": dem_ds.attrs.get("institution", "LPERFECT"),
                "history": f"{datetime.now(timezone.utc).isoformat()}: domain created",
                "references": "Hi-WeFAI / LPERFECT",
                "Conventions": "CF-1.10",
            }
        )
        progress.update(1)

        progress.set_description("Writing NetCDF")
        ds_out.to_netcdf(output_path)
        dem_ds.close()
        progress.update(1)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)  # CLI parser with module docstring.
    parser.add_argument("--dem", required=True, help="Path to DEM NetCDF")
    parser.add_argument("--cn", help="Path to CN NetCDF (optional)")
    parser.add_argument("--d8", help="Path to D8 NetCDF (optional)")
    parser.add_argument("--mask", help="Path to channel/river mask NetCDF (optional)")
    parser.add_argument("--dem-var", default="dem", help="DEM variable name")
    parser.add_argument("--cn-var", default="cn", help="CN variable name")
    parser.add_argument("--d8-var", default="d8", help="D8 variable name")
    parser.add_argument("--mask-var", default="channel_mask", help="Mask variable name")
    parser.add_argument("--longitude-name", default="longitude", help="Longitude coordinate name")
    parser.add_argument("--latitude-name", default="latitude", help="Latitude coordinate name")
    parser.add_argument("--bbox", nargs=4, type=float, metavar=("min_lon", "min_lat", "max_lon", "max_lat"))
    parser.add_argument(
        "--resolution",
        nargs="+",
        type=float,
        metavar=("DX", "DY"),
        help="Target resolution (dx [dy]) in degrees or meters",
    )
    parser.add_argument("--output", required=True, help="Output NetCDF path")
    args = parser.parse_args()  # Parse CLI arguments.

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
        longitude_name=args.longitude_name,
        latitude_name=args.latitude_name,
        bbox=args.bbox,
        resolution=args.resolution,
    )


if __name__ == "__main__":
    main()  # Entry point for CLI execution.
