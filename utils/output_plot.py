#!/usr/bin/env python3
"""
utils/plot_flood_depth.py

A robust, workflow-friendly plotting utility for LPERFECT flood depth outputs.

What it does (in one go):
- Loads an LPERFECT flood depth NetCDF (flood_depth(time, latitude, longitude)).
- Loads a domain NetCDF containing a DEM (dem(latitude, longitude) or similar).
- Checks grid alignment and orientation (monotonic lat/lon, extents).
- Optionally regrids flood depth to the DEM grid (or vice versa) using xarray interpolation.
- Builds a hillshade from the DEM and overlays flood depth on top.
- Handles FillValue/NaNs cleanly and supports threshold masking.
- Supports percentile-based clipping for better contrast.
- Can generate a single PNG or one PNG per time step (batch/animation-ready).
- Can optionally overlay vector boundaries from GeoJSON/Shapefile (if geopandas is installed).
- Uses logging for traceable runs.

Dependencies:
- Required: matplotlib, numpy, xarray
- Recommended: netCDF4 (backend), scipy (for some interpolation paths)
- Optional overlays: geopandas, shapely
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Optional, Tuple, Any, List

import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
import matplotlib.patheffects as pe

LOG = logging.getLogger("plot_flood_depth")


# -----------------------------
# Small helpers (pedagogical)
# -----------------------------

def _to_float(x: Any) -> Optional[float]:
    """Convert a scalar to float, returning None on failure."""
    try:
        if x is None:
            return None
        return float(x)
    except Exception:
        return None


def _is_monotonic_1d(a: np.ndarray) -> bool:
    """True if 1D array is strictly increasing or strictly decreasing."""
    if a.ndim != 1 or a.size < 2:
        return False
    d = np.diff(a.astype(float))
    return bool(np.all(d > 0) or np.all(d < 0))


def _infer_coord_name(ds: xr.Dataset, candidates: List[str]) -> Optional[str]:
    """Return the first candidate that exists as a coord or variable in ds."""
    for c in candidates:
        if c in ds.coords or c in ds.variables:
            return c
    return None


def _get_fill_value(da: xr.DataArray) -> Optional[float]:
    """Read _FillValue / missing_value from attributes if present."""
    for k in ("_FillValue", "missing_value"):
        if k in da.attrs:
            v = _to_float(da.attrs.get(k))
            if v is not None:
                return v
    return None


def _mask_fill_and_nan(da: xr.DataArray) -> xr.DataArray:
    """Replace FillValue with NaN and ensure float dtype."""
    out = da.astype(float)
    fill = _get_fill_value(out)
    if fill is not None:
        out = out.where(out != fill)
    return out


def _percentile_vmax(data: np.ndarray, p: float) -> float:
    """Compute percentile-based vmax on finite values."""
    flat = data[np.isfinite(data)]
    if flat.size == 0:
        return 0.0
    return float(np.percentile(flat, p))


def _pick_var(ds: xr.Dataset, preferred: str, fallback: List[str]) -> xr.DataArray:
    """Pick preferred variable, else fallback list, else raise."""
    if preferred in ds.variables:
        return ds[preferred]
    for n in fallback:
        if n in ds.variables:
            LOG.warning("Variable '%s' not found, using '%s' instead.", preferred, n)
            return ds[n]
    raise KeyError(f"Could not find '{preferred}' (or fallback) in dataset variables: {list(ds.variables)}")


def _subset_bbox(ds: xr.Dataset, lat_name: str, lon_name: str, bbox: Tuple[float, float, float, float]) -> xr.Dataset:
    """Subset dataset to a lon/lat bounding box; raise if empty."""
    min_lon, min_lat, max_lon, max_lat = bbox
    lat_vals = ds[lat_name].values
    lon_vals = ds[lon_name].values
    lat_slice = slice(min_lat, max_lat) if lat_vals[0] < lat_vals[-1] else slice(max_lat, min_lat)
    lon_slice = slice(min_lon, max_lon) if lon_vals[0] < lon_vals[-1] else slice(max_lon, min_lon)
    ds_sub = ds.sel({lat_name: lat_slice, lon_name: lon_slice})
    if ds_sub[lat_name].size == 0 or ds_sub[lon_name].size == 0:
        raise ValueError(f"BBox {bbox} returned an empty selection for dataset coords ({lat_name}, {lon_name}).")
    return ds_sub


def _pick_label_column(gdf, requested: Optional[str]) -> Optional[str]:
    """Return a column name to use for labels (requested, common defaults, or first string-like)."""
    if requested:
        if requested in gdf.columns:
            return requested
        LOG.warning("Requested label column '%s' not found in overlay; falling back to auto-detection.", requested)

    for cand in ("name", "Name", "NAME", "nome", "NOME", "NOME_COM", "COMUNE"):
        if cand in gdf.columns:
            return cand

    for col in gdf.columns:
        if col == gdf.geometry.name:
            continue
        if gdf[col].dtype == object:
            return col
        if np.issubdtype(gdf[col].dtype, np.number):
            return col
    return None


# -----------------------------
# Regridding
# -----------------------------

def maybe_regrid_to_match(
    flood_da: xr.DataArray,
    dem_da: xr.DataArray,
    flood_lat: str,
    flood_lon: str,
    dem_lat: str,
    dem_lon: str,
    mode: str,
) -> Tuple[xr.DataArray, xr.DataArray]:
    """Optionally interpolate flood or DEM to match the other grid."""
    if mode == "none":
        return flood_da, dem_da

    if mode == "flood_to_dem":
        LOG.info("Regridding flood -> DEM grid (xarray.interp, linear).")
        flood_rg = flood_da.interp({flood_lat: dem_da[dem_lat], flood_lon: dem_da[dem_lon]}, method="linear")
        flood_rg = flood_rg.rename({flood_lat: dem_lat, flood_lon: dem_lon})
        return flood_rg, dem_da

    if mode == "dem_to_flood":
        LOG.info("Regridding DEM -> flood grid (xarray.interp, linear).")
        dem_rg = dem_da.interp({dem_lat: flood_da[flood_lat], dem_lon: flood_da[flood_lon]}, method="linear")
        dem_rg = dem_rg.rename({dem_lat: flood_lat, dem_lon: flood_lon})
        return flood_da, dem_rg

    raise ValueError(f"Unknown regrid mode: {mode}")


# -----------------------------
# Optional vector overlay
# -----------------------------

def load_vectors(path: str):
    """Load GeoJSON/Shapefile using geopandas (optional dependency)."""
    try:
        import geopandas as gpd  # type: ignore
    except Exception as e:
        raise RuntimeError("geopandas not installed. Install: pip install geopandas") from e

    gdf = gpd.read_file(path)
    if gdf.crs is not None:
        gdf = gdf.to_crs("EPSG:4326")
    return gdf


# -----------------------------
# Plotting
# -----------------------------

def plot_one(
    flood_da: xr.DataArray,
    dem_da: xr.DataArray,
    lat_name: str,
    lon_name: str,
    title: str,
    out_png: Optional[str],
    cmap_flood: str,
    flood_alpha: float,
    vmin: float,
    vmax: Optional[float],
    vmax_percentile: Optional[float],
    threshold: Optional[float],
    log_scale: bool,
    hillshade_azdeg: float,
    hillshade_altdeg: float,
    hillshade_vert_exag: float,
    overlay_vectors: Optional[str],
    vector_linewidth: float,
    vector_alpha: float,
    overlay_label_field: Optional[str],
    overlay_label_size: float,
    dpi: int,
) -> None:
    """Render one map: DEM hillshade + flood overlay."""
    flood = _mask_fill_and_nan(flood_da)
    dem = _mask_fill_and_nan(dem_da)

    if threshold is not None:
        flood = flood.where(flood > threshold)

    # Always hide non-positive flood depths (transparent).
    flood = flood.where(flood > 0)

    flood_np = np.asarray(flood.values, dtype=float)

    if vmax is None and vmax_percentile is not None:
        vmax = _percentile_vmax(flood_np, vmax_percentile)

    if vmax is None or not np.isfinite(vmax) or vmax <= 0:
        fallback = max(vmin, 1e-6 if log_scale else 0.0)
        vmax = fallback + (1e-3 if log_scale else 1e-3)
        LOG.warning("Using fallback vmax=%s for plotting (data after masking may be empty).", vmax)

    ls = LightSource(azdeg=hillshade_azdeg, altdeg=hillshade_altdeg)
    # LightSource.shade expects a Colormap, not a string (matplotlib>=3.9 raises TypeError)
    shade = ls.shade(
        np.asarray(dem.values, dtype=float),
        cmap=plt.get_cmap("Greys"),
        blend_mode="overlay",
        vert_exag=hillshade_vert_exag,
    )

    lon = np.asarray(dem[lon_name].values, dtype=float)
    lat = np.asarray(dem[lat_name].values, dtype=float)
    if lon.ndim != 1 or lat.ndim != 1:
        raise ValueError("This script expects 1D latitude/longitude coordinates (rectilinear grid).")

    extent = [float(lon.min()), float(lon.max()), float(lat.min()), float(lat.max())]

    fig, ax = plt.subplots(figsize=(10, 8))
    ax.imshow(shade, extent=extent, origin="lower")

    if log_scale:
        from matplotlib.colors import LogNorm
        norm = LogNorm(vmin=max(vmin, 1e-6), vmax=vmax)
        cs = ax.pcolormesh(flood[lon_name], flood[lat_name], flood, cmap=cmap_flood, shading="auto", alpha=flood_alpha, norm=norm)
        cbar = fig.colorbar(cs, ax=ax, label="Flood depth (m) [log]")
    else:
        cs = ax.pcolormesh(flood[lon_name], flood[lat_name], flood, cmap=cmap_flood, shading="auto", alpha=flood_alpha, vmin=vmin, vmax=vmax)
        cbar = fig.colorbar(cs, ax=ax, label="Flood depth (m)")

    if overlay_vectors:
        try:
            gdf = load_vectors(overlay_vectors)
            gdf.boundary.plot(ax=ax, linewidth=vector_linewidth, alpha=vector_alpha)

            label_col = _pick_label_column(gdf, overlay_label_field)
            if label_col:
                for _, row in gdf.iterrows():
                    geom = row.geometry
                    label = row[label_col]
                    if geom is None or geom.is_empty or label is None or (isinstance(label, (float, np.floating)) and np.isnan(label)):
                        continue
                    centroid = geom.centroid
                    ax.text(
                        centroid.x,
                        centroid.y,
                        str(label),
                        ha="center",
                        va="center",
                        fontsize=overlay_label_size,
                        color="black",
                        alpha=vector_alpha,
                        path_effects=[pe.withStroke(linewidth=1.5, foreground="white")],
                    )
            else:
                LOG.info("Overlay provided but no suitable label column found; skipping labels.")
        except Exception as e:
            LOG.error("Vector overlay failed (%s). Continuing without vectors.", e)

    ax.set_title(title)
    ax.set_xlabel("Longitude")
    ax.set_ylabel("Latitude")
    plt.tight_layout()

    if out_png:
        Path(out_png).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_png, dpi=dpi)
        LOG.info("Saved figure: %s", out_png)
        plt.close(fig)
    else:
        plt.show()


# -----------------------------
# Main (CLI)
# -----------------------------

def main() -> int:
    ap = argparse.ArgumentParser(description="Plot LPERFECT flood depth over DEM hillshade.")
    ap.add_argument("--flood", required=True, help="LPERFECT flood depth NetCDF path.")
    ap.add_argument("--domain", required=True, help="Domain NetCDF path containing DEM.")
    ap.add_argument("--out", default=None, help="Output PNG path. If omitted, show interactive window.")
    ap.add_argument("--out-dir", default=None, help="If set with --all-times, save frames here.")
    ap.add_argument("--all-times", action="store_true", help="Render one PNG per time index.")
    ap.add_argument("--time-index", type=int, default=0, help="Time index (when not using --all-times).")

    ap.add_argument("--flood-var", default="flood_depth", help="Flood depth variable name.")
    ap.add_argument("--dem-var", default="dem", help="DEM variable name.")
    ap.add_argument("--lat-name", default=None, help="Latitude coord name (auto if omitted).")
    ap.add_argument("--lon-name", default=None, help="Longitude coord name (auto if omitted).")

    ap.add_argument("--regrid", choices=["none", "flood_to_dem", "dem_to_flood"], default="flood_to_dem",
                    help="Align grids via interpolation. Default: flood_to_dem.")

    ap.add_argument("--title", default=None, help="Custom plot title.")
    ap.add_argument("--cmap-flood", default="Blues", help="Colormap for flood depth.")
    ap.add_argument("--alpha", type=float, default=0.7, help="Flood overlay alpha (0..1).")
    ap.add_argument("--vmin", type=float, default=0.0, help="Linear min flood depth.")
    ap.add_argument("--vmax", type=float, default=None, help="Linear max flood depth.")
    ap.add_argument("--vmax-percentile", type=float, default=99.5, help="Percentile vmax if --vmax not set.")
    ap.add_argument("--threshold", type=float, default=None, help="Mask flood depth <= threshold (m).")
    ap.add_argument("--log-scale", action="store_true", help="Use log scaling for flood depth colors.")

    ap.add_argument("--azdeg", type=float, default=315.0, help="Hillshade azimuth (deg).")
    ap.add_argument("--altdeg", type=float, default=45.0, help="Hillshade altitude (deg).")
    ap.add_argument("--vert-exag", type=float, default=1.0, help="Hillshade vertical exaggeration.")

    ap.add_argument("--overlay", "--overlay-vector", dest="overlay", default=None,
                    help="GeoJSON/Shapefile path to overlay (optional).")
    ap.add_argument("--vector-linewidth", type=float, default=1.0)
    ap.add_argument("--vector-alpha", type=float, default=0.9)
    ap.add_argument("--overlay-label-field", default=None,
                    help="Column to use as label for overlay features (auto-detect if omitted).")
    ap.add_argument("--overlay-label-size", type=float, default=8.0, help="Font size for overlay labels.")

    ap.add_argument("--bbox", nargs=4, type=float, metavar=("MIN_LON", "MIN_LAT", "MAX_LON", "MAX_LAT"),
                    help="Clip both DEM and flood data to a lon/lat bounding box.")

    ap.add_argument("--dpi", type=int, default=150, help="PNG DPI when saving.")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")

    LOG.info("Opening flood dataset: %s", args.flood)
    flood_ds = xr.open_dataset(args.flood)

    LOG.info("Opening domain dataset: %s", args.domain)
    domain_ds = xr.open_dataset(args.domain)

    flood_var = _pick_var(flood_ds, args.flood_var, ["flood_depth", "water_depth", "depth"])
    dem_var = _pick_var(domain_ds, args.dem_var, ["dem", "elevation", "h", "topography"])

    flood_lat = args.lat_name or _infer_coord_name(flood_ds, ["latitude", "lat", "y"])
    flood_lon = args.lon_name or _infer_coord_name(flood_ds, ["longitude", "lon", "x"])
    dem_lat = args.lat_name or _infer_coord_name(domain_ds, ["latitude", "lat", "y"])
    dem_lon = args.lon_name or _infer_coord_name(domain_ds, ["longitude", "lon", "x"])

    if not flood_lat or not flood_lon:
        raise ValueError("Could not infer flood lat/lon names. Use --lat-name/--lon-name.")
    if not dem_lat or not dem_lon:
        raise ValueError("Could not infer DEM lat/lon names. Use --lat-name/--lon-name.")

    # Sort by latitude if decreasing, so we keep a consistent north-up map
    if flood_ds[flood_lat].values[0] > flood_ds[flood_lat].values[-1]:
        flood_ds = flood_ds.sortby(flood_ds[flood_lat])
    if domain_ds[dem_lat].values[0] > domain_ds[dem_lat].values[-1]:
        domain_ds = domain_ds.sortby(domain_ds[dem_lat])

    if args.bbox:
        bbox = tuple(args.bbox)
        LOG.info("Applying bounding box (min_lon, min_lat, max_lon, max_lat) = %s", bbox)
        flood_ds = _subset_bbox(flood_ds, flood_lat, flood_lon, bbox)
        domain_ds = _subset_bbox(domain_ds, dem_lat, dem_lon, bbox)

    flood_var = flood_ds[flood_var.name]
    dem_var = domain_ds[dem_var.name]

    # Sanity checks for coordinates
    if not _is_monotonic_1d(np.asarray(flood_ds[flood_lat].values)):
        raise ValueError(f"Flood coord '{flood_lat}' is not strictly monotonic 1D.")
    if not _is_monotonic_1d(np.asarray(flood_ds[flood_lon].values)):
        raise ValueError(f"Flood coord '{flood_lon}' is not strictly monotonic 1D.")
    if not _is_monotonic_1d(np.asarray(domain_ds[dem_lat].values)):
        raise ValueError(f"DEM coord '{dem_lat}' is not strictly monotonic 1D.")
    if not _is_monotonic_1d(np.asarray(domain_ds[dem_lon].values)):
        raise ValueError(f"DEM coord '{dem_lon}' is not strictly monotonic 1D.")

    dem_da = dem_var
    if dem_lat not in dem_da.dims or dem_lon not in dem_da.dims:
        raise ValueError(f"DEM dims {dem_da.dims} do not include ({dem_lat}, {dem_lon}).")

    time_dim = "time" if "time" in flood_var.dims else None

    def make_title(ti: int) -> str:
        if args.title:
            return args.title
        return f"LPERFECT flood depth – time index {ti}" if time_dim else "LPERFECT flood depth"

    if args.all_times:
        if not args.out_dir:
            raise ValueError("--all-times requires --out-dir.")
        if not time_dim:
            indices = [0]
            LOG.warning("No 'time' dimension found; plotting one frame.")
        else:
            indices = list(range(int(flood_ds.dims.get(time_dim, 1))))
    else:
        indices = [args.time_index]

    for ti in indices:
        if time_dim:
            if ti < 0 or ti >= int(flood_ds.dims.get(time_dim, 1)):
                raise ValueError(f"time index {ti} out of range.")
            flood_da = flood_var.isel({time_dim: ti})
        else:
            flood_da = flood_var

        if flood_lat not in flood_da.dims or flood_lon not in flood_da.dims:
            raise ValueError(f"Flood dims {flood_da.dims} do not include ({flood_lat}, {flood_lon}).")

        flood_aligned, dem_aligned = maybe_regrid_to_match(
            flood_da=flood_da,
            dem_da=dem_da,
            flood_lat=flood_lat,
            flood_lon=flood_lon,
            dem_lat=dem_lat,
            dem_lon=dem_lon,
            mode=args.regrid,
        )

        if args.all_times:
            out_png = str(Path(args.out_dir) / f"flood_depth_t{ti:03d}.png")
        else:
            out_png = args.out

        # Choose coordinate names for plotting depending on regrid choice
        if args.regrid == "flood_to_dem":
            lat_name, lon_name = dem_lat, dem_lon
        elif args.regrid == "dem_to_flood":
            lat_name, lon_name = flood_lat, flood_lon
        else:
            # no regrid: use DEM coords for hillshade extent; requires both grids comparable
            lat_name, lon_name = dem_lat, dem_lon

        plot_one(
            flood_da=flood_aligned,
            dem_da=dem_aligned,
            lat_name=lat_name,
            lon_name=lon_name,
            title=make_title(ti),
            out_png=out_png,
            cmap_flood=args.cmap_flood,
            flood_alpha=args.alpha,
            vmin=args.vmin,
            vmax=args.vmax,
            vmax_percentile=args.vmax_percentile,
            threshold=args.threshold,
            log_scale=args.log_scale,
            hillshade_azdeg=args.azdeg,
            hillshade_altdeg=args.altdeg,
            hillshade_vert_exag=args.vert_exag,
            overlay_vectors=args.overlay,
            vector_linewidth=args.vector_linewidth,
            vector_alpha=args.vector_alpha,
            overlay_label_field=args.overlay_label_field,
            overlay_label_size=args.overlay_label_size,
            dpi=args.dpi,
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
