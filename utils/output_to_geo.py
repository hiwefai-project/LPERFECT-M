#!/usr/bin/env python3
"""
utils/output_to_geo.py

LPERFECT NetCDF (cdl/output.cdl) -> GeoJSON enrichment:
- flood_depth_mean, flood_depth_max
- risk_index_mean, risk_index_max
- flood_depth_pct_gt_thr  (% area where flood_depth > threshold)
- risk_index_class        (R1..R4, PAI/DPC-style classes)

Italian scheme used for classes:
R1 basso, R2 medio/moderato, R3 elevato, R4 molto elevato. :contentReference[oaicite:1]{index=1}

Dependencies:
  pip install xarray netCDF4 numpy pyproj shapely
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, Tuple, Optional, List

import numpy as np
import xarray as xr
from tqdm import tqdm
from pyproj import Transformer
from shapely.geometry import shape, Polygon, MultiPolygon, Point, LineString, MultiLineString
from shapely.ops import transform as shp_transform
from shapely.prepared import prep

# Optional fast zonal stats (rasterization). If unavailable, we fall back to exact (slow) geometry-cell intersections.
try:
    import rasterio
    from rasterio import features as rio_features
    from rasterio.transform import Affine
except Exception:  # pragma: no cover
    rasterio = None
    rio_features = None
    Affine = None

LOG = logging.getLogger("output_to_geo")

LAT_NAME = "latitude"
LON_NAME = "longitude"
TIME_NAME = "time"
FLOOD_NAME = "flood_depth"
RISK_NAME = "risk_index"


def _try_init_mpi(disable_mpi: bool):
    """Return (comm, rank, world_size) or (None, 0, 1) when MPI unavailable/disabled."""
    if disable_mpi:
        return None, 0, 1

    try:
        from mpi4py import MPI
    except Exception:
        return None, 0, 1

    comm = MPI.COMM_WORLD
    world_size = max(1, comm.Get_size())
    rank = comm.Get_rank()
    if world_size <= 1:
        return None, rank, 1
    return comm, rank, world_size


def _calc_thread_workers(requested: int) -> int:
    if requested is None or requested == 0:
        try:
            return max(1, (os.cpu_count() or 1) // 2)
        except Exception:
            return 1
    return max(1, int(requested))


def _is_finite(x: float) -> bool:
    return isinstance(x, (int, float)) and math.isfinite(float(x))


def _as_float(v: Any) -> float:
    try:
        return float(np.asarray(v).item())
    except Exception:
        return float("nan")


def _ensure_1d(da: xr.DataArray, name: str) -> np.ndarray:
    arr = np.asarray(da.values)
    if arr.ndim != 1:
        raise ValueError(f"Expected 1D {name}, got shape={arr.shape}")
    return arr


def nearest_ij(lat1d: np.ndarray, lon1d: np.ndarray, lat: float, lon: float) -> Tuple[int, int]:
    ilat = int(np.argmin(np.abs(lat1d - lat)))
    ilon = int(np.argmin(np.abs(lon1d - lon)))
    return ilat, ilon


def grid_edges_from_centers(vals: np.ndarray) -> np.ndarray:
    vals = np.asarray(vals, dtype=float)
    if vals.size < 2:
        raise ValueError("Need at least 2 points to infer grid edges.")
    d = np.diff(vals)
    if not (np.all(d > 0) or np.all(d < 0)):
        raise ValueError("Latitude/longitude must be strictly monotonic for edge inference.")
    mids = vals[:-1] + 0.5 * d
    first = vals[0] - 0.5 * d[0]
    last = vals[-1] + 0.5 * d[-1]
    return np.concatenate([[first], mids, [last]])


def cell_polygon_lonlat(lon_edges: np.ndarray, lat_edges: np.ndarray, ilon: int, ilat: int) -> Polygon:
    x0 = float(lon_edges[ilon])
    x1 = float(lon_edges[ilon + 1])
    y0 = float(lat_edges[ilat])
    y1 = float(lat_edges[ilat + 1])
    minx, maxx = (x0, x1) if x0 <= x1 else (x1, x0)
    miny, maxy = (y0, y1) if y0 <= y1 else (y1, y0)
    return Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])


def bbox_to_index_range(edges: np.ndarray, vmin: float, vmax: float) -> Tuple[int, int]:
    lo, hi = (vmin, vmax) if vmin <= vmax else (vmax, vmin)
    domain_min = min(float(edges[0]), float(edges[-1]))
    domain_max = max(float(edges[0]), float(edges[-1]))
    lo = max(lo, domain_min)
    hi = min(hi, domain_max)
    if lo > hi:
        return 0, -1

    es = np.sort(edges)
    left_edge_idx = int(np.searchsorted(es, lo, side="right") - 1)
    right_edge_idx = int(np.searchsorted(es, hi, side="left"))
    i0 = max(0, min(len(edges) - 2, left_edge_idx))
    i1 = max(0, min(len(edges) - 2, right_edge_idx - 1))
    return (i0, i1) if i0 <= i1 else (0, -1)


def _is_fill(val: float, fill: Any) -> bool:
    if fill is None:
        return False
    try:
        f = float(fill)
    except Exception:
        return False
    return _is_finite(val) and _is_finite(f) and val == f


def classify_risk_r1_r4(
    risk_value: Optional[float],
    thresholds: Optional[Tuple[float, float, float]] = None,
) -> Optional[str]:
    """
    Compute Italian-style risk class label: R1..R4.

    - If thresholds are provided (t1,t2,t3), we classify:
        <=t1 -> R1, <=t2 -> R2, <=t3 -> R3, else -> R4
    - Else, we try to interpret the model output robustly:
        * if it looks like 1..4 already -> round and map
        * else if it looks like 0..1 -> use default quartiles (0.25,0.5,0.75)
    """
    if risk_value is None or (not _is_finite(risk_value)):
        return None

    rv = float(risk_value)

    if thresholds is not None:
        t1, t2, t3 = thresholds
        if rv <= t1:
            return "R1"
        if rv <= t2:
            return "R2"
        if rv <= t3:
            return "R3"
        return "R4"

    # Case A: already on ~[1..4]
    if 0.5 <= rv <= 4.5:
        k = int(round(rv))
        if 1 <= k <= 4:
            return f"R{k}"

    # Case B: normalized [0..1]
    if 0.0 <= rv <= 1.0:
        if rv <= 0.25:
            return "R1"
        if rv <= 0.50:
            return "R2"
        if rv <= 0.75:
            return "R3"
        return "R4"

    return None


def _load_geojson_feature_collection(path: str) -> Dict[str, Any]:
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(f"Input GeoJSON not found: {path}")

    try:
        raw_text = p.read_text(encoding="utf-8")
    except OSError as exc:
        raise OSError(f"Failed to read GeoJSON '{path}': {exc}") from exc

    if not raw_text.strip():
        raise ValueError(f"Input GeoJSON '{path}' is empty.")

    try:
        gj = json.loads(raw_text)
    except json.JSONDecodeError as exc:
        raise ValueError(
            f"Invalid GeoJSON in '{path}': {exc.msg} (line {exc.lineno} column {exc.colno})"
        ) from exc

    if gj.get("type") != "FeatureCollection":
        raise ValueError("Input GeoJSON must be a FeatureCollection.")

    feats = gj.get("features")
    if feats is None:
        raise ValueError("Input GeoJSON FeatureCollection must contain a 'features' array.")
    if not isinstance(feats, list):
        raise ValueError("GeoJSON 'features' must be a list.")

    return gj



def _grid_is_regular(edges: np.ndarray, rtol: float = 1e-6, atol: float = 0.0) -> bool:
    """Return True if edge spacing is (approximately) constant."""
    e = np.asarray(edges, dtype=float)
    if e.size < 3:
        return True
    d = np.diff(e)
    dm = float(np.mean(np.abs(d)))
    if dm == 0.0:
        return False
    return bool(np.allclose(np.abs(d), dm, rtol=rtol, atol=atol))


def _build_affine_from_edges(lon_edges: np.ndarray, lat_edges: np.ndarray) -> Tuple["Affine", bool]:
    """
    Build a rasterio Affine transform from 1D lon/lat edges.

    Returns (transform, flip_lat):
      - transform maps pixel coordinates (col,row) to lon/lat.
      - flip_lat indicates whether data must be flipped along latitude to match
        rasterio row order (north->south).
    """
    if Affine is None:
        raise RuntimeError("rasterio is not available")

    lon_e = np.asarray(lon_edges, dtype=float)
    lat_e = np.asarray(lat_edges, dtype=float)

    # Use constant resolution assumption (regular grid).
    if not _grid_is_regular(lon_e) or not _grid_is_regular(lat_e):
        raise ValueError("Non-regular lon/lat grid: cannot rasterize safely.")

    xres = float(np.mean(np.abs(np.diff(lon_e))))
    yres = float(np.mean(np.abs(np.diff(lat_e))))

    west = float(min(lon_e[0], lon_e[-1]))
    north = float(max(lat_e[0], lat_e[-1]))

    # rasterio expects row=0 at north, col=0 at west.
    transform = Affine.translation(west, north) * Affine.scale(xres, -yres)

    # If latitude edges are ascending (south->north), the data's first row is south.
    flip_lat = True if lat_e[1] > lat_e[0] else False
    return transform, flip_lat


def compute_area_stats_fast(
    geom_4326,
    *,
    transform: "Affine",
    out_shape: Tuple[int, int],
    flood_np: np.ndarray,
    risk_np: np.ndarray,
    fill_flood: Any,
    fill_risk: Any,
    lat_centers_for_rows: np.ndarray,
    depth_thr_m: Optional[float],
    all_touched: bool = False,
    area_weighting: str = "coslat",
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Very fast zonal stats via rasterization on the model grid.

    - Assumes the model grid is a regular lon/lat grid defined by lon/lat edges.
    - Uses a binary mask created by rasterio.features.rasterize.
    - Mean is area-weighted using a fast approximation:
        * area_weighting="coslat": weights ~ cos(latitude) (good for lon/lat grids)
        * area_weighting="none": unweighted mean over selected cells
    """
    if rio_features is None:
        raise RuntimeError("rasterio is not available")
    if geom_4326.is_empty:
        return None, None, None, None, None

    mask = rio_features.rasterize(
        [(geom_4326, 1)],
        out_shape=out_shape,
        transform=transform,
        fill=0,
        dtype="uint8",
        all_touched=bool(all_touched),
    ).astype(bool)

    if not np.any(mask):
        return None, None, None, None, None

    flood = flood_np
    risk = risk_np

    f_valid = np.isfinite(flood)
    r_valid = np.isfinite(risk)
    if fill_flood is not None and _is_finite(fill_flood):
        f_valid &= (flood != float(fill_flood))
    if fill_risk is not None and _is_finite(fill_risk):
        r_valid &= (risk != float(fill_risk))

    both_valid = mask & f_valid & r_valid

    fmean = rmean = fmax = rmax = None
    if np.any(both_valid):
        if area_weighting == "coslat":
            w_row = np.cos(np.deg2rad(lat_centers_for_rows)).astype(np.float64)
            w = w_row[:, None]
            w_sel = w[both_valid]
            ws = float(np.sum(w_sel))
            if ws > 0.0:
                fmean = float(np.sum(flood[both_valid].astype(np.float64) * w_sel) / ws)
                rmean = float(np.sum(risk[both_valid].astype(np.float64) * w_sel) / ws)
            else:
                fmean = float(np.nanmean(flood[both_valid]))
                rmean = float(np.nanmean(risk[both_valid]))
        else:
            fmean = float(np.nanmean(flood[both_valid]))
            rmean = float(np.nanmean(risk[both_valid]))

        fmax = float(np.nanmax(flood[both_valid]))
        rmax = float(np.nanmax(risk[both_valid]))

    fpct = None
    if depth_thr_m is not None:
        d_valid = mask & f_valid
        if np.any(d_valid):
            if area_weighting == "coslat":
                w_row = np.cos(np.deg2rad(lat_centers_for_rows)).astype(np.float64)
                w = w_row[:, None]
                denom = float(np.sum(w[d_valid]))
                if denom > 0.0:
                    numer = float(np.sum(w[d_valid & (flood > float(depth_thr_m))]))
                    fpct = 100.0 * (numer / denom)
            else:
                denom = int(np.count_nonzero(d_valid))
                if denom > 0:
                    numer = int(np.count_nonzero(d_valid & (flood > float(depth_thr_m))))
                    fpct = 100.0 * (numer / float(denom))

    return fmean, fmax, rmean, rmax, fpct


def compute_area_stats_exact(
    geom_4326,
    lat_edges: np.ndarray,
    lon_edges: np.ndarray,
    flood2d: xr.DataArray,
    risk2d: xr.DataArray,
    fill_flood: Any,
    fill_risk: Any,
    to_area_from4326: Transformer,
    depth_thr_m: Optional[float],
) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float], Optional[float]]:
    """
    Exact (but very slow) area-weighted zonal stats by intersecting each overlapped grid cell.
    Kept for validation/debugging or when rasterization is impossible.
    """
    geom_area = shp_transform(lambda x, y, z=None: to_area_from4326.transform(x, y), geom_4326)
    if geom_area.is_empty:
        return None, None, None, None, None

    geom_area_p = prep(geom_area)

    minx, miny, maxx, maxy = geom_4326.bounds
    ilon0, ilon1 = bbox_to_index_range(lon_edges, minx, maxx)
    ilat0, ilat1 = bbox_to_index_range(lat_edges, miny, maxy)
    if ilon1 < ilon0 or ilat1 < ilat0:
        return None, None, None, None, None

    sum_area = 0.0
    flood_sum = 0.0
    risk_sum = 0.0
    flood_max: Optional[float] = None
    risk_max: Optional[float] = None

    valid_area_depth = 0.0
    gt_area_depth = 0.0

    for ilat_i in range(ilat0, ilat1 + 1):
        for ilon_i in range(ilon0, ilon1 + 1):
            cell_ll = cell_polygon_lonlat(lon_edges, lat_edges, ilon_i, ilat_i)
            cell_area_poly = shp_transform(lambda x, y, z=None: to_area_from4326.transform(x, y), cell_ll)

            if cell_area_poly.is_empty or not geom_area_p.intersects(cell_area_poly):
                continue

            inter = geom_area.intersection(cell_area_poly)
            if inter.is_empty:
                continue

            a = float(inter.area)
            if a <= 0:
                continue

            d = _as_float(flood2d.isel({LAT_NAME: ilat_i, LON_NAME: ilon_i}).values)
            r = _as_float(risk2d.isel({LAT_NAME: ilat_i, LON_NAME: ilon_i}).values)

            d_valid = _is_finite(d) and not _is_fill(d, fill_flood)
            r_valid = _is_finite(r) and not _is_fill(r, fill_risk)

            if d_valid and r_valid:
                sum_area += a
                flood_sum += d * a
                risk_sum += r * a
                flood_max = d if flood_max is None else max(flood_max, d)
                risk_max = r if risk_max is None else max(risk_max, r)

            if depth_thr_m is not None and d_valid:
                valid_area_depth += a
                if d > depth_thr_m:
                    gt_area_depth += a

    if sum_area <= 0:
        flood_mean = None
        risk_mean = None
    else:
        flood_mean = flood_sum / sum_area
        risk_mean = risk_sum / sum_area

    flood_pct = None
    if depth_thr_m is not None and valid_area_depth > 0:
        flood_pct = 100.0 * (gt_area_depth / valid_area_depth)

    return flood_mean, flood_max, risk_mean, risk_max, flood_pct



def _process_feature(
    idx: int,
    feat: Dict[str, Any],
    args: argparse.Namespace,
    *,
    lat: np.ndarray,
    lon: np.ndarray,
    lat_edges: np.ndarray,
    lon_edges: np.ndarray,
    flood2d: xr.DataArray,
    risk2d: xr.DataArray,
    fill_flood: Any,
    fill_risk: Any,
    to4326: Transformer,
    to_area: Transformer,
    back_to_4326: Transformer,
    risk_thr_tuple: Optional[Tuple[float, float, float]],
    # Fast zonal stats context (optional)
    stats_mode: str,
    raster_transform: Optional["Affine"],
    raster_out_shape: Optional[Tuple[int, int]],
    flood_np: Optional[np.ndarray],
    risk_np: Optional[np.ndarray],
    lat_rows: Optional[np.ndarray],
) -> Dict[str, Any]:
    """
    Process a single feature and return a small dict with results.
    Dict keys: idx, feature, updated, skipped.
    """
    feature = feat
    g = feature.get("geometry")
    if not g:
        return {"idx": idx, "feature": feature, "updated": False, "skipped": True}

    geom_src = shape(g)
    if geom_src.is_empty:
        return {"idx": idx, "feature": feature, "updated": False, "skipped": True}

    props = feature.get("properties")
    if props is None or not isinstance(props, dict):
        props = {}
        feature["properties"] = props

    if args.geojson_epsg != 4326:
        geom_4326 = shp_transform(lambda x, y, z=None: to4326.transform(x, y), geom_src)
    else:
        geom_4326 = geom_src

    if geom_4326.is_empty:
        return {"idx": idx, "feature": feature, "updated": False, "skipped": True}

    # ----- POINT: nearest sample -----
    if isinstance(geom_4326, Point) and args.point_buffer_m <= 0:
        lonp, latp = float(geom_4326.x), float(geom_4326.y)
        ilat, ilon = nearest_ij(lat, lon, lat=latp, lon=lonp)

        d = _as_float(flood2d.isel({LAT_NAME: ilat, LON_NAME: ilon}).values)
        r = _as_float(risk2d.isel({LAT_NAME: ilat, LON_NAME: ilon}).values)

        if _is_fill(d, fill_flood):
            d = float("nan")
        if _is_fill(r, fill_risk):
            r = float("nan")

        d_out = None if not _is_finite(d) else d
        r_out = None if not _is_finite(r) else r

        props[args.prop_flood_mean] = d_out
        props[args.prop_flood_max] = d_out
        props[args.prop_risk_mean] = r_out
        props[args.prop_risk_max] = r_out

        if args.depth_threshold_m is not None:
            props[args.prop_flood_pct] = None

        props[args.prop_risk_class] = classify_risk_r1_r4(r_out, risk_thr_tuple)
        props[args.prop_mode] = "point"
        if args.add_grid_idx:
            props["_lperfect_ilat"] = ilat
            props["_lperfect_ilon"] = ilon

        return {"idx": idx, "feature": feature, "updated": True, "skipped": False}

    # ----- POINT buffered to area -----
    if isinstance(geom_4326, Point) and args.point_buffer_m > 0:
        geom_area = shp_transform(lambda x, y, z=None: to_area.transform(x, y), geom_4326).buffer(args.point_buffer_m)
        geom_4326 = shp_transform(lambda x, y, z=None: back_to_4326.transform(x, y), geom_area)

    # ----- LINE buffered to corridor -----
    if isinstance(geom_4326, (LineString, MultiLineString)):
        if args.line_buffer_m > 0:
            geom_area = shp_transform(lambda x, y, z=None: to_area.transform(x, y), geom_4326).buffer(args.line_buffer_m)
            geom_4326 = shp_transform(lambda x, y, z=None: back_to_4326.transform(x, y), geom_area)
        else:
            c = geom_4326.centroid
            lonp, latp = float(c.x), float(c.y)
            ilat, ilon = nearest_ij(lat, lon, lat=latp, lon=lonp)
            d = _as_float(flood2d.isel({LAT_NAME: ilat, LON_NAME: ilon}).values)
            r = _as_float(risk2d.isel({LAT_NAME: ilat, LON_NAME: ilon}).values)
            if _is_fill(d, fill_flood):
                d = float("nan")
            if _is_fill(r, fill_risk):
                r = float("nan")

            d_out = None if not _is_finite(d) else d
            r_out = None if not _is_finite(r) else r

            props[args.prop_flood_mean] = d_out
            props[args.prop_flood_max] = d_out
            props[args.prop_risk_mean] = r_out
            props[args.prop_risk_max] = r_out
            if args.depth_threshold_m is not None:
                props[args.prop_flood_pct] = None

            props[args.prop_risk_class] = classify_risk_r1_r4(r_out, risk_thr_tuple)
            props[args.prop_mode] = "line_centroid"
            return {"idx": idx, "feature": feature, "updated": True, "skipped": False}

    # ----- AREA: Polygon/MultiPolygon (and buffered point/line) -----
    if isinstance(geom_4326, (Polygon, MultiPolygon)):
        fmean = fmax = rmean = rmax = fpct = None
        use_fast = (stats_mode in ("auto", "fast")) and (raster_transform is not None) and (raster_out_shape is not None) and (flood_np is not None) and (risk_np is not None) and (lat_rows is not None)
        if use_fast:
            try:
                fmean, fmax, rmean, rmax, fpct = compute_area_stats_fast(
                    geom_4326,
                    transform=raster_transform,
                    out_shape=raster_out_shape,
                    flood_np=flood_np,
                    risk_np=risk_np,
                    fill_flood=fill_flood,
                    fill_risk=fill_risk,
                    lat_centers_for_rows=lat_rows,
                    depth_thr_m=args.depth_threshold_m,
                    all_touched=args.all_touched,
                    area_weighting=args.area_weighting,
                )
            except Exception as exc:
                # Fall back to exact if fast mode fails (e.g., non-regular grid).
                LOG.debug("Fast zonal stats failed for feature %s: %s; falling back to exact.", idx, exc)
                if stats_mode == "fast":
                    raise
                fmean, fmax, rmean, rmax, fpct = compute_area_stats_exact(
                    geom_4326=geom_4326,
                    lat_edges=lat_edges,
                    lon_edges=lon_edges,
                    flood2d=flood2d,
                    risk2d=risk2d,
                    fill_flood=fill_flood,
                    fill_risk=fill_risk,
                    to_area_from4326=to_area,
                    depth_thr_m=args.depth_threshold_m,
                )
        else:
            fmean, fmax, rmean, rmax, fpct = compute_area_stats_exact(
                geom_4326=geom_4326,
                lat_edges=lat_edges,
                lon_edges=lon_edges,
                flood2d=flood2d,
                risk2d=risk2d,
                fill_flood=fill_flood,
                fill_risk=fill_risk,
                to_area_from4326=to_area,
                depth_thr_m=args.depth_threshold_m,
            )

        props[args.prop_flood_mean] = fmean
        props[args.prop_flood_max] = fmax
        props[args.prop_risk_mean] = rmean
        props[args.prop_risk_max] = rmax

        if args.depth_threshold_m is not None:
            props[args.prop_flood_pct] = fpct
            props.setdefault("_depth_threshold_m", args.depth_threshold_m)

        base_for_class = rmax if args.risk_class_mode == "max" else rmean
        props[args.prop_risk_class] = classify_risk_r1_r4(base_for_class, risk_thr_tuple)

        props[args.prop_mode] = "area"
        return {"idx": idx, "feature": feature, "updated": True, "skipped": False}

    # Other geometries -> centroid fallback
    c = geom_4326.centroid
    lonp, latp = float(c.x), float(c.y)
    ilat, ilon = nearest_ij(lat, lon, lat=latp, lon=lonp)
    d = _as_float(flood2d.isel({LAT_NAME: ilat, LON_NAME: ilon}).values)
    r = _as_float(risk2d.isel({LAT_NAME: ilat, LON_NAME: ilon}).values)
    if _is_fill(d, fill_flood):
        d = float("nan")
    if _is_fill(r, fill_risk):
        r = float("nan")

    d_out = None if not _is_finite(d) else d
    r_out = None if not _is_finite(r) else r

    props[args.prop_flood_mean] = d_out
    props[args.prop_flood_max] = d_out
    props[args.prop_risk_mean] = r_out
    props[args.prop_risk_max] = r_out
    if args.depth_threshold_m is not None:
        props[args.prop_flood_pct] = None

    props[args.prop_risk_class] = classify_risk_r1_r4(r_out, risk_thr_tuple)
    props[args.prop_mode] = "centroid"
    return {"idx": idx, "feature": feature, "updated": True, "skipped": False}


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Enrich GeoJSON with LPERFECT flood depth/risk (mean+max + % area above depth threshold + risk class)."
    )
    ap.add_argument("--nc", required=True)
    ap.add_argument("--geojson-in", required=True)
    ap.add_argument("--geojson-out", required=True)
    ap.add_argument("--time-index", type=int, default=0)

    ap.add_argument("--geojson-epsg", type=int, default=4326)
    ap.add_argument("--area-epsg", type=int, default=6933)

    ap.add_argument("--line-buffer-m", type=float, default=0.0)
    ap.add_argument("--point-buffer-m", type=float, default=0.0)

    ap.add_argument("--depth-threshold-m", type=float, default=None)

    # ---- zonal statistics performance options ----
    ap.add_argument(
        "--stats-mode",
        choices=["auto", "fast", "exact"],
        default="auto",
        help="Zonal stats mode. 'fast' uses rasterization on the regular lon/lat grid (requires rasterio); "
             "'exact' uses polygon-cell intersections (very slow); 'auto' tries fast then falls back to exact.",
    )
    ap.add_argument(
        "--all-touched",
        action="store_true",
        help="When using fast rasterization, include all cells touched by the geometry boundary (faster but less conservative).",
    )
    ap.add_argument(
        "--area-weighting",
        choices=["coslat", "none"],
        default="coslat",
        help="Weighting for fast means on lon/lat grids. 'coslat' approximates true cell areas; 'none' is unweighted.",
    )


    # ---- risk class options ----
    ap.add_argument("--prop-risk-class", default="risk_index_class",
                    help="Property name for classified risk (R1..R4).")
    ap.add_argument("--risk-class-mode", choices=["max", "mean"], default="max",
                    help="Which aggregated risk value to classify (default: max, conservative).")
    ap.add_argument(
        "--risk-class-thresholds",
        default=None,
        help="Optional comma-separated thresholds t1,t2,t3 to map risk->R1..R4 (<=t1 R1, <=t2 R2, <=t3 R3, else R4).",
    )

    ap.add_argument("--prop-flood-mean", default="flood_depth_mean")
    ap.add_argument("--prop-flood-max", default="flood_depth_max")
    ap.add_argument("--prop-risk-mean", default="risk_index_mean")
    ap.add_argument("--prop-risk-max", default="risk_index_max")
    ap.add_argument("--prop-flood-pct", default="flood_depth_pct_gt_thr")
    ap.add_argument("--prop-mode", default="risk_mode")

    ap.add_argument("--add-grid-idx", action="store_true")
    ap.add_argument("--threads", type=int, default=0, help="Number of threads per rank (0 = auto).")
    ap.add_argument("--disable-mpi", action="store_true", help="Disable MPI even if mpi4py is available.")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    args = ap.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s: %(message)s")

    # Parse thresholds if provided
    risk_thr_tuple: Optional[Tuple[float, float, float]] = None
    if args.risk_class_thresholds:
        parts = [p.strip() for p in args.risk_class_thresholds.split(",") if p.strip()]
        if len(parts) != 3:
            raise ValueError("--risk-class-thresholds must have exactly 3 values: t1,t2,t3")
        risk_thr_tuple = (float(parts[0]), float(parts[1]), float(parts[2]))

    to4326 = Transformer.from_crs(f"EPSG:{args.geojson_epsg}", "EPSG:4326", always_xy=True)
    to_area = Transformer.from_crs("EPSG:4326", f"EPSG:{args.area_epsg}", always_xy=True)
    back_to_4326 = Transformer.from_crs(f"EPSG:{args.area_epsg}", "EPSG:4326", always_xy=True)

    comm, rank, world_size = _try_init_mpi(args.disable_mpi)
    threads_per_rank = _calc_thread_workers(args.threads)
    if rank == 0:
        LOG.info("Parallel config: ranks=%d threads-per-rank=%d (MPI %s)", world_size, threads_per_rank,
                 "enabled" if comm else "disabled")

    LOG.info("Opening NetCDF: %s", args.nc)
    ds = xr.open_dataset(args.nc)
    for v in (LAT_NAME, LON_NAME, TIME_NAME, FLOOD_NAME, RISK_NAME):
        if v not in ds.variables:
            raise ValueError(f"Missing variable '{v}' in NetCDF. Found: {list(ds.variables)}")

    lat = _ensure_1d(ds[LAT_NAME], LAT_NAME)
    lon = _ensure_1d(ds[LON_NAME], LON_NAME)
    lat_edges = grid_edges_from_centers(lat)
    lon_edges = grid_edges_from_centers(lon)

    tsize = int(ds.sizes.get(TIME_NAME, 0))
    if tsize <= 0:
        raise ValueError("Invalid time dimension size.")
    if not (0 <= args.time_index < tsize):
        raise ValueError(f"--time-index out of range (time size={tsize})")

    flood2d = ds[FLOOD_NAME].isel({TIME_NAME: args.time_index})
    risk2d = ds[RISK_NAME].isel({TIME_NAME: args.time_index})

    # Materialize 2D slices once (avoid per-cell xarray indexing in hot loops).
    # Keep as float32 for cache efficiency; computations use float64 where needed.
    flood_np = np.asarray(flood2d.values, dtype=np.float32)
    risk_np = np.asarray(risk2d.values, dtype=np.float32)

    # ---- Fast zonal stats context (rasterization on the lon/lat grid) ----
    raster_transform = None
    raster_out_shape = None
    lat_rows = None

    if args.stats_mode in ("auto", "fast"):
        if rio_features is None:
            if args.stats_mode == "fast":
                raise RuntimeError("stats-mode=fast requires 'rasterio' (pip install rasterio).")
            LOG.info("rasterio not available -> falling back to exact zonal stats.")
        else:
            try:
                raster_transform, flip_lat = _build_affine_from_edges(lon_edges, lat_edges)
                raster_out_shape = (int(flood_np.shape[0]), int(flood_np.shape[1]))

                # Align arrays to raster row order (north->south).
                if flip_lat:
                    flood_np = np.flipud(flood_np)
                    risk_np = np.flipud(risk_np)
                    lat_rows = np.asarray(lat[::-1], dtype=np.float64)
                else:
                    lat_rows = np.asarray(lat, dtype=np.float64)

                LOG.info("Using fast zonal stats (rasterize) on %s x %s grid.", raster_out_shape[0], raster_out_shape[1])
            except Exception as exc:
                if args.stats_mode == "fast":
                    raise
                LOG.info("Cannot enable fast zonal stats (%s) -> falling back to exact.", exc)
                raster_transform = None
                raster_out_shape = None
                lat_rows = None


    fill_flood = ds[FLOOD_NAME].attrs.get("_FillValue", None)
    fill_risk = ds[RISK_NAME].attrs.get("_FillValue", None)

    LOG.info("Loading GeoJSON: %s", args.geojson_in)
    gj = _load_geojson_feature_collection(args.geojson_in) if rank == 0 else None

    if comm is not None:
        gj = comm.bcast(gj, root=0)

    feats = gj.get("features") or []

    my_indices = list(range(rank, len(feats), world_size))
    LOG.info("Rank %d handling %d features", rank, len(my_indices))

    results: List[Dict[str, Any]] = []
    progress_disable = rank != 0

    if threads_per_rank > 1 and len(my_indices) > 0:
        with ThreadPoolExecutor(max_workers=threads_per_rank) as ex:
            futs = {
                ex.submit(
                    _process_feature,
                    idx,
                    feats[idx],
                    args,
                    lat=lat,
                    lon=lon,
                    lat_edges=lat_edges,
                    lon_edges=lon_edges,
                    flood2d=flood2d,
                    risk2d=risk2d,
                    fill_flood=fill_flood,
                    fill_risk=fill_risk,
                    to4326=to4326,
                    to_area=to_area,
                    back_to_4326=back_to_4326,
                    risk_thr_tuple=risk_thr_tuple,
                    stats_mode=args.stats_mode,
                    raster_transform=raster_transform,
                    raster_out_shape=raster_out_shape,
                    flood_np=flood_np,
                    risk_np=risk_np,
                    lat_rows=lat_rows,
                ): idx
                for idx in my_indices
            }
            for fut in tqdm(as_completed(futs), total=len(futs), desc="Enriching features", unit="feature", disable=progress_disable):
                results.append(fut.result())
    else:
        for idx in tqdm(my_indices, desc="Enriching features", unit="feature", disable=progress_disable):
            results.append(
                _process_feature(
                    idx,
                    feats[idx],
                    args,
                    lat=lat,
                    lon=lon,
                    lat_edges=lat_edges,
                    lon_edges=lon_edges,
                    flood2d=flood2d,
                    risk2d=risk2d,
                    fill_flood=fill_flood,
                    fill_risk=fill_risk,
                    to4326=to4326,
                    to_area=to_area,
                    back_to_4326=back_to_4326,
                    risk_thr_tuple=risk_thr_tuple,
                    stats_mode=args.stats_mode,
                    raster_transform=raster_transform,
                    raster_out_shape=raster_out_shape,
                    flood_np=flood_np,
                    risk_np=risk_np,
                    lat_rows=lat_rows,
                )
            )

    gathered = [results]
    if comm is not None:
        gathered = comm.gather(results, root=0)
        comm.Barrier()
        if rank != 0:
            return 0

    merged_results: List[Dict[str, Any]] = []
    for chunk in gathered:
        merged_results.extend(chunk)

    for res in merged_results:
        feats[res["idx"]] = res["feature"]

    updated = sum(1 for r in merged_results if r["updated"])
    skipped = sum(1 for r in merged_results if r["skipped"])

    LOG.info("Updated=%d skipped=%d total=%d", updated, skipped, len(feats))

    with open(args.geojson_out, "w", encoding="utf-8") as f:
        json.dump(gj, f, ensure_ascii=False, indent=2)

    LOG.info("Wrote: %s", args.geojson_out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
