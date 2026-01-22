#!/usr/bin/env python3
"""Convert weather-station CSV data into CF-compliant rainfall NetCDF."""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from datetime import datetime, timezone
import importlib.util
import logging
import math
from pathlib import Path
from typing import Iterable, Iterator, Sequence, Optional

import numpy as np
import xarray as xr


CF_CONVENTIONS = "CF-1.10"
TIME_UNITS = "hours since 1900-01-01 00:00:0.0"
LOGGER = logging.getLogger(__name__)
METERS_PER_DEGREE_LAT = 111_320.0


@dataclass
class StationAggregate:
    """Aggregate rain-rate samples for a station within a time interval."""

    station_id: str
    longitude: float
    latitude: float
    total: float
    count: int

    def add(self, value: float) -> None:
        """Add a rain-rate sample to the aggregate."""
        self.total += value
        self.count += 1

    @property
    def mean(self) -> float:
        """Return the average rain rate for the interval."""
        return self.total / self.count if self.count > 0 else float("nan")


def _setup_logging(level: str) -> None:
    """Configure logging for console output."""
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="%(asctime)s %(levelname)s %(message)s",
    )


def _require_pykrige() -> None:
    """Ensure PyKrige is available before running kriging interpolation."""
    if importlib.util.find_spec("pykrige") is None:
        raise ImportError(
            "PyKrige is required for kriging interpolation. "
            "Install it with 'pip install pykrige'."
        )



def _require_pyproj_shapely() -> None:
    """Ensure pyproj and shapely are available for metric buffering/masking."""
    if importlib.util.find_spec("pyproj") is None:
        raise ImportError(
            "pyproj is required for metric buffering/masking. "
            "Install it with 'pip install pyproj'."
        )
    if importlib.util.find_spec("shapely") is None:
        raise ImportError(
            "shapely is required for metric buffering/masking. "
            "Install it with 'pip install shapely'."
        )


def _auto_utm_epsg(longitudes: np.ndarray, latitudes: np.ndarray) -> int:
    """Pick a suitable UTM EPSG code based on the mean station location."""
    mean_lon = float(np.mean(longitudes))
    mean_lat = float(np.mean(latitudes))
    zone = int(math.floor((mean_lon + 180.0) / 6.0) + 1)
    zone = min(max(zone, 1), 60)
    if mean_lat >= 0.0:
        return 32600 + zone  # WGS84 / UTM Northern Hemisphere
    return 32700 + zone  # WGS84 / UTM Southern Hemisphere


def _project_lonlat_to_xy(
    lons: np.ndarray,
    lats: np.ndarray,
    epsg: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Project lon/lat arrays to metric x/y arrays using the provided EPSG."""
    _require_pyproj_shapely()
    from pyproj import Transformer

    transformer = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
    x, y = transformer.transform(lons, lats)
    return np.asarray(x, dtype=float), np.asarray(y, dtype=float)


def _build_hull_buffer_mask(
    station_lons: np.ndarray,
    station_lats: np.ndarray,
    grid_lons: np.ndarray,
    grid_lats: np.ndarray,
    buffer_meters: float,
    epsg: Optional[int] = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    """Return (mask2d, x2d, y2d, used_epsg) for the convex-hull+buffer domain."""
    if buffer_meters < 0:
        raise ValueError("Buffer meters must be non-negative.")
    _require_pyproj_shapely()
    from shapely.geometry import MultiPoint, Point
    from shapely.prepared import prep

    if epsg is None:
        epsg = _auto_utm_epsg(station_lons, station_lats)

    # Project stations and grid to metric coordinates.
    station_x, station_y = _project_lonlat_to_xy(station_lons, station_lats, epsg)
    lon_mesh, lat_mesh = np.meshgrid(grid_lons, grid_lats)
    grid_x, grid_y = _project_lonlat_to_xy(lon_mesh.ravel(), lat_mesh.ravel(), epsg)
    grid_x2d = grid_x.reshape(lon_mesh.shape)
    grid_y2d = grid_y.reshape(lon_mesh.shape)

    # Build convex hull and apply metric buffer.
    hull = MultiPoint(np.column_stack([station_x, station_y])).convex_hull
    buffered = hull.buffer(buffer_meters)

    prepared = prep(buffered)

    # Build mask. Prefer vectorized predicates when available.
    try:
        # Shapely 2.x
        from shapely import contains_xy  # type: ignore

        mask_flat = contains_xy(buffered, grid_x, grid_y)
        mask2d = np.asarray(mask_flat, dtype=bool).reshape(lon_mesh.shape)
    except Exception:
        # Generic (works on Shapely 1.8+): loop over points.
        mask2d = np.zeros(lon_mesh.shape, dtype=bool)
        it = np.nditer(mask2d, flags=["multi_index"], op_flags=["readwrite"])
        for cell in it:
            j, i = it.multi_index
            cell[...] = prepared.contains(
                Point(float(grid_x2d[j, i]), float(grid_y2d[j, i]))
            )

    return mask2d, grid_x2d, grid_y2d, int(epsg)



def _parse_time(value: str) -> datetime:
    """Parse ISO-8601 timestamps and normalize to UTC."""
    text = value.strip()
    if text.endswith("Z"):
        text = text[:-1] + "+00:00"
    dt = datetime.fromisoformat(text)
    if dt.tzinfo is None:
        dt = dt.replace(tzinfo=timezone.utc)
    return dt.astimezone(timezone.utc)


def _time_to_seconds(dt: datetime) -> float:
    """Convert a datetime to seconds since 1900-01-01 UTC."""
    base = datetime(1900, 1, 1, tzinfo=timezone.utc)
    return (dt - base).total_seconds()


def _seconds_to_hours(seconds: float) -> float:
    """Convert seconds since 1900-01-01 to hours for CF time units."""
    return seconds / 3600.0


def _input_paths(value: str) -> list[Path]:
    """Resolve input paths from a file or directory argument."""
    path = Path(value)
    if not path.exists():
        raise FileNotFoundError(f"Input path not found: {path}")
    if path.is_dir():
        csv_files = sorted([child for child in path.iterdir() if child.is_file()])
        if not csv_files:
            raise ValueError(f"No CSV files found in directory: {path}")
        return csv_files
    if path.is_file():
        return [path]
    raise ValueError(f"Unsupported input path: {path}")


def _load_csv_rows(paths: Sequence[Path]) -> list[dict[str, str]]:
    """Load CSV rows from one or more CSV files."""
    rows: list[dict[str, str]] = []
    for csv_path in paths:
        with csv_path.open("r", newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                raise ValueError(f"CSV file has no header: {csv_path}")
            rows.extend(reader)
    if not rows:
        raise ValueError("No CSV rows found in input files.")
    return rows


def _iter_records(rows: Iterable[dict[str, str]]) -> Iterator[tuple[str, datetime, float, float, float]]:
    """Yield parsed records from CSV row dictionaries."""
    required = {
        "WEATHER_STATION_ID",
        "ISO_8601_TIMESTAMP",
        "LONGITUDE",
        "LATITUDE",
        "RAINRATE_VALUE",
    }
    for row in rows:
        missing = required - set(row)
        if missing:
            raise ValueError(f"CSV row missing columns: {sorted(missing)}")
        station_id = row["WEATHER_STATION_ID"].strip()
        timestamp = _parse_time(row["ISO_8601_TIMESTAMP"])
        longitude = float(row["LONGITUDE"])
        latitude = float(row["LATITUDE"])
        rainrate = float(row["RAINRATE_VALUE"])
        yield station_id, timestamp, longitude, latitude, rainrate


def _build_time_bins(
    records: Iterable[tuple[str, datetime, float, float, float]],
    interval_seconds: float,
) -> tuple[float, dict[int, dict[str, StationAggregate]]]:
    """Group station samples into interval bins with station-level averages."""
    if interval_seconds <= 0:
        raise ValueError("Interval seconds must be positive.")
    entries = list(records)
    if not entries:
        raise ValueError("No valid station records to process.")
    seconds_since = np.array([_time_to_seconds(item[1]) for item in entries], dtype=float)
    start_seconds = math.floor(seconds_since.min() / interval_seconds) * interval_seconds
    bins: dict[int, dict[str, StationAggregate]] = {}
    for station_id, timestamp, longitude, latitude, rainrate in entries:
        if not math.isfinite(rainrate):
            LOGGER.warning("Skipping non-finite rainrate for station %s", station_id)
            continue
        seconds = _time_to_seconds(timestamp)
        bin_index = int((seconds - start_seconds) // interval_seconds)
        bucket = bins.setdefault(bin_index, {})
        aggregate = bucket.get(station_id)
        if aggregate is None:
            aggregate = StationAggregate(
                station_id=station_id,
                longitude=longitude,
                latitude=latitude,
                total=0.0,
                count=0,
            )
            bucket[station_id] = aggregate
        if not math.isclose(aggregate.longitude, longitude) or not math.isclose(aggregate.latitude, latitude):
            LOGGER.warning("Station %s has changing coordinates; using first occurrence.", station_id)
        aggregate.add(rainrate)
    if not bins:
        raise ValueError("All rain-rate samples were invalid or filtered out.")
    return start_seconds, bins


def _load_domain(domain_path: str) -> xr.Dataset:
    """Load the domain NetCDF for latitude/longitude coordinates."""
    resolved = Path(domain_path)
    if not resolved.exists():
        raise FileNotFoundError(f"Domain NetCDF not found: {resolved}")
    return xr.open_dataset(resolved)


def _domain_coords(domain: xr.Dataset) -> tuple[np.ndarray, np.ndarray]:
    """Extract latitude and longitude arrays from the domain dataset."""
    if "latitude" not in domain or "longitude" not in domain:
        raise ValueError("Domain dataset must define latitude and longitude variables.")
    lat_values = np.asarray(domain["latitude"].values, dtype=float)
    lon_values = np.asarray(domain["longitude"].values, dtype=float)
    if lat_values.ndim != 1 or lon_values.ndim != 1:
        raise ValueError("Domain latitude and longitude must be 1D arrays.")
    return lat_values, lon_values


def _domain_crs_attrs(domain: xr.Dataset) -> dict[str, str]:
    """Extract CRS attributes from the domain, if present."""
    if "crs" not in domain:
        return {}
    return {str(key): str(value) for key, value in domain["crs"].attrs.items()}


def _krige_points(
    station_x: np.ndarray,
    station_y: np.ndarray,
    station_values: np.ndarray,
    points_x: np.ndarray,
    points_y: np.ndarray,
    variogram_model: str,
) -> np.ndarray:
    """Perform ordinary kriging for arbitrary points in metric coordinates."""
    _require_pykrige()
    from pykrige.ok import OrdinaryKriging

    ok = OrdinaryKriging(
        station_x,
        station_y,
        station_values,
        variogram_model=variogram_model,
        verbose=False,
        enable_plotting=False,
    )
    z, _ = ok.execute("points", points_x, points_y)
    return np.asarray(z, dtype=float)


def _idw_points(
    station_x: np.ndarray,
    station_y: np.ndarray,
    station_values: np.ndarray,
    points_x: np.ndarray,
    points_y: np.ndarray,
    power: float = 2.0,
    min_distance_meters: float = 1.0,
) -> np.ndarray:
    """Inverse-distance weighting for arbitrary points in metric coordinates."""
    station_x = np.asarray(station_x, dtype=float)
    station_y = np.asarray(station_y, dtype=float)
    station_values = np.asarray(station_values, dtype=float)

    px = np.asarray(points_x, dtype=float)
    py = np.asarray(points_y, dtype=float)

    # Compute distances to each station (n_points, n_stations).
    dx = px[:, None] - station_x[None, :]
    dy = py[:, None] - station_y[None, :]
    dist = np.hypot(dx, dy)

    # Handle exact matches.
    exact = dist == 0.0
    out = np.empty(px.shape[0], dtype=float)
    if exact.any():
        # For each point, if it matches any station, take the first match value.
        first_match = exact.argmax(axis=1)
        has_match = exact.any(axis=1)
        out[has_match] = station_values[first_match[has_match]]
        # For non-matching points proceed with IDW.
        use = ~has_match
    else:
        use = np.ones(px.shape[0], dtype=bool)

    if use.any():
        safe = np.maximum(dist[use], min_distance_meters)
        w = 1.0 / np.power(safe, power)
        wsum = np.sum(w, axis=1)
        out[use] = np.sum(w * station_values[None, :], axis=1) / wsum

    return out


def _interpolate_interval(
    aggregates: dict[str, StationAggregate],
    grid_longitudes: np.ndarray,
    grid_latitudes: np.ndarray,
    grid_x2d: np.ndarray,
    grid_y2d: np.ndarray,
    domain_mask2d: Optional[np.ndarray],
    epsg: int,
    variogram_model: str,
    fill_value: float,
    outside_value: float,
) -> np.ndarray:
    """Interpolate station rain rates for a single time interval.

    Implements:
      - Convex-hull + metric buffer masking (if domain_mask2d provided)
      - Ordinary kriging on log1p(rain_rate) with back-transform
      - Fallbacks: 2 stations -> IDW, 1 station -> constant, 0 -> fill_value
    """
    out = np.full((grid_latitudes.size, grid_longitudes.size), outside_value, dtype=float)

    if domain_mask2d is None:
        mask = np.ones_like(out, dtype=bool)
    else:
        if domain_mask2d.shape != out.shape:
            raise ValueError("Domain mask shape does not match the domain grid.")
        mask = domain_mask2d

    if not mask.any():
        return out.astype(np.float32)

    if not aggregates:
        out[mask] = fill_value
        return out.astype(np.float32)

    lons = np.array([agg.longitude for agg in aggregates.values()], dtype=float)
    lats = np.array([agg.latitude for agg in aggregates.values()], dtype=float)
    vals = np.array([agg.mean for agg in aggregates.values()], dtype=float)

    valid = np.isfinite(vals)
    lons = lons[valid]
    lats = lats[valid]
    vals = vals[valid]

    if lons.size == 0:
        out[mask] = fill_value
        return out.astype(np.float32)

    # Project station coordinates into the same metric CRS used for the grid.
    # We infer the CRS by reusing the grid_x2d/grid_y2d generation path:
    # it guarantees station/grid were projected with the same EPSG in main().
    # Here we approximate by projecting using local meters-per-degree only if needed.
    # NOTE: main() passes projected station coords by re-projecting lon/lat to the same EPSG.
    # To avoid duplicating logic, we project again using a local automatic UTM.
    station_x, station_y = _project_lonlat_to_xy(lons, lats, epsg)

    # Points to interpolate (only inside mask).
    px = grid_x2d[mask].ravel()
    py = grid_y2d[mask].ravel()

    if lons.size == 1:
        out[mask] = float(vals[0])
        return out.astype(np.float32)

    if lons.size == 2:
        LOGGER.warning("Only two stations for interpolation; using IDW fallback.")
        out_vals = _idw_points(station_x, station_y, vals, px, py)
        out[mask] = np.where(np.isfinite(out_vals), out_vals, fill_value)
        return out.astype(np.float32)

    # Kriging on transformed values.
    vals_pos = np.clip(vals, 0.0, None)
    z = np.log1p(vals_pos)

    z_hat = _krige_points(station_x, station_y, z, px, py, variogram_model)
    r_hat = np.expm1(z_hat)
    r_hat = np.where(np.isfinite(r_hat), r_hat, fill_value)
    r_hat = np.clip(r_hat, 0.0, None)

    out[mask] = r_hat
    return out.astype(np.float32)


def _build_rain_dataset(
    time_hours: np.ndarray,
    rain_rate: np.ndarray,
    latitudes: np.ndarray,
    longitudes: np.ndarray,
    crs_attrs: dict[str, str],
    fill_value: float,
    source_name: str,
    institution: str,
    source: str,
    history: str,
) -> xr.Dataset:
    """Construct the CF-compliant rainfall dataset."""
    dataset = xr.Dataset()
    dataset["time"] = ("time", time_hours.astype(np.float64))
    dataset["time"].attrs.update(
        {
            "description": "Time",
            "long_name": "time",
            "units": TIME_UNITS,
        }
    )
    dataset["latitude"] = ("latitude", latitudes.astype(np.float32))
    dataset["latitude"].attrs.update(
        {
            "description": "Latitude",
            "long_name": "latitude",
            "units": "degrees_north",
        }
    )
    dataset["longitude"] = ("longitude", longitudes.astype(np.float32))
    dataset["longitude"].attrs.update(
        {
            "description": "Longitude",
            "long_name": "longitude",
            "units": "degrees_east",
        }
    )
    dataset["crs"] = xr.DataArray(0, attrs=crs_attrs)
    dataset["rain_rate"] = (
        ("time", "latitude", "longitude"),
        rain_rate.astype(np.float32),
    )
    dataset["rain_rate"].attrs.update(
        {
            "long_name": "rainfall_rate",
            "standard_name": "rainfall_rate",
            "units": "mm h-1",
            "grid_mapping": "crs",
            "_FillValue": fill_value,
        }
    )
    dataset.attrs.update(
        {
            "Conventions": CF_CONVENTIONS,
            "title": f"{source_name} rainfall forcing",
            "institution": institution,
            "source": source,
            "history": history,
        }
    )
    return dataset


def _build_history() -> str:
    """Create a history string for the NetCDF output."""
    timestamp = datetime.now(timezone.utc).isoformat()
    return f"{timestamp}: produced/ingested"


def main() -> None:
    """CLI entry point."""
    parser = argparse.ArgumentParser(
        description="Convert weather station CSV data into CF-compliant rainfall NetCDF.",
    )
    parser.add_argument(
        "input",
        help="CSV file or directory with CSV files.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output NetCDF path.",
    )
    parser.add_argument(
        "--domain",
        required=True,
        help="Domain NetCDF path (matches cdl/domain.cdl).",
    )
    parser.add_argument(
        "--interval",
        type=float,
        required=True,
        help="Time interval in seconds for grouping station data.",
    )
    parser.add_argument(
        "--source-name",
        default="Stations",
        help="Source name for the output title.",
    )
    parser.add_argument(
        "--institution",
        default="Unknown",
        help="Institution metadata for the output NetCDF.",
    )
    parser.add_argument(
        "--source",
        default="stations",
        help="Source metadata string.",
    )
    parser.add_argument(
        "--grid-mapping-name",
        default=None,
        help="Override grid_mapping_name for the CRS variable.",
    )
    parser.add_argument(
        "--epsg",
        default=None,
        help="Override EPSG code for the CRS variable (e.g., EPSG:4326).",
    )
    parser.add_argument(
        "--fill-value",
        type=float,
        default=-9999.0,
        help="Fill value for missing rain rates.",
    )
    parser.add_argument(
        "--variogram-model",
        default="spherical",
        help="PyKrige variogram model (e.g., linear, spherical, exponential).",
    )
    parser.add_argument(
        "--buffer-meters",
        type=float,
        default=None,
        help="Limit interpolation to the station convex hull buffered by this many meters (metric buffer).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        help="Logging level (default: INFO).",
    )
    args = parser.parse_args()

    _setup_logging(args.log_level)

    LOGGER.info("Reading input CSVs from %s", args.input)
    input_paths = _input_paths(args.input)
    rows = _load_csv_rows(input_paths)

    LOGGER.info("Parsing weather station records")
    records = list(_iter_records(rows))

    LOGGER.info("Grouping records into %s-second intervals", args.interval)
    start_seconds, bins = _build_time_bins(records, args.interval)

    LOGGER.info("Loading domain grid from %s", args.domain)
    domain = _load_domain(args.domain)
    latitudes, longitudes = _domain_coords(domain)
    crs_attrs = _domain_crs_attrs(domain)
    # Pre-compute station hull+buffer mask (and metric grid coordinates) in a UTM CRS.
    station_lonlat: dict[str, tuple[float, float]] = {}
    for station_id, _, lon, lat, _ in records:
        station_lonlat.setdefault(station_id, (lon, lat))
    station_lons = np.array([item[0] for item in station_lonlat.values()], dtype=float)
    station_lats = np.array([item[1] for item in station_lonlat.values()], dtype=float)
    if station_lons.size == 0:
        raise ValueError("No station coordinates found in input records.")

    epsg_used = _auto_utm_epsg(station_lons, station_lats)

    if args.buffer_meters is not None:
        LOGGER.info("Building convex-hull+buffer mask (buffer=%.2f m, EPSG:%d)", args.buffer_meters, epsg_used)
        domain_mask2d, grid_x2d, grid_y2d, epsg_used = _build_hull_buffer_mask(
            station_lons=station_lons,
            station_lats=station_lats,
            grid_lons=longitudes,
            grid_lats=latitudes,
            buffer_meters=args.buffer_meters,
            epsg=epsg_used,
        )
    else:
        LOGGER.info("No buffer provided; interpolating on the full domain (EPSG:%d for metric calculations).", epsg_used)
        lon_mesh, lat_mesh = np.meshgrid(longitudes, latitudes)
        grid_x, grid_y = _project_lonlat_to_xy(lon_mesh.ravel(), lat_mesh.ravel(), epsg_used)
        grid_x2d = grid_x.reshape(lon_mesh.shape)
        grid_y2d = grid_y.reshape(lon_mesh.shape)
        domain_mask2d = None

    if args.grid_mapping_name is not None:
        crs_attrs["grid_mapping_name"] = args.grid_mapping_name
    if args.epsg is not None:
        crs_attrs["epsg_code"] = args.epsg

    LOGGER.info("Interpolating %s time intervals", len(bins))
    time_indices = sorted(bins.keys())
    time_hours: list[float] = []
    rain_fields: list[np.ndarray] = []
    for index in time_indices:
        interval_seconds = start_seconds + index * args.interval
        time_hours.append(_seconds_to_hours(interval_seconds))
        field = _interpolate_interval(
            bins[index],
            longitudes,
            latitudes,
            grid_x2d,
            grid_y2d,
            domain_mask2d,
            epsg_used,
            args.variogram_model,
            args.fill_value,
            outside_value=0.0,
        )
        rain_fields.append(field)

    rain_rate = np.stack(rain_fields, axis=0)
    dataset = _build_rain_dataset(
        time_hours=np.array(time_hours, dtype=float),
        rain_rate=rain_rate,
        latitudes=latitudes,
        longitudes=longitudes,
        crs_attrs=crs_attrs,
        fill_value=args.fill_value,
        source_name=args.source_name,
        institution=args.institution,
        source=args.source,
        history=_build_history(),
    )

    output_path = Path(args.output)
    LOGGER.info("Writing rainfall NetCDF to %s", output_path)
    dataset.to_netcdf(output_path)


if __name__ == "__main__":
    main()
