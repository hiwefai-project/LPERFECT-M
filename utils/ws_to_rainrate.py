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
from typing import Iterable, Iterator, Sequence

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


def _station_bounding_box(
    records: Sequence[tuple[str, datetime, float, float, float]],
    buffer_meters: float,
) -> tuple[float, float, float, float]:
    """Compute a lon/lat bounding box around stations with a meter buffer."""
    if buffer_meters < 0:
        raise ValueError("Buffer meters must be non-negative.")
    # Collect station coordinates so we can compute extrema.
    lons = np.array([item[2] for item in records], dtype=float)
    lats = np.array([item[3] for item in records], dtype=float)
    if lons.size == 0 or lats.size == 0:
        raise ValueError("No station coordinates available for bounding box.")
    # Convert meter buffer to degrees using mean latitude scaling.
    mean_lat_radians = math.radians(float(np.mean(lats)))
    cosine_lat = max(math.cos(mean_lat_radians), 1e-6)
    buffer_degrees_lat = buffer_meters / METERS_PER_DEGREE_LAT
    buffer_degrees_lon = buffer_meters / (METERS_PER_DEGREE_LAT * cosine_lat)
    # Expand the station extent by the buffer in both directions.
    min_lon = float(np.min(lons) - buffer_degrees_lon)
    max_lon = float(np.max(lons) + buffer_degrees_lon)
    min_lat = float(np.min(lats) - buffer_degrees_lat)
    max_lat = float(np.max(lats) + buffer_degrees_lat)
    return min_lon, max_lon, min_lat, max_lat


def _krige_grid(
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    values: np.ndarray,
    grid_longitudes: np.ndarray,
    grid_latitudes: np.ndarray,
    variogram_model: str,
) -> np.ndarray:
    """Perform ordinary kriging over the domain grid."""
    _require_pykrige()
    from pykrige.ok import OrdinaryKriging

    kriging = OrdinaryKriging(
        longitudes,
        latitudes,
        values,
        variogram_model=variogram_model,
        verbose=False,
        enable_plotting=False,
    )
    grid, _ = kriging.execute("grid", grid_longitudes, grid_latitudes)
    return np.asarray(grid, dtype=float)


def _idw_grid(
    longitudes: np.ndarray,
    latitudes: np.ndarray,
    values: np.ndarray,
    grid_longitudes: np.ndarray,
    grid_latitudes: np.ndarray,
    power: float = 2.0,
    min_distance_meters: float = 1.0,
) -> np.ndarray:
    """Perform inverse-distance weighting over the domain grid."""
    # Prepare a latitude/longitude mesh for vectorized distance math.
    lon_mesh = grid_longitudes[np.newaxis, :]
    lat_mesh = grid_latitudes[:, np.newaxis]
    # Allocate arrays for weighted sums and total weights.
    weighted_sum = np.zeros((grid_latitudes.size, grid_longitudes.size), dtype=float)
    weight_total = np.zeros_like(weighted_sum)
    # Walk each station so we can accumulate IDW contributions.
    for lon, lat, value in zip(longitudes, latitudes, values):
        # Convert longitudinal degrees into meters for this station's latitude.
        meters_per_degree_lon = METERS_PER_DEGREE_LAT * max(math.cos(math.radians(lat)), 1e-6)
        # Convert degree offsets into meter offsets along each axis.
        dx = (lon_mesh - lon) * meters_per_degree_lon
        dy = (lat_mesh - lat) * METERS_PER_DEGREE_LAT
        # Compute straight-line distances in meters.
        distances = np.hypot(dx, dy)
        # Identify grid points that coincide with the station location.
        exact_mask = distances == 0.0
        # Force exact matches to the station value.
        if exact_mask.any():
            weighted_sum[exact_mask] = value
            weight_total[exact_mask] = np.inf
        # Apply a minimum distance to avoid divide-by-zero in weights.
        safe_distances = np.maximum(distances, min_distance_meters)
        # Compute inverse-distance weights for the station.
        weights = 1.0 / np.power(safe_distances, power)
        # Skip updates where an exact match was already set.
        weights = np.where(weight_total == np.inf, 0.0, weights)
        # Accumulate weighted values and weights.
        weighted_sum += weights * value
        weight_total += weights
    # Compute the weighted average for non-exact cells.
    interpolated = np.divide(
        weighted_sum,
        weight_total,
        out=np.zeros_like(weighted_sum),
        where=(weight_total > 0) & (weight_total != np.inf),
    )
    # Restore exact station values where we locked them in.
    interpolated = np.where(weight_total == np.inf, weighted_sum, interpolated)
    return interpolated


def _interpolate_interval(
    aggregates: dict[str, StationAggregate],
    grid_longitudes: np.ndarray,
    grid_latitudes: np.ndarray,
    variogram_model: str,
    fill_value: float,
    bounding_box: tuple[float, float, float, float] | None,
    outside_value: float,
) -> np.ndarray:
    """Interpolate station rain rates for a single time interval."""
    if bounding_box is not None:
        # Limit the kriging grid to the buffered station bounding box.
        min_lon, max_lon, min_lat, max_lat = bounding_box
        lon_mask = (grid_longitudes >= min_lon) & (grid_longitudes <= max_lon)
        lat_mask = (grid_latitudes >= min_lat) & (grid_latitudes <= max_lat)
        if not lon_mask.any() or not lat_mask.any():
            return np.full((grid_latitudes.size, grid_longitudes.size), outside_value, dtype=float)
        # Extract the bounded longitude/latitude subgrid for kriging.
        sub_longitudes = grid_longitudes[lon_mask]
        sub_latitudes = grid_latitudes[lat_mask]
    else:
        lon_mask = None
        lat_mask = None
        sub_longitudes = grid_longitudes
        sub_latitudes = grid_latitudes
    if not aggregates:
        # No stations in this interval: fill the kriging area with the fill value.
        interval_grid = np.full((sub_latitudes.size, sub_longitudes.size), fill_value, dtype=float)
    else:
        lons = np.array([agg.longitude for agg in aggregates.values()], dtype=float)
        lats = np.array([agg.latitude for agg in aggregates.values()], dtype=float)
        vals = np.array([agg.mean for agg in aggregates.values()], dtype=float)
        # Filter out non-finite station values before interpolation.
        valid_mask = np.isfinite(vals)
        lons = lons[valid_mask]
        lats = lats[valid_mask]
        vals = vals[valid_mask]
        if lons.size == 0:
            # Fall back to the fill value when all station values are invalid.
            interval_grid = np.full(
                (sub_latitudes.size, sub_longitudes.size),
                fill_value,
                dtype=float,
            )
        elif lons.size == 1:
            LOGGER.warning("Only one station for interpolation; using station value fill.")
            # Fill the kriging area with the single station value.
            interval_grid = np.full(
                (sub_latitudes.size, sub_longitudes.size),
                float(vals[0]),
                dtype=float,
            )
        elif lons.size == 2:
            LOGGER.warning("Only two stations for interpolation; using IDW fallback.")
            # Use inverse-distance weighting to provide a spatial gradient.
            interval_grid = _idw_grid(
                lons,
                lats,
                vals,
                sub_longitudes,
                sub_latitudes,
            )
        else:
            # Run kriging on the reduced grid and replace NaNs with the fill value.
            interval_grid = _krige_grid(
                lons,
                lats,
                vals,
                sub_longitudes,
                sub_latitudes,
                variogram_model,
            )
            interval_grid = np.where(np.isfinite(interval_grid), interval_grid, fill_value)
    if bounding_box is None:
        return interval_grid.astype(np.float32)
    # Insert the kriged subgrid into the full domain with outside values set to zero.
    full_grid = np.full((grid_latitudes.size, grid_longitudes.size), outside_value, dtype=float)
    full_grid[np.ix_(lat_mask, lon_mask)] = interval_grid
    return full_grid.astype(np.float32)


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
        help="Apply a station-based bounding box buffered by this many meters.",
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

    bounding_box = None
    if args.buffer_meters is not None:
        # Build a bounding box around all station locations with a meter buffer.
        bounding_box = _station_bounding_box(records, args.buffer_meters)
        LOGGER.info(
            "Applying station buffer %.2f m -> bbox lon[%.4f, %.4f], lat[%.4f, %.4f]",
            args.buffer_meters,
            bounding_box[0],
            bounding_box[1],
            bounding_box[2],
            bounding_box[3],
        )

    LOGGER.info("Loading domain grid from %s", args.domain)
    domain = _load_domain(args.domain)
    latitudes, longitudes = _domain_coords(domain)
    crs_attrs = _domain_crs_attrs(domain)
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
            args.variogram_model,
            args.fill_value,
            bounding_box,
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
