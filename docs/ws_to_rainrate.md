# Weather station CSV to rainfall NetCDF

`utils/ws_to_rainrate.py` converts weather-station rain-rate CSV files into a CF-compliant rainfall NetCDF file following `cdl/rain_time_dependent.cdl`. The script groups station observations by a fixed time interval, averages rain-rate samples per station, and spatially interpolates the resulting station field onto the model domain using 2D ordinary kriging.

## What the script produces

The output NetCDF includes:

- `rain_rate(time, latitude, longitude)` in **mm h-1**.
- A `crs` variable with CF grid-mapping metadata.
- `time` in **hours since 1900-01-01 00:00:0.0**.
- Global attributes describing the source, institution, and history.

## Requirements

- `pykrige` is required for kriging interpolation.
- `xarray` and `numpy` are used for NetCDF writing and metadata handling.

## Input format

Each CSV row must contain the following columns:

```
WEATHER_STATION_ID, ISO_8601_TIMESTAMP, LONGITUDE, LATITUDE, RAINRATE_VALUE
```

All timestamps should be ISO-8601 strings (UTC recommended, e.g., `2025-01-01T00:00:00Z`).

## Command-line usage

```bash
python utils/ws_to_rainrate.py \
  path/to/stations.csv \
  --output path/to/rain_time_dependent.nc \
  --domain path/to/domain.nc \
  --interval 600 \
  --source-name "Stations" \
  --institution "Your Institute" \
  --source "station gauge network" \
  --variogram-model spherical \
  --fill-value -9999 \
  --log-level INFO
```

### Arguments

- `input`: CSV file or directory containing CSV files (alphabetical order).
- `--output`: Output NetCDF path.
- `--domain`: Domain NetCDF (matching `cdl/domain.cdl`).
- `--interval`: Time interval in seconds used for grouping station samples.
- `--source-name`: Title prefix for the output dataset.
- `--institution`: Institution metadata.
- `--source`: Source metadata string.
- `--grid-mapping-name`: Override the CF `grid_mapping_name` attribute.
- `--epsg`: Override the EPSG code (e.g., `EPSG:4326`).
- `--fill-value`: Fill value for missing rain rates.
- `--variogram-model`: PyKrige variogram model (e.g., `linear`, `spherical`, `exponential`).
- `--log-level`: Logging level for console output.

## Example

To convert a folder of station CSVs into a time-dependent rainfall NetCDF aligned to the domain grid:

```bash
python utils/ws_to_rainrate.py \
  data/stations/ \
  --output data/rain/stations_rain.nc \
  --domain data/domain/domain.nc \
  --interval 900 \
  --source-name "Regional Gauges" \
  --institution "Hydro-Meteorological Center" \
  --source "gauges" \
  --variogram-model exponential
```

The resulting `stations_rain.nc` can be used as a rainfall forcing input in LPERFECT.
