# `wr_to_rain.py`: Weather radar VMI GeoTIFF to rainfall NetCDF

`utils/wr_to_rain.py` converts a single-band weather radar VMI GeoTIFF into a CF-compliant rainfall NetCDF file that can be used as a rain forcing input for LPERFECT. The converter normalizes raster dimensions, attaches CF metadata, and writes a `rain_rate` variable with a time axis.

## What the script produces

The output NetCDF includes:

- `rain_rate(time, latitude, longitude)` in **mm h-1**.
- A `crs` variable with CF grid mapping metadata.
- `time` in **hours since 1900-01-01 00:00:00 UTC**.
- Global attributes describing the source, institution, and history.

## Requirements

- `rioxarray` and `rasterio` are required to read GeoTIFF inputs.
- `xarray` and `numpy` are used for NetCDF writing and metadata handling.

If you see an error about missing `rioxarray`, install the missing dependencies before running the script.

## Inputs and behavior

- **Input raster:** a single-band GeoTIFF (band 1 is used if multiple exist).
- **Dimensions:** accepts `x/y` or `latitude/longitude` and normalizes to `latitude/longitude`.
- **Fill values:** uses `_FillValue` or `missing_value` from the raster when present; otherwise uses the CLI `--fill-value` fallback.
- **Time:** use `--time` to provide an ISO-8601 timestamp. If omitted, the current UTC time is used.

## Command-line usage

```bash
python utils/wr_to_rain.py \
  --input PATH/TO/input.tif \
  --output PATH/TO/output.nc \
  --time 2024-06-01T12:00:00Z \
  --source-name "Radar" \
  --institution "Your Institute" \
  --source "Radar VMI" \
  --grid-mapping-name latitude_longitude \
  --epsg EPSG:4326 \
  --fill-value -9999
```

### Arguments

- `--input`: Path to the input VMI GeoTIFF.
- `--output`: Path for the output NetCDF.
- `--time`: ISO-8601 timestamp (defaults to current UTC if omitted).
- `--source-name`: Title prefix for the output dataset (defaults to `Radar`).
- `--institution`: Institution metadata (defaults to `Unknown`).
- `--source`: Source metadata (defaults to `radar`).
- `--grid-mapping-name`: Override the CF `grid_mapping_name` attribute.
- `--epsg`: Override the EPSG code (e.g., `EPSG:4326`).
- `--fill-value`: Fallback fill value when the raster provides none.

## Example

Convert a radar GeoTIFF into a rain forcing file with explicit metadata:

```bash
python utils/wr_to_rain.py \
  --input data/radar/vmi_20240601_1200.tif \
  --output data/rain/rain_20240601_1200.nc \
  --time 2024-06-01T12:00:00Z \
  --source-name "Regional Radar" \
  --institution "Hydro-Meteorological Center" \
  --source "VMI radar product" \
  --epsg EPSG:4326
```

The resulting `rain_20240601_1200.nc` is ready to be referenced in your LPERFECT configuration as a rainfall forcing dataset.
