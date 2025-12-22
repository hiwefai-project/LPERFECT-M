# Weather radar VMI GeoTIFF to rainfall NetCDF

`utils/wr_to_rain.py` converts a single-band weather radar VMI GeoTIFF (reflectivity in dBZ) into a CF-compliant rainfall NetCDF file that can be used as a rain forcing input for LPERFECT. The converter normalizes raster dimensions, applies a Z–R relationship to compute rain rates, attaches CF metadata, and writes a `rain_rate` variable with a time axis.

## What the script produces

The output NetCDF includes:

- `rain_rate(time, latitude, longitude)` in **mm h-1**.
- A `crs` variable with CF grid mapping metadata.
- `time` in **hours since 1900-01-01 00:00:0.0**.
- Global attributes describing the source, institution, and history.

## Requirements

- `rioxarray` and `rasterio` are required to read GeoTIFF inputs.
- `xarray` and `numpy` are used for NetCDF writing and metadata handling.

If you see an error about missing `rioxarray`, install the missing dependencies before running the script.

## Inputs and behavior

- **Input rasters:** comma-separated list of single-band GeoTIFFs (band 1 is used if multiple exist).
- **Dimensions:** accepts `x/y` or `latitude/longitude` and normalizes to `latitude/longitude`.
- **Fill values:** uses `_FillValue` or `missing_value` from the raster when present; otherwise uses the CLI `--fill-value` fallback.
- **Time:** use `--time` to provide an ISO-8601 timestamp for the first raster. The script steps forward by `--dt` seconds for each input.
- **Domain:** provide a domain NetCDF (`--domain`) that defines the target latitude/longitude grid.
- **Z–R parameters:** uses `Z = a * R^b` to convert dBZ to rain rate with configurable `a` and `b`.

## Converting dBZ to rain rate

The conversion follows the standard Z–R relationship:

```
Z = a * R^b
R = (Z / a)^(1 / b)
```

with `Z = 10^(dBZ / 10)`. The defaults (`a=200`, `b=1.6`) correspond to a Marshall–Palmer-style relationship often used for stratiform precipitation. You can supply alternate coefficients for local climatology or radar calibration via `--z-r-a` and `--z-r-b` (see the cited reference for guidance on appropriate ranges and tuned values).

## Command-line usage

```bash
python utils/wr_to_rain.py \
  --input PATH/TO/input_0000.tif,PATH/TO/input_0005.tif \
  --output PATH/TO/output.nc \
  --time 2024-06-01T12:00:00Z \
  --dt 300 \
  --domain PATH/TO/domain.nc \
  --source-name "Radar" \
  --institution "Your Institute" \
  --source "Radar VMI" \
  --grid-mapping-name latitude_longitude \
  --epsg EPSG:4326 \
  --fill-value -9999 \
  --z-r-a 200 \
  --z-r-b 1.6 \
  --log-level INFO
```

### Arguments

- `--input`: Comma-separated list of input VMI GeoTIFFs, in time order.
- `--output`: Path for the output NetCDF.
- `--time`: ISO-8601 timestamp for the first raster.
- `--dt`: Seconds between consecutive inputs.
- `--domain`: Domain NetCDF (matching `cdl/domain.cdl`) used for regridding.
- `--source-name`: Title prefix for the output dataset (defaults to `Radar`).
- `--institution`: Institution metadata (defaults to `Unknown`).
- `--source`: Source metadata (defaults to `radar`).
- `--grid-mapping-name`: Override the CF `grid_mapping_name` attribute.
- `--epsg`: Override the EPSG code (e.g., `EPSG:4326`).
- `--fill-value`: Fallback fill value when the raster provides none.
- `--z-r-a`: Z–R parameter `a` in `Z=a*R^b` (defaults to `200.0`).
- `--z-r-b`: Z–R parameter `b` in `Z=a*R^b` (defaults to `1.6`).
- `--log-level`: Logging level for console output (default: `INFO`).

## Example

Convert a radar GeoTIFF into a rain forcing file with explicit metadata:

```bash
python utils/wr_to_rain.py \
  --input data/radar/vmi_20240601_1200.tif,data/radar/vmi_20240601_1205.tif \
  --output data/rain/rain_20240601_1200.nc \
  --time 2024-06-01T12:00:00Z \
  --dt 300 \
  --domain data/domain/domain.nc \
  --source-name "Regional Radar" \
  --institution "Hydro-Meteorological Center" \
  --source "VMI radar product" \
  --epsg EPSG:4326 \
  --z-r-a 300 \
  --z-r-b 1.4
```

The resulting `rain_20240601_1200.nc` is ready to be referenced in your LPERFECT configuration as a rainfall forcing dataset.
