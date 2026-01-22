# meteo@uniparthenope to rain rate

This script converts **one or many** WRF-derived NetCDF files (e.g. `wrf5_d02_*.nc`) into a **single merged time series** NetCDF compliant with a simplified `cdl/rain_time_dependent.cdl` schema.

It is intended for hydrological / impact models that require **rain rate** instead of accumulated precipitation.

It now always **regrids the output onto a reference domain NetCDF** built from `cdl/domain.cdl`, so the merged time series matches your model grid.

---

## What it does

### Input (WRF-derived)
WRF post-processing products often provide:
- **Accumulated precipitation** over a time window (commonly 1 hour)
- Units: **mm**

The rain accumulation must be provided as **`DELTA_RAIN`**.

### Output (`rain_time_dependent`)
The output contains:
- `rain_rate(time, latitude, longitude)` in **kg m⁻² s⁻¹**
- Latitude/longitude and metadata are taken from the supplied **domain NetCDF** (matching `cdl/domain.cdl`).

### New: merge multiple times into one time series
You can pass multiple file paths or globs. The script will:
1. Convert each file to `rain_rate`
2. Concatenate along `time`
3. Sort by time (default)
4. Deduplicate identical times (default, keeps first)
5. **Regrid every timestep onto the reference domain grid**

---

## Physical conversion

\[
\text{rain\_rate} = \frac{\text{accumulated\_rain (mm)}}{\Delta t \cdot 3600}
\]

Using **1 mm = 1 kg m⁻²** (liquid water equivalent).

Set the accumulation window with:
- `--accum-hours` (default 1.0)

---

## Input requirements for merging

Each input file must expose latitude/longitude coordinates so it can be regridded.
All timesteps are interpolated (or reprojected) onto the **domain grid you pass with `--domain`**, so the source grids can differ but must contain valid spatial metadata.

---

## Installation

```bash
pip install numpy xarray netCDF4
```

---

## Usage

> **Required**: pass the reference domain NetCDF (`--domain path/to/domain.nc`), built from `cdl/domain.cdl`.

### Single file
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z1200.nc \
  --out rain_time_dependent_20251221Z1200.nc \
  --domain domain.nc
```

### Multiple files using a glob
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z*.nc \
  --out rain_time_series_20251221.nc \
  --domain domain.nc
```

### Multiple explicit files
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z1200.nc wrf5_d02_20251221Z1300.nc wrf5_d02_20251221Z1400.nc \
  --out rain_time_series_20251221.nc \
  --domain domain.nc
```

### Output to a directory
If `--out` points to a directory (existing or not), the script writes
`rain_rate_merged.nc` inside it.

```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z*.nc \
  --out data/ \
  --domain domain.nc
```

### Non-hourly accumulations
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z*.nc \
  --out rain_time_series.nc \
  --domain domain.nc \
  --accum-hours 3
```

---

## CLI reference

| Option | Description |
|------|-------------|
| `--in PATH [PATH ...]` | One or more input NetCDFs. Can be globs. |
| `--out PATH` | Output merged NetCDF (file path or directory) |
| `--rain-var NAME` | Rain accumulation variable name (**only `DELTA_RAIN` is supported**) |
| `--accum-hours H` | Accumulation window in hours (default 1.0) |
| `--domain PATH` | Domain NetCDF defining the target grid (compliant with `cdl/domain.cdl`) |
| `--no-sort` | Disable sorting by time |
| `--no-dedupe` | Disable duplicate time removal |
| `--log-level LEVEL` | `DEBUG`, `INFO`, `WARNING`, `ERROR` |

---

## Output details

### Variable
`rain_rate(time, latitude, longitude)`

Attributes include:
- `standard_name = lwe_precipitation_rate`
- `units = kg m-2 s-1`
- `_FillValue` propagated from input (or default `1e37`)
- Latitude/longitude copied from the supplied domain NetCDF

### Global attributes
Includes:
- `input_files`
- `input_rain_variable`
- `accumulation_period_hours`
- `reference_domain` (absolute path to the domain NetCDF)

---

## Quick verification

```python
import xarray as xr
ds = xr.open_dataset("rain_time_series.nc")
print(ds.time.values[:5])
ds.rain_rate.isel(time=0).plot()
```

---

## Troubleshooting

- **Reference grid missing**: pass the correct `--domain` NetCDF that matches `cdl/domain.cdl`.
- **Rain var missing**: the script only supports `DELTA_RAIN`; check your inputs.
- **Duplicate times**: by default duplicates are removed; use `--no-dedupe` if needed.
