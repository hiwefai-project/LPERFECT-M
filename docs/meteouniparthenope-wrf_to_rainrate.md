# meteo@uniparthenope to rain rate

This script converts **one or many** WRF-derived NetCDF files (e.g. `wrf5_d02_*.nc`) into a **single merged time series** NetCDF compliant with a simplified `cdl/rain_time_dependent.cdl` schema.

It is intended for hydrological / impact models that require **rain rate** instead of accumulated precipitation.

---

## What it does

### Input (WRF-derived)
WRF post-processing products often provide:
- **Accumulated precipitation** over a time window (commonly 1 hour)
- Units: **mm**

The rain accumulation is typically in:
- `RAIN_DELTA` *(as you mentioned)* or
- `DELTA_RAIN` *(as shown in the CDL)*

### Output (`rain_time_dependent`)
The output contains:
- `rain_rate(time, latitude, longitude)` in **kg m⁻² s⁻¹**

### New: merge multiple times into one time series
You can pass multiple file paths or globs. The script will:
1. Convert each file to `rain_rate`
2. Concatenate along `time`
3. Sort by time (default)
4. Deduplicate identical times (default, keeps first)

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

All input files must share the same grid:
- same latitude values
- same longitude values

If any file differs, the script stops with:
- “Grid shape mismatch” or “Grid coordinate mismatch”

---

## Installation

```bash
pip install numpy xarray netCDF4
```

---

## Usage

### Single file
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z1200.nc \
  --out rain_time_dependent_20251221Z1200.nc
```

### Multiple files using a glob
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z*.nc \
  --out rain_time_series_20251221.nc
```

### Multiple explicit files
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z1200.nc wrf5_d02_20251221Z1300.nc wrf5_d02_20251221Z1400.nc \
  --out rain_time_series_20251221.nc
```

### Explicit rain variable name
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z*.nc \
  --out rain_time_series.nc \
  --rain-var RAIN_DELTA
```

### Non-hourly accumulations
```bash
python convert_wrf_rain_to_rain_time_dependent.py \
  --in wrf5_d02_20251221Z*.nc \
  --out rain_time_series.nc \
  --accum-hours 3
```

---

## CLI reference

| Option | Description |
|------|-------------|
| `--in PATH [PATH ...]` | One or more input NetCDFs. Can be globs. |
| `--out PATH` | Output merged NetCDF |
| `--rain-var NAME` | Rain accumulation variable name (`RAIN_DELTA` / `DELTA_RAIN`) |
| `--accum-hours H` | Accumulation window in hours (default 1.0) |
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

### Global attributes
Includes:
- `input_files`
- `input_rain_variable`
- `accumulation_period_hours`

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

- **Grid mismatch**: ensure all inputs come from the same WRF domain/grid.
- **Rain var missing**: pass `--rain-var RAIN_DELTA` or `--rain-var DELTA_RAIN`.
- **Duplicate times**: by default duplicates are removed; use `--no-dedupe` if needed.
