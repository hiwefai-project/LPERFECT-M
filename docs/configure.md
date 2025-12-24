# LPERFECT configuration guide

This guide explains **all configuration options** for the LPERFECT model, covering both
command-line arguments and the JSON configuration file. Defaults listed here come from
`lperfect/config.py`.

## Configuration sources and precedence

LPERFECT builds its runtime configuration in this order:

1. **Defaults** from `lperfect/config.py`.
2. **JSON config file** (merged on top of defaults, recursively).
3. **CLI overrides** (only a few operational settings).

That means you can start from a minimal JSON file and only override what you need. CLI
arguments always win over JSON values for the same setting.

## Command-line interface

Run `python main.py --help` to see the available flags. The CLI currently offers:

| Argument | Default | Description | Example |
| --- | --- | --- | --- |
| `--config` | `config.json` | Path to the JSON configuration file. | `python main.py --config configs/flood.json` |
| `--log-level` | `INFO` | Logging level (`DEBUG`, `INFO`, `WARNING`, `ERROR`). | `python main.py --log-level DEBUG` |
| `--restart-in` | `None` | Resume from a restart NetCDF. Overrides `restart.in`. | `python main.py --restart-in restart_state.nc` |
| `--restart-out` | `None` | Path to write restart NetCDF. Overrides `restart.out`. | `python main.py --restart-out restart_state.nc` |
| `--out-nc` | `None` | Output NetCDF path. Overrides `output.out_netcdf`. | `python main.py --out-nc flood_depth.nc` |
| `--device` | `None` | Compute device override (`cpu` or `gpu`). | `python main.py --device gpu` |

> Tip: You can combine the CLI and JSON file. For example, keep a stable `config.json`
> and vary only the output path with `--out-nc` for batch runs.

## JSON configuration file

The JSON file mirrors the default configuration structure. A complete example
(with defaults) looks like:

```json
{
  "domain": {
    "mode": "netcdf",
    "domain_nc": "domain.nc",
    "varmap": {
      "dem": "dem",
      "d8": "d8",
      "cn": "cn",
      "channel_mask": "channel_mask",
      "x": "longitude",
      "y": "latitude"
    }
  },
  "model": {
    "start_time": "2025-12-21T00:00:00Z",
    "T_s": 7200,
    "dt_s": 5,
    "encoding": "esri",
    "ia_ratio": 0.2,
    "particle_vol_m3": 0.25,
    "travel_time_s": 5,
    "travel_time_channel_s": 1,
    "outflow_sink": true,
    "log_every": 10
  },
  "rain": {
    "schema": {
      "time_var": "time",
      "lat_var": "latitude",
      "lon_var": "longitude",
      "rain_var": "rain_rate",
      "crs_var": "crs",
      "time_units": "hours since 1900-01-01 00:00:0.0",
      "rate_units": "mm h-1",
      "require_cf": true,
      "require_time_dim": true
    },
    "sources": {
      "rain": {
        "kind": "netcdf",
        "path": "rain_time_dependent.nc",
        "var": "rain_rate",
        "time_var": "time",
        "select": "nearest",
        "mode": "intensity_mmph",
        "weight": 1.0
      }
    }
  },
  "risk": {
    "enabled": true,
    "balance": 0.55,
    "p_low": 5.0,
    "p_high": 95.0
  },
  "restart": {
    "in": null,
    "out": "restart_state.nc",
    "every": 120,
    "strict_grid_check": true
  },
  "output": {
    "out_netcdf": "flood_depth.nc",
    "Conventions": "CF-1.10",
    "title": "LPERFECT flood depth + hydrogeological risk index",
    "institution": "UniParthenope"
  },
  "compute": {
    "device": "cpu"
  }
}
```

You can shorten this file by including only the fields you want to override.

---

## `domain` section

Controls where the spatial domain is loaded from and how variables are mapped.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `domain.mode` | `netcdf` | Domain input format. Currently only NetCDF is supported. | `"mode": "netcdf"` |
| `domain.domain_nc` | `domain.nc` | Path to the domain NetCDF file. | `"domain_nc": "data/domain.nc"` |
| `domain.varmap.dem` | `dem` | Variable name for DEM (elevation). | `"dem": "elevation"` |
| `domain.varmap.d8` | `d8` | Variable name for D8 flow directions. | `"d8": "flowdir"` |
| `domain.varmap.cn` | `cn` | Variable name for Curve Number grid. | `"cn": "curve_number"` |
| `domain.varmap.channel_mask` | `channel_mask` | Variable name for channel mask (optional). | `"channel_mask": "river_mask"` |
| `domain.varmap.x` | `x` | Coordinate name for the x/longitude axis. | `"x": "longitude"` |
| `domain.varmap.y` | `y` | Coordinate name for the y/latitude axis. | `"y": "latitude"` |

> Note: The provided CDL templates use `latitude`/`longitude` coordinate names. If you
> follow those templates, set `domain.varmap.x = "longitude"` and
> `domain.varmap.y = "latitude"` (or rename the coordinates in your NetCDF).

**Example:**

```json
{
  "domain": {
    "domain_nc": "inputs/domain_italy.nc",
    "varmap": {
      "dem": "elev",
      "d8": "d8",
      "cn": "cn",
      "x": "longitude",
      "y": "latitude"
    }
  }
}
```

---

## `compute` section

Select the array backend for compute-heavy operations. GPU support uses CuPy when available,
and gracefully falls back to CPU if CuPy is missing.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `compute.device` | `cpu` | Compute device (`cpu` or `gpu`). | `"device": "gpu"` |
| `compute.shared_memory.enabled` | `false` | Enable shared-memory parallelism inside each rank. | `"enabled": true` |
| `compute.shared_memory.workers` | `null` | Worker threads (defaults to CPU cores when `null`). | `"workers": 8` |
| `compute.shared_memory.min_particles_per_worker` | `5000` | Minimum local particles before parallelizing (prevents overhead on small runs). | `20000` |
| `compute.shared_memory.chunk_size` | `65536` | Particle chunk size per task (tune for memory/cache). | `131072` |

**Example:**

```json
{
  "compute": {
    "device": "gpu",
    "shared_memory": {
      "enabled": true,
      "workers": 8,
      "min_particles_per_worker": 20000,
      "chunk_size": 131072
    }
  }
}
```

---

## `model` section

Controls the core simulation parameters.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `model.start_time` | `2025-12-21T00:00:00Z` | ISO-8601 timestamp (UTC recommended). Used to select rain data when `rain.sources.*.select = "nearest"`. Set to `null` to disable time-based selection. | `"start_time": "2025-01-01T12:00:00Z"` |
| `model.T_s` | `7200` | Total simulation duration in seconds. | `"T_s": 14400` |
| `model.dt_s` | `5` | Timestep in seconds. | `"dt_s": 10` |
| `model.encoding` | `esri` | D8 encoding (`esri` or `cw0_7`). | `"encoding": "cw0_7"` |
| `model.ia_ratio` | `0.2` | Initial abstraction ratio for SCS-CN runoff. | `"ia_ratio": 0.05` |
| `model.particle_vol_m3` | `0.25` | Target particle volume for Lagrangian routing (m³). | `"particle_vol_m3": 0.1` |
| `model.travel_time_s` | `5` | Hillslope travel time per D8 hop (seconds). | `"travel_time_s": 10` |
| `model.travel_time_channel_s` | `1` | Channel travel time per D8 hop (seconds). | `"travel_time_channel_s": 2` |
| `model.outflow_sink` | `true` | Drop particles that flow out of the domain. | `"outflow_sink": false` |
| `model.log_every` | `10` | Log diagnostic output every N steps (`0` to disable). | `"log_every": 20` |

**Example:**

```json
{
  "model": {
    "start_time": "2025-03-03T06:00:00Z",
    "T_s": 21600,
    "dt_s": 10,
    "encoding": "esri",
    "ia_ratio": 0.15,
    "particle_vol_m3": 0.2,
    "travel_time_s": 8,
    "travel_time_channel_s": 2,
    "outflow_sink": true,
    "log_every": 50
  }
}
```

---

## `rain` section

Controls rainfall sources and blending. Rainfall sources are combined as a weighted sum
each timestep. Weights do not need to sum to 1, but relative values matter.

### `rain.schema`

`rain.schema` sets **CF-1.10 defaults** and validation behavior for time-dependent rain
forcing files (see `cdl/rain_time_dependent.cdl`).

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `time_var` | `time` | Name of the CF time coordinate. | `"time_var": "time"` |
| `lat_var` | `latitude` | Name of the latitude coordinate. | `"lat_var": "latitude"` |
| `lon_var` | `longitude` | Name of the longitude coordinate. | `"lon_var": "longitude"` |
| `rain_var` | `rain_rate` | Name of the rainfall rate variable. | `"rain_var": "rain_rate"` |
| `crs_var` | `crs` | Name of the CF grid mapping variable. | `"crs_var": "crs"` |
| `time_units` | `hours since 1900-01-01 00:00:0.0` | CF time units string. | `"time_units": "hours since 1900-01-01 00:00:0.0"` |
| `rate_units` | `mm h-1` | Expected rain rate units. | `"rate_units": "mm h-1"` |
| `require_cf` | `true` | Enforce CF metadata checks for rain inputs. | `"require_cf": true` |
| `require_time_dim` | `true` | Require a time dimension on rain_rate. | `"require_time_dim": true` |

### `rain.sources`

`rain.sources` is a mapping of **source name → configuration**. Each source has:

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `kind` | `netcdf` | Source type: `netcdf` or `scalar`. | `"kind": "scalar"` |
| `weight` | `0.0` | Weight used in the blending sum. | `"weight": 0.3` |
| `mode` | `intensity_mmph` | Data units: `intensity_mmph` (mm/hour). CDL-compliant rain files use `rain_rate` with units `mm h-1`. | `"mode": "intensity_mmph"` |
| `path` | `null` | NetCDF path (required when `kind = netcdf`). | `"path": "radar.nc"` |
| `var` | `rain_rate` | NetCDF variable name (required when `kind = netcdf`). | `"var": "rain_rate"` |
| `time_var` | `time` | Name of the time coordinate. Used when data is 3D. | `"time_var": "valid_time"` |
| `lat_var` | `latitude` | Latitude coordinate name. | `"lat_var": "lat"` |
| `lon_var` | `longitude` | Longitude coordinate name. | `"lon_var": "lon"` |
| `crs_var` | `crs` | Grid mapping variable name. | `"crs_var": "crs"` |
| `time_units` | `hours since 1900-01-01 00:00:0.0` | CF time units. | `"time_units": "hours since 1900-01-01 00:00:0.0"` |
| `rate_units` | `mm h-1` | Rain rate units. | `"rate_units": "mm h-1"` |
| `require_cf` | `true` | Enforce CF validation for this source. | `"require_cf": true` |
| `require_time_dim` | `true` | Require time dimension for this source. | `"require_time_dim": true` |
| `select` | `nearest` | How to choose a time slice: `nearest` (by timestamp) or `step` (by index). | `"select": "step"` |
| `value` | `null` | Scalar intensity/depth (required when `kind = scalar`). | `"value": 2.5` |

**NetCDF source example:**

```json
{
  "rain": {
    "sources": {
      "rain": {
        "kind": "netcdf",
        "path": "rain_time_dependent.nc",
        "var": "rain_rate",
        "time_var": "time",
        "select": "nearest",
        "mode": "intensity_mmph",
        "weight": 1.0
      }
    }
  }
}
```

**Scalar source example (uniform rainfall):**

```json
{
  "rain": {
    "sources": {
      "uniform": {
        "kind": "scalar",
        "value": 1.2,
        "mode": "intensity_mmph",
        "weight": 1.0
      }
    }
  }
}
```

**Notes:**

- A 2D variable is treated as time-invariant; a 3D variable must have a time dimension.
- `select = "nearest"` requires `model.start_time` (otherwise time selection falls back to step index).
- All rain fields must match the domain grid shape.
- With `require_cf = true`, inputs must comply with `cdl/rain_time_dependent.cdl`.

---

## `risk` section

Controls the hydrogeological risk index computation. The risk index combines normalized
runoff and flow accumulation.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `risk.enabled` | `true` | Whether to compute `risk_index` in the output. | `"enabled": false` |
| `risk.balance` | `0.55` | Weighting between runoff (1.0) and flow accumulation (0.0). Clipped to `[0, 1]`. | `"balance": 0.7` |
| `risk.p_low` | `5.0` | Lower percentile for robust normalization. | `"p_low": 2.0` |
| `risk.p_high` | `95.0` | Upper percentile for robust normalization. | `"p_high": 98.0` |

**Example:**

```json
{
  "risk": {
    "enabled": true,
    "balance": 0.6,
    "p_low": 5.0,
    "p_high": 95.0
  }
}
```

---

## `restart` section

Controls restart input/output and checkpoint cadence.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `restart.in` | `null` | Restart NetCDF to resume from. | `"in": "restart_state.nc"` |
| `restart.out` | `restart_state.nc` | Restart NetCDF output path. | `"out": "outputs/restart.nc"` |
| `restart.every` | `120` | Write restart every N steps (`0` disables periodic writes). | `"every": 60` |
| `restart.strict_grid_check` | `true` | Validate that restart grids match the domain in MPI mode. | `"strict_grid_check": false` |

**Example:**

```json
{
  "restart": {
    "in": null,
    "out": "restart_state.nc",
    "every": 100,
    "strict_grid_check": true
  }
}
```

---

## `output` section

Controls final NetCDF output metadata and paths.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `output.out_netcdf` | `flood_depth.nc` | Final output NetCDF path. | `"out_netcdf": "outputs/flood_depth.nc"` |
| `output.Conventions` | `CF-1.10` | CF metadata convention string. | `"Conventions": "CF-1.8"` |
| `output.title` | `LPERFECT flood depth + hydrogeological risk index` | Global title attribute. | `"title": "LPERFECT flood run"` |
| `output.institution` | `UniParthenope` | Global institution attribute. | `"institution": "My Lab"` |

**Example:**

```json
{
  "output": {
    "out_netcdf": "flood_depth.nc",
    "Conventions": "CF-1.10",
    "title": "LPERFECT flood depth",
    "institution": "UniParthenope"
  }
}
```

---

## Minimal configuration example

A minimal JSON file that only changes rain sources and output path could be:

```json
{
  "rain": {
    "sources": {
      "rain": {
        "kind": "netcdf",
        "path": "rain_time_dependent.nc",
        "var": "rain_rate",
        "mode": "intensity_mmph",
        "weight": 1.0
      }
    }
  },
  "output": {
    "out_netcdf": "results/flood_depth.nc"
  }
}
```

This will still inherit all other defaults (domain location, timing, restart settings, etc.)
from `lperfect/config.py`.
