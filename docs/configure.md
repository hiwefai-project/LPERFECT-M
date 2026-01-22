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
| `--travel-time-mode` | `None` | Override `model.travel_time_mode` (`fixed` or `auto`). | `python main.py --travel-time-mode auto` |
| `--travel-time-hill-vel` | `None` | Hillslope velocity (m/s) used when `--travel-time-mode auto`. | `python main.py --travel-time-mode auto --travel-time-hill-vel 0.8` |
| `--travel-time-channel-vel` | `None` | Channel velocity (m/s) used when `--travel-time-mode auto`. | `python main.py --travel-time-mode auto --travel-time-channel-vel 2.0` |
| `--travel-time-min` | `None` | Minimum hop time (s) when `--travel-time-mode auto`. | `python main.py --travel-time-mode auto --travel-time-min 0.5` |
| `--travel-time-max` | `None` | Maximum hop time (s) when `--travel-time-mode auto`. | `python main.py --travel-time-mode auto --travel-time-max 900` |
| `--outflow-geojson` | `None` | GeoJSON path for logging sea/lake outflow hit points. Overrides `output.outflow_geojson`. | `python main.py --outflow-geojson outputs/outflow.geojson` |
| `--runoff-only-risk` | `False` | Disable Lagrangian transport and compute the risk index from runoff only. Overrides `model.runoff_only_risk`. | `python main.py --runoff-only-risk` |
| `--parallel-metrics` | `False` | Enable parallelization metrics collection. Overrides `metrics.parallelization.enabled`. | `python main.py --parallel-metrics` |
| `--parallel-metrics-output` | `None` | Path to write GPT-friendly metrics JSON. Overrides `metrics.parallelization.output`. | `python main.py --parallel-metrics --parallel-metrics-output runs/metrics.json` |
| `--ai-metrics` | `False` | Enable GPT-ready hydrology + compute metrics. Overrides `metrics.assistant.enabled`. | `python main.py --ai-metrics` |
| `--ai-metrics-output` | `None` | Path to write hydrology + compute metrics JSON. Overrides `metrics.assistant.output`. | `python main.py --ai-metrics --ai-metrics-output runs/ai_metrics.json` |

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
    "travel_time_mode": "fixed",
    "travel_time_auto": {
      "hillslope_velocity_ms": 0.5,
      "channel_velocity_ms": 1.5,
      "min_s": 0.25,
      "max_s": 3600.0
    },
    "runoff_only_risk": false,
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
    "save_every_s": 0,
    "rotate_every_s": 0,
    "outflow_geojson": null,
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

## `domain` / `domains` sections

Controls where the spatial domain is loaded from and how variables are mapped. You can:

1. Provide a single `domain` object (backwards compatible).
2. Provide a `domains` array to run multiple nested grids with explicit parent/child relationships (e.g., national 90 m, regional 30 m, city 10 m). Each domain keeps the same hierarchical / heterogeneous parallelization scheme already used by LPERFECT, and coarse domains must feed and receive particle fluxes from their children at every time step to honor the two-way nesting assumption.

Common keys for each domain object:

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `mode` | `netcdf` | Domain input format. Currently only NetCDF is supported. | `"mode": "netcdf"` |
| `domain_nc` | `domain.nc` | Path to the domain NetCDF file. | `"domain_nc": "data/domain_30m.nc"` |
| `varmap.dem` | `dem` | Variable name for DEM (elevation). | `"dem": "elevation"` |
| `varmap.d8` | `d8` | Variable name for D8 flow directions. | `"d8": "flowdir"` |
| `varmap.cn` | `cn` | Variable name for Curve Number grid. | `"cn": "curve_number"` |
| `varmap.channel_mask` | `channel_mask` | Variable name for channel mask (optional). | `"channel_mask": "river_mask"` |
| `varmap.x` | `longitude` | Coordinate name for the x/longitude axis. | `"x": "longitude"` |
| `varmap.y` | `latitude` | Coordinate name for the y/latitude axis. | `"y": "latitude"` |
| `name` | `domain_1`, … | Human-readable label used for logging and automatic file suffixes. | `"name": "city_10m"` |
| `output.out_netcdf` | inherited | Output NetCDF path for this domain; if omitted in multi-domain runs it is auto-suffixed by the domain name. | `"out_netcdf": "flood_depth_30m.nc"` |
| `restart.out` | inherited | Restart file to write for this domain; auto-suffixed when running multiple domains unless explicitly set. | `"out": "restart_10m.nc"` |
| `restart.in` | inherited | Restart file to read for this domain; auto-suffixed when running multiple domains unless explicitly set. | `"in": "restart_10m.nc"` |

> Note: The provided CDL templates use `latitude`/`longitude` coordinate names. If you
> follow those templates, set `varmap.x = "longitude"` and
> `varmap.y = "latitude"` (or rename the coordinates in your NetCDF).

### `metrics.parallelization`

Controls optional metrics that capture wall-clock timings, throughput, and migration ratios for MPI + threading + GPU runs. These metrics are formatted to be GPT-friendly for downstream optimization workflows.

| Key | Default | Description |
| --- | --- | --- |
| `enabled` | `false` | Turn on metrics collection. |
| `output` | `null` | Path to write the metrics JSON. When `null`, the JSON is logged only. |
| `max_samples` | `256` | Maximum per-step samples to retain; larger runs are down-sampled evenly. |
| `format` | `detailed` | JSON verbosity: `detailed` pretty-prints, `compact` minifies for smaller files. |

Example:

```json
{
  "metrics": {
    "parallelization": {
      "enabled": true,
      "output": "runs/parallel_metrics.json",
      "max_samples": 128
    }
  }
}
```

### `metrics.assistant`

Produces a compact GPT-friendly JSON report with computational settings and hydrological outcomes so you can ask an AI to summarize run quality, diagnose issues, or propose tuning ideas. See `docs/ai.md` for end-to-end usage examples.

| Key | Default | Description |
| --- | --- | --- |
| `enabled` | `false` | Turn on AI-assistant metrics. |
| `output` | `null` | Path to write the metrics JSON. When `null`, the JSON is logged only. |
| `format` | `detailed` | JSON verbosity: `detailed` pretty-prints, `compact` minifies for smaller files. |

Example:

```json
{
  "metrics": {
    "assistant": {
      "enabled": true,
      "output": "runs/hydrology_compute_metrics.json"
    }
  }
}
```

**Single-domain example:**

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

**Nested domains example (90 m → 30 m → 10 m):**

```json
{
  "domain": { "varmap": { "x": "longitude", "y": "latitude" } },
  "domains": [
    {
      "name": "national_90m",
      "domain_nc": "domain_90m.nc",
      "output": { "out_netcdf": "flood_depth_90m.nc" },
      "restart": { "out": "restart_90m.nc" },
      "parent": "root"
    },
    {
      "name": "regional_30m",
      "domain_nc": "domain_30m.nc",
      "output": { "out_netcdf": "flood_depth_30m.nc" },
      "restart": { "out": "restart_30m.nc" },
      "parent": "national_90m"
    },
    {
      "name": "city_10m",
      "domain_nc": "domain_10m.nc",
      "output": { "out_netcdf": "flood_depth_10m.nc" },
      "restart": { "out": "restart_10m.nc" },
      "parent": "regional_30m"
    }
  ]
}
```

Domains are validated so each `parent` is either `"root"` or the name of another
domain in the list, mirroring the hierarchy used by `utils/make_domain.py`.
LPERFECT executes domains in coarse-to-fine order to support the two-way nesting
assumption (fine grids exchange inflow/outflow with their parents at each step).

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
| `compute.mpi.enabled` | `null` | Enable/disable MPI explicitly (`true`, `false`, or `null` for auto-detect when launched under MPI). | `"enabled": true` |
| `compute.mpi.decomposition` | `auto` | MPI row-slab decomposition strategy (`auto` or `balanced`). | `"decomposition": "balanced"` |
| `compute.mpi.min_rows_per_rank` | `1` | Minimum grid rows per rank (prevents tiny slabs). | `"min_rows_per_rank": 8` |
| `compute.mpi.migration_mode` | `agg_nonblocking` | Particle migration algorithm (`agg_nonblocking` or `legacy`). | `"migration_mode": "agg_nonblocking"` |
| `compute.mpi.timing_every_steps` | `50` | Emit MPI timing summaries every N steps (0 disables). | `"timing_every_steps": 25` |
| `compute.mpi.overlap_migration` | `true` | Overlap migration communication with local work when possible. | `"overlap_migration": false` |
| `compute.mpi.balance.every_steps` | `0` | Force MPI rebalancing every N steps (0 disables). | `"every_steps": 250` |
| `compute.mpi.balance.every_sim_s` | `0` | Force MPI rebalancing every N simulated seconds (0 disables). | `"every_sim_s": 1800` |
| `compute.mpi.balance.auto` | `false` | Enable automatic load-balance when particle counts drift. | `"auto": true` |
| `compute.mpi.balance.imbalance_threshold` | `2.0` | Trigger auto-balance when max/min particle ratio exceeds this value. | `"imbalance_threshold": 1.5` |

### MPI + shared-memory (hybrid) configuration

To run LPERFECT-M with **MPI across ranks** and **shared-memory threads inside each rank**,
configure both the `compute.mpi` and `compute.shared_memory` blocks. MPI distributes the
domain by row slabs, while shared-memory workers parallelize particle operations within
each rank. This hybrid mode is ideal for multi-core nodes where you want fewer MPI ranks
and more per-rank threads.

**Checklist**

1. Install `mpi4py` and ensure an MPI runtime (OpenMPI/MPICH) is available.
2. Launch with `mpirun` or `mpiexec` using the desired number of ranks.
3. Set `compute.shared_memory.enabled = true` and pick a reasonable `workers` count **per rank**.
4. Optional: tune `compute.mpi.decomposition` and `compute.mpi.balance` for load balance.

**Hybrid MPI + threads example (CPU-only):**

```json
{
  "compute": {
    "device": "cpu",
    "mpi": {
      "enabled": true,
      "decomposition": "balanced",
      "min_rows_per_rank": 8,
      "balance": {
        "every_steps": 0,
        "every_sim_s": 0,
        "auto": true,
        "imbalance_threshold": 1.5
      }
    },
    "shared_memory": {
      "enabled": true,
      "workers": 6,
      "min_particles_per_worker": 20000,
      "chunk_size": 131072
    }
  }
}
```

**Launch command (4 ranks × 6 threads each):**

```bash
mpirun -np 4 python main.py --config config.json
```

**Notes**

- `compute.mpi.enabled = null` (default) auto-enables MPI only when launched under an MPI
  launcher with more than one rank; set it to `false` to force serial behavior even when
  running under `mpirun`.
- `compute.shared_memory.workers` is **per rank**; ensure `ranks × workers` fits the
  physical cores on each node to avoid oversubscription.
- For small runs, keep `min_particles_per_worker` higher to avoid thread overhead.

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
| `model.travel_time_mode` | `fixed` | Travel time calculation mode: `fixed` uses the scalar values above, `auto` derives hop times from cell area and user-provided velocities. | `"travel_time_mode": "auto"` |
| `model.travel_time_auto.hillslope_velocity_ms` | `0.5` | Mean hillslope velocity (m/s) used when `travel_time_mode = "auto"`. | `1.0` |
| `model.travel_time_auto.channel_velocity_ms` | `1.5` | Mean channel velocity (m/s) used when `travel_time_mode = "auto"`. | `2.5` |
| `model.travel_time_auto.min_s` | `0.25` | Lower clamp for automatically derived travel times (seconds). | `0.5` |
| `model.travel_time_auto.max_s` | `3600.0` | Upper clamp for automatically derived travel times (seconds). | `900.0` |
| `model.runoff_only_risk` | `false` | Skip particle transport and compute the risk index directly from cumulative runoff. | `"runoff_only_risk": true` |
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
    "travel_time_mode": "auto",
    "travel_time_auto": {
      "hillslope_velocity_ms": 0.75,
      "channel_velocity_ms": 2.0,
      "min_s": 0.5,
      "max_s": 1800.0
    },
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
| `select` | `previous` | How to choose a time slice: `previous`/`floor` (use the latest time ≤ simulation time), `nearest` (by timestamp midpoint), or `step` (by index). | `"select": "nearest"` |
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
              "select": "previous",
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
- `select = "previous"` (default) holds each rain rate constant until the next timestamp, which avoids mid-interval switches when your model `dt_s` is smaller than the rain timestep.
- `select = "nearest"` requires `model.start_time` (otherwise time selection falls back to step index) and switches at the midpoint between rain timestamps.
- All rain fields must match the domain grid shape.
- With `require_cf = true`, inputs must comply with `cdl/rain_time_dependent.cdl`.

---

## `risk` section

Controls the hydrogeological risk index computation. The risk index combines normalized
runoff and flow accumulation.

| Key | Default | Description | Example |
| --- | --- | --- | --- |
| `risk.enabled` | `true` | Whether to compute `risk_index` in the output. | `"enabled": false` |
| `risk.balance` | `0.55` | Weighting between normalized runoff (1.0) and the multiplicative coupling of runoff and flow accumulation (0.0). Clipped to `[0, 1]`. | `"balance": 0.7` |
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
| `output.save_every_s` | `0` | Append a new time slice to the same NetCDF every N simulated seconds (`0` disables periodic writes). | `"save_every_s": 3600` |
| `output.rotate_every_s` | `0` | Write a brand-new NetCDF every N simulated seconds, using the `out_netcdf` basename plus `_0000.nc`, `_0001.nc`, ... | `"rotate_every_s": 1800` |
| `output.outflow_geojson` | `null` | When set, write a GeoJSON listing cells where particles exit to sea/lakes, with particle counts per save interval. | `"outflow_geojson": "outputs/outflow_hits.geojson"` |
| `output.variables` | `all` | Output variables to write (`all` or list: `flood_depth`, `risk_index`, `inundation_mask`, `flood_depth_max`, `inundation_mask_max`). | `"variables": ["flood_depth", "risk_index"]` |
| `output.Conventions` | `CF-1.10` | CF metadata convention string. | `"Conventions": "CF-1.8"` |
| `output.title` | `LPERFECT flood depth + hydrogeological risk index` | Global title attribute. | `"title": "LPERFECT flood run"` |
| `output.institution` | `UniParthenope` | Global institution attribute. | `"institution": "My Lab"` |
| `output.inundation_threshold_m` | `0.01` | Depth (m) used to derive `inundation_mask` and `inundation_mask_max`. | `"inundation_threshold_m": 0.05` |
| `output.fill_value` | `-9999.0` | `_FillValue` applied to float outputs (`flood_depth`, `risk_index`, `flood_depth_max`). | `"fill_value": -32767.0` |

> Configure **either** `save_every_s` **or** `rotate_every_s`. If both are set, rotation takes precedence and the model logs a warning. The final state is always written even if it does not land exactly on the requested cadence.
> At the end of each run, a simulation quality report is logged automatically (no configuration needed), covering mass balance and hydrological checks.

**Example:**

```json
{
  "output": {
    "out_netcdf": "flood_depth.nc",
    "save_every_s": 3600,
    "outflow_geojson": "outputs/outflow_hits.geojson",
    "variables": [
      "flood_depth",
      "risk_index",
      "inundation_mask",
      "flood_depth_max",
      "inundation_mask_max"
    ],
    "Conventions": "CF-1.10",
    "title": "LPERFECT flood depth",
    "institution": "UniParthenope"
  }
}
```

When `outflow_geojson` is set (or `--outflow-geojson` is provided), the model emits a GeoJSON FeatureCollection where each point marks a cell whose advection path exits to the sea or a lake. The `particles_per_interval` property lists particle counts per save interval, and `total_particles` aggregates the full run.

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
