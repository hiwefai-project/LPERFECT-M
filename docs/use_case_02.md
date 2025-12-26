# Use Case: March 14th, 2025 – Severe rainfall event in Italy

This guide walks through a 24-hour LPERFECT run (00:00–24:00 UTC on 14 March 2025) using the **Civil Protection Department Weather Radar** mosaics. It covers data download, rainfall preprocessing with `utils/wr_to_rainrate.py`, model configuration, execution, and a quick visualization of the flood depth output using the domain DEM as a basemap.

## 1. Prepare local folders
```bash
mkdir -p data/radar/20250314
```

## 2. Download the domain template
Fetch the domain grid that matches the radar mosaic projection and resolution:
```bash
wget -O data/domain.nc \
  "https://data.meteo.uniparthenope.it/instruments/rdr0/tutorial/lperfectm/domain.nc"
```

## 3. Retrieve the radar GeoTIFFs (10-minute steps)
Radar images are published every 600 seconds at
```
https://data.meteo.uniparthenope.it/instruments/rdr0/YYYY/MM/DD/rdr0_d01_YYYYMMDDZhhmm_VMI.tiff
```
For the 00:00–24:00 UTC window on 2025-03-14 the timestamps span 00:00 through 23:50 (144 files). Download them in chronological order:
```bash
base_url="https://data.meteo.uniparthenope.it/instruments/rdr0/2025/03/14"
for hh in {00..23}; do
  for mm in 00 10 20 30 40 50; do
    hhmm="${hh}${mm}"
    url="${base_url}/rdr0_d01_20250314Z${hhmm}_VMI.tiff"
    wget -q --show-progress -O "data/radar/20250314/${hhmm}.tiff" "$url"
  done
done
```

## 4. Convert radar reflectivity to rain rate NetCDF
Use the provided utility to stack the 10-minute GeoTIFFs into a CF-compliant rainfall forcing aligned to `data/domain.nc`.
```bash
python utils/wr_to_rainrate.py \
  --input $(for hh in {00..23}; do for mm in 00 10 20 30 40 50; do printf "data/radar/20250314/%02d%s.tiff," "$hh" "$mm"; done; done | sed 's/,$//') \
  --output data/20250314Z0000_radar.nc \
  --time 2025-03-14T00:00:00Z \
  --dt 600 \
  --domain data/domain.nc \
  --source-name "Civil Protection Department Weather Radar" \
  --institution "Italian Department of Civil Protection" \
  --source "https://data.meteo.uniparthenope.it/instruments/rdr0"
```
- `--time` is the timestamp of the first image (00:00 UTC).
- `--dt 600` reflects the 10-minute cadence.
- The script reads reflectivity (dBZ), converts it to rain rate (mm/h) with the default Z–R relationship (`Z = 200 * R^1.6`), reprojects/interpolates to the domain grid, and writes `rain_rate(time, latitude, longitude)` with CF metadata.

## 5. Create the simulation configuration
Save the following as `config_use_case_02.json` (or update your config accordingly). It mirrors
the structure of `config.json.sample`, using the `domains` array even for a single domain and
retaining the output/restart fields exposed in the sample file:
```json
{
  "domain": {
    "mode": "netcdf",
    "varmap": {
      "dem": "dem",
      "d8": "d8",
      "cn": "cn",
      "channel_mask": "channel_mask",
      "x": "longitude",
      "y": "latitude"
    }
  },
  "domains": [
    {
      "name": "italy_20250314",
      "domain_nc": "data/domain.nc",
      "output": {
        "out_netcdf": "data/20250314Z0000_flood_depth.nc",
        "save_every_s": 0,
        "rotate_every_s": 3600,
        "outflow_geojson": null,
        "Conventions": "CF-1.10",
        "title": "LPERFECT flood depth + hydrogeological risk index",
        "institution": "UniParthenope"
      },
      "restart": {
        "in": null,
        "out": "data/20250314Z0000_restart_state.nc",
        "every": 120,
        "strict_grid_check": true
      }
    }
  ],
  "model": {
    "start_time": "2025-03-14T00:00:00Z",
    "T_s": 21600,
    "dt_s": 5,
    "encoding": "esri",
    "ia_ratio": 0.2,
    "particle_vol_m3": 0.25,
    "travel_time_s": 5,
    "travel_time_channel_s": 1,
    "travel_time_mode": "auto",
    "travel_time_auto": {
      "hillslope_velocity_ms": 0.5,
      "channel_velocity_ms": 1.5,
      "min_s": 0.25,
      "max_s": 3600.0
    },
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
        "path": "data/20250314Z0000_radar.nc",
        "var": "rain_rate",
        "time_var": "time",
        "select": "previous",
        "mode": "intensity_mmph",
        "weight": 1.0
      }
    }
  },
  "risk": {
    "enabled": true,
    "balance": 0.15,
    "p_low": 5.0,
    "p_high": 95.0
  },
  "restart": {
    "in": null,
    "out": "data/20250314Z0000_restart_state.nc",
    "every": 120,
    "strict_grid_check": true
  },
  "output": {
    "out_netcdf": "data/20250314Z0000_flood_depth.nc",
    "save_every_s": 0,
    "rotate_every_s": 3600,
    "outflow_geojson": null,
    "Conventions": "CF-1.10",
    "title": "LPERFECT flood depth + hydrogeological risk index",
    "institution": "UniParthenope"
  },
  "compute": {
    "device": "cpu",
    "shared_memory": {
      "enabled": true,
      "workers": 8,
      "min_particles_per_worker": 20000
    }
  }
}
```

## 6. Run the model with planned restarts (6-hour legs)
The configuration above sets `model.T_s` to 21,600 seconds (6 hours) so each run finishes after one leg.
Use `--restart-out` to capture the end-of-leg state and `--restart-in` to resume from it. Repeat until you
complete the 24-hour window (four legs):
```bash
# Leg 1: 00:00–06:00 UTC
python main.py --config config_use_case_02.json \
  --out-nc data/20250314Z0000_flood_depth.nc \
  --restart-out data/restart/20250314Z0000_t06h.nc

# Leg 2: 06:00–12:00 UTC
python main.py --config config_use_case_02.json \
  --restart-in data/restart/20250314Z0000_t06h.nc \
  --restart-out data/restart/20250314Z0000_t12h.nc \
  --out-nc data/20250314Z0000_flood_depth.nc

# Leg 3: 12:00–18:00 UTC
python main.py --config config_use_case_02.json \
  --restart-in data/restart/20250314Z0000_t12h.nc \
  --restart-out data/restart/20250314Z0000_t18h.nc \
  --out-nc data/20250314Z0000_flood_depth.nc

# Leg 4: 18:00–24:00 UTC
python main.py --config config_use_case_02.json \
  --restart-in data/restart/20250314Z0000_t18h.nc \
  --out-nc data/20250314Z0000_flood_depth.nc
```
Outputs per leg:
- `data/20250314Z0000_flood_depth.nc` (updated after each leg with `flood_depth` and, if enabled, `risk_index`).
- `data/restart/20250314Z0000_tXXh.nc` snapshots written at the end of each 6-hour leg, plus the
  intra-leg restart frequency defined by `restart.every`.

## 7. Visualize flood depth with the DEM as basemap
Use the bundled plotting utility to overlay flood depth on a DEM hillshade. The example below renders the first time step to a PNG; omit `--out` to open an interactive window.
```bash
python utils/output_plot.py \
  --flood data/20250314Z0000_flood_depth.nc \
  --domain data/domain.nc \
  --time-index 0 \
  --out data/20250314Z0000_flood_depth.png \
  --title "LPERFECT flood depth – 2025-03-14 00:00 UTC"
```
- Optional refinements:
  - Crop the map to a region (e.g., Campania) with `--bbox 13.7 39.9 15.9 41.6`.
  - Overlay and label boundaries when `geopandas` is installed, e.g. `--overlay data/boundaries/regions.geojson --overlay-label-field name`.

## 8. Submit a Slurm batch job (HPC)
To run the 6-hour-leg workflow on an HPC cluster managed by Slurm, create `run_lperfect_use_case_02.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=lperfect-italy-20250314
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=26:00:00
#SBATCH --partition=normal
#SBATCH --output=logs/lperfect_%j.out
#SBATCH --error=logs/lperfect_%j.err

module load python/3.11
source /path/to/venv/bin/activate

set -euo pipefail
mkdir -p logs data/restart

CONFIG=config_use_case_02.json
OUT=data/20250314Z0000_flood_depth.nc

declare -a LEGS=(06 12 18 24)
RESTART_IN=""

for LEG in "${LEGS[@]}"; do
  if [[ -z "$RESTART_IN" ]]; then
    python main.py --config "$CONFIG" \
      --out-nc "$OUT" \
      --restart-out "data/restart/20250314Z0000_t${LEG}h.nc"
  else
    python main.py --config "$CONFIG" \
      --restart-in "$RESTART_IN" \
      --restart-out "data/restart/20250314Z0000_t${LEG}h.nc" \
      --out-nc "$OUT"
  fi
  RESTART_IN="data/restart/20250314Z0000_t${LEG}h.nc"
done
```
Submit with `sbatch run_lperfect_use_case_02.slurm`. Adjust `#SBATCH` resources, module loads, and virtual environment paths to match your HPC environment.

## 9. Recap
1. Download `data/domain.nc`.
2. Pull 00:00–23:50 UTC radar GeoTIFFs for 2025-03-14 from the Civil Protection Department Weather Radar mosaics.
3. Convert them to `data/20250314Z0000_radar.nc` with `utils/wr_to_rainrate.py`.
4. Run the 24-hour window as four 6-hour legs with planned restarts.
5. Plot `flood_depth` over the DEM to inspect impacted areas.
