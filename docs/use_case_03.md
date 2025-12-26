# Use Case: September 23rd, 2025 – Severe rainfall event in Campania Region, Italy

This walkthrough reproduces a 12-hour LPERFECT run (03:00–15:00 UTC on 23 September 2025) blending **Civil Protection Department Weather Radar** mosaics with **WRF-5 km (d02)** and **WRF-1 km (d03)** guidance from the University of Naples Parthenope. It shows how to fetch the inputs, convert WRF accumulations to rainfall with `utils/meteouniparthenope-wrf_to_rainrate.py`, configure the model, run two 6-hour legs with restarts, and visualize the output.

## 1. Prepare local folders
```bash
mkdir -p data/radar/20250923 data/wrf/20250923/d02 data/wrf/20250923/d03 data/restart
```

## 2. Download the domain template
Use the radar-aligned domain grid provided for LPERFECT tutorials:
```bash
wget -O data/domain.nc \
  "https://data.meteo.uniparthenope.it/instruments/rdr0/tutorial/lperfectm/domain.nc"
```

## 3. Retrieve the radar GeoTIFFs (10-minute steps)
Radar images are published every 600 seconds at
```
https://data.meteo.uniparthenope.it/instruments/rdr0/YYYY/MM/DD/rdr0_d01_YYYYMMDDZhhmm_VMI.tiff
```
For the 03:00–15:00 UTC window on 2025-09-23 the timestamps span 03:00 through 14:50 (72 files). Download them in chronological order:
```bash
base_url="https://data.meteo.uniparthenope.it/instruments/rdr0/2025/09/23"
for hh in {03..14}; do
  for mm in 00 10 20 30 40 50; do
    hhmm="${hh}${mm}"
    url="${base_url}/rdr0_d01_20250923Z${hhmm}_VMI.tiff"
    wget -q --show-progress -O "data/radar/20250923/${hhmm}.tiff" "$url"
  done
done
```

## 4. Convert radar reflectivity to rain rate NetCDF
Stack the GeoTIFFs into a CF-compliant rainfall forcing aligned to `data/domain.nc`.
```bash
python utils/wr_to_rainrate.py \
  --input $(for hh in {03..14}; do for mm in 00 10 20 30 40 50; do printf "data/radar/20250923/%02d%s.tiff," "$hh" "$mm"; done; done | sed 's/,$//') \
  --output data/20250923Z0300_radar.nc \
  --time 2025-09-23T03:00:00Z \
  --dt 600 \
  --domain data/domain.nc \
  --source-name "Civil Protection Department Weather Radar" \
  --institution "Italian Department of Civil Protection" \
  --source "https://data.meteo.uniparthenope.it/instruments/rdr0"
```
- `--time` is the timestamp of the first image (03:00 UTC).
- `--dt 600` reflects the 10-minute cadence.

## 5. Download WRF guidance (hourly files)
WRF archives are hosted at:
- **Italy 5 km (d02):** `https://data.meteo.uniparthenope.it/files/wrf5/d02/archive/YYYY/MM/DD/wrf5_d02_YYYYMMDDZhh00.nc`
- **Campania 1 km (d03):** `https://data.meteo.uniparthenope.it/files/wrf5/d03/archive/YYYY/MM/DD/wrf5_d03_YYYYMMDDZhh00.nc`

For 03:00–15:00 UTC on 2025-09-23, grab the hourly accumulations (13 files per domain):
```bash
YMD=20250923
for hh in {03..15}; do
  wget -q --show-progress \
    -O "data/wrf/20250923/d02/wrf5_d02_${YMD}Z${hh}00.nc" \
    "https://data.meteo.uniparthenope.it/files/wrf5/d02/archive/2025/09/23/wrf5_d02_${YMD}Z${hh}00.nc"
  wget -q --show-progress \
    -O "data/wrf/20250923/d03/wrf5_d03_${YMD}Z${hh}00.nc" \
    "https://data.meteo.uniparthenope.it/files/wrf5/d03/archive/2025/09/23/wrf5_d03_${YMD}Z${hh}00.nc"
done
```

## 6. Convert WRF accumulations to rain rate on the domain grid
Use `utils/meteouniparthenope-wrf_to_rainrate.py` to merge and regrid the WRF files into two rainfall NetCDFs. The script reads hourly accumulated rain (`DELTA_RAIN`), converts it to rain rate (kg m-2 s-1), sorts/deduplicates time, and reprojects to `data/domain.nc`.
```bash
# 5 km Italy (d02)
python utils/meteouniparthenope-wrf_to_rainrate.py \
  --in data/wrf/20250923/d02/wrf5_d02_20250923Z*.nc \
  --out data/20250923Z0300_wrf5_d02.nc \
  --accum-hours 1.0 \
  --domain data/domain.nc \
  --log-level INFO

# 1 km Campania (d03)
python utils/meteouniparthenope-wrf_to_rainrate.py \
  --in data/wrf/20250923/d03/wrf5_d03_20250923Z*.nc \
  --out data/20250923Z0300_wrf5_d03.nc \
  --accum-hours 1.0 \
  --domain data/domain.nc \
  --log-level INFO
```

## 7. Create the simulation configuration
Save the following as `config_use_case_03.json` (or adjust your config). It blends radar with both WRF sources (weights can be tuned; here radar dominates) and runs two 6-hour legs over the 12-hour window.
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
      "name": "campania_20250923",
      "domain_nc": "data/domain.nc",
      "output": {
        "out_netcdf": "data/20250923Z0300_flood_depth.nc",
        "save_every_s": 3600,
        "rotate_every_s": 3600,
        "outflow_geojson": null,
        "Conventions": "CF-1.10",
        "title": "LPERFECT flood depth + hydrogeological risk index",
        "institution": "UniParthenope"
      },
      "restart": {
        "in": null,
        "out": "data/20250923Z0300_restart_state.nc",
        "every": 120,
        "strict_grid_check": true
      }
    }
  ],
  "model": {
    "start_time": "2025-09-23T03:00:00Z",
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
      "radar": {
        "kind": "netcdf",
        "path": "data/20250923Z0300_radar.nc",
        "var": "rain_rate",
        "time_var": "time",
        "select": "previous",
        "mode": "intensity_mmph",
        "weight": 0.6
      },
      "wrf_d02": {
        "kind": "netcdf",
        "path": "data/20250923Z0300_wrf5_d02.nc",
        "var": "rain_rate",
        "time_var": "time",
        "select": "previous",
        "mode": "intensity_kgm2s",
        "weight": 0.2
      },
      "wrf_d03": {
        "kind": "netcdf",
        "path": "data/20250923Z0300_wrf5_d03.nc",
        "var": "rain_rate",
        "time_var": "time",
        "select": "previous",
        "mode": "intensity_kgm2s",
        "weight": 0.2
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
    "out": "data/20250923Z0300_restart_state.nc",
    "every": 120,
    "strict_grid_check": true
  },
  "output": {
    "out_netcdf": "data/20250923Z0300_flood_depth.nc",
    "save_every_s": 3600,
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

## 8. Run the model with restarts (two 6-hour legs)
The configuration uses `model.T_s = 21600` (6 hours). Run twice to cover 03:00–15:00 UTC:
```bash
# Leg 1: 03:00–09:00 UTC
python main.py --config config_use_case_03.json \
  --out-nc data/20250923Z0300_flood_depth.nc \
  --restart-out data/restart/20250923Z0300_t06h.nc

# Leg 2: 09:00–15:00 UTC
python main.py --config config_use_case_03.json \
  --restart-in data/restart/20250923Z0300_t06h.nc \
  --out-nc data/20250923Z0300_flood_depth.nc
```

## 9. Visualize flood depth with the DEM as basemap
Render the first time slice to a PNG (omit `--out` for interactive mode):
```bash
python utils/output_plot.py \
  --flood data/20250923Z0300_flood_depth.nc \
  --domain data/domain.nc \
  --time-index 0 \
  --out data/20250923Z0300_flood_depth.png \
  --bbox 13.7 39.9 15.9 41.6 \
  --overlay data/boundaries/campania_municipalities.geojson \
  --overlay-label-field name \
  --title "LPERFECT flood depth – 2025-09-23 03:00 UTC"
```
The bounding box above frames Regione Campania (lon/lat order), speeding up plotting and zooming into the impacted area. Labels are placed at feature centroids when the overlay column exists; omit `--overlay-label-field` to auto-detect common name fields.

## 10. Submit a Slurm batch job (HPC)
To run the two-leg workflow on a Slurm cluster, create `run_lperfect_use_case_03.slurm`:
```bash
#!/bin/bash
#SBATCH --job-name=lperfect-campania-20250923
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --partition=normal
#SBATCH --output=logs/lperfect_%j.out
#SBATCH --error=logs/lperfect_%j.err

module load python/3.11
source /path/to/venv/bin/activate

set -euo pipefail
mkdir -p logs data/restart

CONFIG=config_use_case_03.json
OUT=data/20250923Z0300_flood_depth.nc

# Leg 1: 03–09 UTC
python main.py --config "$CONFIG" \
  --out-nc "$OUT" \
  --restart-out data/restart/20250923Z0300_t06h.nc

# Leg 2: 09–15 UTC
python main.py --config "$CONFIG" \
  --restart-in data/restart/20250923Z0300_t06h.nc \
  --out-nc "$OUT"
```
Submit with `sbatch run_lperfect_use_case_03.slurm`, adjusting resources, modules, and virtual environment paths for your cluster.

## 11. Recap
1. Download `data/domain.nc`.
2. Pull 03:00–14:50 UTC radar GeoTIFFs for 2025-09-23 and convert them to `data/20250923Z0300_radar.nc`.
3. Download 5 km (d02) and 1 km (d03) WRF hourly accumulations for 03:00–15:00 UTC and convert them with `utils/meteouniparthenope-wrf_to_rainrate.py` to domain-aligned rainfall NetCDFs.
4. Blend radar + WRF sources in `config_use_case_03.json` and run two 6-hour legs with restarts.
5. Plot `flood_depth` over the DEM to inspect impacted Campania areas.
