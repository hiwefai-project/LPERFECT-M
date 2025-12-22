# Getting Started with LPERFECT-M
*A step-by-step guide to run LPERFECT-M with real terrain and rainfall data*

LPERFECT-M is a Lagrangian, MPI-parallel hydrological model designed for fast runoff and flood simulation on gridded domains.  
This guide walks you through the complete workflow, from real GIS data to a working flood simulation.

---

## 1. Prerequisites

### Software requirements
- Python ≥ 3.9
- MPI implementation (OpenMPI or MPICH) for parallel runs
- GDAL (for raster preprocessing)

Python dependencies are listed in `requirements.txt`.

---

## 2. Install LPERFECT-M

```bash
git clone https://github.com/hiwefai-project/LPERFECT-M
cd LPERFECT-M

python -m venv .venv
source .venv/bin/activate

pip install -r requirements.txt
```

---

## 3. Required Inputs

LPERFECT-M requires two main inputs:

1. **Domain NetCDF**
2. **Rainfall NetCDF** (one or more)

### 3.1 Domain NetCDF

The domain NetCDF must contain the following 2D variables (dimensions `latitude, longitude`):

| Variable | Description | Units |
|--------|-------------|-------|
| `dem`  | Digital elevation model | meters |
| `d8`   | Flow direction (D8) | encoded |
| `cn`   | Curve Number | 0–100 |
| `channel_mask` *(optional)* | River network mask | 0/1 |

It must also include 1D coordinate variables:
- `latitude(latitude)` – latitude coordinate
- `longitude(longitude)` – longitude coordinate

And a CF grid mapping variable:
- `crs` with `grid_mapping_name`, `epsg_code`, `semi_major_axis`, and `inverse_flattening`
- Each spatial variable should include `grid_mapping = "crs"`

> Coordinate names can be remapped via `domain.varmap` in `config.json` if your files
> use `x/y` or other naming conventions.

---

## 4. Preparing the Domain Data

### 4.1 Choose Grid and Projection

- Use a **projected CRS** (e.g. UTM)
- Choose a consistent resolution (e.g. 50 m, 100 m)
- All rasters **must share the same grid**

---

### 4.2 Prepare DEM

Reproject and resample the DEM:

```bash
gdalwarp -t_srs EPSG:32633 -tr 100 100 -r bilinear   -dstnodata -9999 dem_raw.tif dem_utm_100m.tif
```

Recommended preprocessing:
- Sink filling
- Removal of spikes
- Consistent nodata handling

---

### 4.3 Compute D8 Flow Direction

Compute D8 flow directions from the DEM using a hydrology tool such as:
- GRASS GIS (`r.watershed`)
- WhiteboxTools
- TauDEM

LPERFECT-M supports two encodings:
- `esri` (powers of two)
- `cw0_7` (clockwise 0–7)

Set the corresponding encoding in `config.json`.

> Note: the CDL template in `/cdl/domain.cdl` encodes ESRI D8 values (1,2,4,8,16,32,64,128).

---

### 4.4 Generate Curve Number (CN) Map

CN values must be in the range **0–100**.

Typical workflow:
1. Start from land cover (e.g. CORINE)
2. Combine with Hydrologic Soil Group (HSG)
3. Apply CN lookup table
4. Regrid to the DEM grid

Values outside `[0,100]` are treated as `0`.

---

### 4.5 Optional: Channel Mask

If available, create a `channel_mask` raster:
- 1 = channel
- 0 = non-channel

This can be derived from flow accumulation thresholds.

---

### 4.6 Create the Domain NetCDF

Combine `dem`, `d8`, `cn` (and optional `channel_mask`) into a single NetCDF file.

Requirements:
- All variables on dimensions `(latitude, longitude)`
- Include `latitude(latitude)` and `longitude(longitude)` coordinates
- Consistent metadata and no misaligned grids

For a ready-made CLI workflow, see [`docs/make_domains.md`](make_domains.md).

---

## 5. Preparing Rainfall Data

### 5.1 Supported Rainfall Formats

Rainfall NetCDFs are **time-varying by default**:
- `(time, latitude, longitude)` per `cdl/rain_time_dependent.cdl`
- Static grids are only accepted if you set `rain.schema.require_time_dim = false`.

Required variable:
- `rain_rate` with units **mm h-1**

Time-varying files must define `time(time)` with units `hours since 1900-01-01 00:00:0.0`.
Rainfall files should include a `crs` grid-mapping variable and a
`grid_mapping = "crs"` attribute on `rain_rate`.

---

### 5.2 Align Rainfall to Domain Grid

Rainfall rasters **must match** the domain grid exactly.

Regrid using the DEM as reference:

```bash
gdalwarp -t_srs EPSG:32633 -tr 100 100 -r bilinear   rain_raw.tif rain_aligned.tif
```

---

### 5.3 Time Handling

Since rainfall has a time dimension:
- Include a proper `time` coordinate
- Choose selection method in the config:
  - `nearest`
  - `step`

Ensure simulation time matches rainfall timestamps.

---

## 6. Configuration

Edit `config.json`:

Key fields:
- `domain.domain_nc` → path to domain NetCDF
- `model.encoding` → `esri` or `cw0_7`
- Rainfall source paths
- Simulation start time and duration

Outputs and restarts can be controlled via CLI options.

---

## 7. Running the Model

Prepare the configuration file.
```bash
cp config.json.sample config.json 
```

### Serial execution
```bash
python main.py --config config.json --out-nc flood.nc
```

### Parallel execution (MPI)
```bash
mpirun -np 8 python main.py --config config.json --out-nc flood.nc
```

Rank 0 handles NetCDF I/O in parallel runs.

---

## 8. Outputs

The model produces a NetCDF file containing fields such as:
- `flood_depth(time, latitude, longitude)`
- `risk_index(time, latitude, longitude)`
- `time(time)` coordinate (hours since 1900-01-01 00:00:0.0)

These can be visualized using:
- Python (`xarray`, `matplotlib`)
- QGIS
- Panoply

---

## 9. Common Pitfalls

| Issue | Cause | Fix |
|-----|------|-----|
| No flooding | CN too low | Verify CN values |
| No flow | Wrong D8 encoding | Check `model.encoding` |
| Misaligned rain | Different grids | Regrid rainfall |
| Time mismatch | Wrong timestamps | Fix `time` coordinate |
| Edge artifacts | Nodata in DEM | Fill or mask properly |

---

## 10. Next Steps

- Enable restart files for long simulations
- Couple with forecast rainfall
- Integrate in larger workflows (e.g. Hi-WeFAI)
- Run ensemble simulations

---

**You are now ready to run LPERFECT-M with real-world data.**
