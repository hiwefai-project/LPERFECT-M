# LPERFECT – Input Data Preparation Guide

This document describes **how to prepare all input datasets** required by **LPERFECT**, in a way that is **fully consistent with the provided NcML specifications**:

- `ncml/domain.ncml`
- `ncml/rain_time_dependent.ncml`
- `ncml/rain_static.ncml`
- `ncml/rain_bundle_optional.ncml`
- `ncml/output_flood_depth.ncml`
- `ncml/restart_state.ncml`

LPERFECT **only uses NetCDF inputs** and follows **CF-1.10 conventions** to ensure interoperability, reproducibility, and long-term maintainability.

---

## 1. Overview of Required Inputs

LPERFECT requires two categories of input data:

1. **Domain data** (static, gridded)
2. **Rainfall forcing data** (static or time-dependent)

All datasets must:
- share the **same spatial grid** (same `latitude`, `longitude`, resolution, CRS),
- be stored as **NetCDF files**,
- follow the variable names and metadata described in the NcML files.

---

## 2. Domain Dataset (`domain.nc`)

**Reference specification:** `ncml/domain.ncml`

The domain dataset defines the **computational terrain and hydrological structure**.

### 2.1 Required Variables

The following 2D variables (dimensions `latitude, longitude`) are required:

| Variable | Description | Units |
|--------|------------|-------|
| `dem` | Digital Elevation Model | meters |
| `d8` | D8 flow direction | encoded |
| `cn` | SCS Curve Number | 0–100 |

Optional:

| Variable | Description |
|--------|------------|
| `channel_mask` | River/channel indicator (0/1) |

### 2.2 Coordinates

The dataset must define **horizontal coordinates**:

- `latitude(latitude)` – latitude coordinate  
- `longitude(longitude)` – longitude coordinate  

Attributes:
- `description`: `Latitude` / `Longitude`
- `units`: `degrees_north` / `degrees_east`

LPERFECT assumes:
- regular grid spacing,
- cell area derived from coordinate spacing.

### 2.3 Coordinate Reference System (CRS)

A **CF-compliant grid mapping variable** is strongly recommended.

Example:
```text
crs:grid_mapping_name = "transverse_mercator"
crs:epsg_code = "EPSG:32633"
```

Each spatial variable (`dem`, `cn`, etc.) must reference it via:
```text
grid_mapping = "crs"
```

The NcML template includes a `crs` variable with CF grid-mapping attributes. If you omit
the CRS variable, ensure your tools can still infer the grid geometry.

### 2.4 DEM Preparation

The DEM should:
- be hydrologically conditioned (pit-filled),
- use a consistent vertical datum,
- contain a `_FillValue` for no-data areas.

No-data cells define **inactive areas** where no particles are generated.

### 2.5 D8 Flow Direction

The D8 raster must be derived from the DEM using tools such as:
- GRASS GIS (`r.watershed`),
- TauDEM,
- RichDEM.

Supported encodings:
- **ESRI**: 1,2,4,8,16,32,64,128
- **Clockwise**: 0–7

The encoding must match `model.encoding` in `config.json`.

### 2.6 Curve Number (CN)

Curve Number values should be derived from:
- land cover data (e.g., CORINE),
- hydrologic soil groups.

Requirements:
- valid range: 0–100,
- floating-point or integer,
- aligned exactly with DEM grid.

---

## 3. Rainfall Forcing Datasets

LPERFECT supports **multiple rainfall sources**, each provided as a NetCDF file.

### 3.1 Time-Dependent Rainfall (`rain_time_dependent.nc`)

**Reference specification:** `ncml/rain_time_dependent.ncml`

Used for:
- radar nowcasting,
- station-based analyses,
- numerical weather prediction forecasts.

#### Required Dimensions

| Dimension | Meaning |
|---------|--------|
| `time` | time steps |
| `latitude` | grid rows |
| `longitude` | grid columns |

#### Required Variables

| Variable | Description | Units |
|--------|------------|-------|
| `time(time)` | time coordinate | hours since 1900-01-01 00:00:0.0 |
| `rain_rate(time,latitude,longitude)` | rainfall rate | mm h-1 |

#### Time Coordinate

The `time` variable must include:
- `description`: `Time`
- `long_name`: `time`
- `units`: `hours since 1900-01-01 00:00:0.0`

LPERFECT supports:
- nearest-time selection,
- step-index-based selection.

### 3.2 Static Rainfall (`rain_static.nc`)

**Reference specification:** `ncml/rain_static.ncml`

Used for:
- design storms,
- synthetic experiments,
- constant rainfall scenarios.

Required variable:

| Variable | Description | Units |
|--------|------------|-------|
| `rain_rate(latitude,longitude)` | rainfall rate | mm h-1 |

No time dimension is present.

### 3.3 Optional Multi-Source Bundle (`rain_bundle_optional.nc`)

**Reference specification:** `ncml/rain_bundle_optional.ncml`

This optional format packages multiple rainfall sources in a single file using:
- `radar_rain_rate(time,latitude,longitude)`
- `station_rain_rate(time,latitude,longitude)`
- `model_rain_rate(time,latitude,longitude)`

Each variable uses units of `mm h-1` and shares the same `time`, `latitude`, and
`longitude` coordinates as the other rainfall inputs.

---

## 4. Spatial Consistency Rules

All input datasets **must be spatially consistent**:

- same grid shape `(latitude,longitude)`,
- same `latitude` and `longitude` coordinate values,
- same CRS,
- same orientation (no flipped axes).

LPERFECT does **not** perform reprojection or resampling internally.

---

## 5. Output and Restart Datasets (Reference)

While LPERFECT writes outputs automatically, the following NcML files describe the
expected structure for post-processing and validation.

### 5.1 Flood Depth Output (`output_flood_depth.nc`)

**Reference specification:** `ncml/output_flood_depth.ncml`

Expected variables:
- `flood_depth(latitude,longitude)` in meters
- `risk_index(latitude,longitude)` dimensionless

### 5.2 Restart State (`restart_state.nc`)

**Reference specification:** `ncml/restart_state.ncml`

Expected fields include:
- `P_cum_mm(latitude,longitude)` and `Q_cum_mm(latitude,longitude)`
- Particle arrays (`particle_r`, `particle_c`, `particle_vol`, `particle_tau`)
- Scalars like `elapsed_s`, `cum_rain_vol_m3`, `cum_runoff_vol_m3`, `cum_outflow_vol_m3`

---

## 6. Quality Checks (Strongly Recommended)

Before running LPERFECT, verify:

- no NaNs in active DEM cells,
- D8 directions are valid everywhere DEM is valid,
- CN values are within [0,100],
- rainfall units are correct (`mm h-1`),
- time axis is monotonic and CF-compliant.

Tools:
- `ncdump -h`
- `xarray.open_dataset()`
- `cdo sinfo`

---

## 6. Example Workflow

1. Prepare DEM → fill sinks
2. Compute D8 flow directions
3. Derive CN raster
4. Export all layers to NetCDF
5. Validate against `domain.ncml`
6. Prepare rainfall NetCDFs
7. Validate against `rain_time_dependent.ncml`
8. Run LPERFECT

---

## 7. Common Pitfalls

- Mixing geographic and projected coordinates
- Incorrect D8 encoding
- Rainfall in mm per timestep instead of mm h⁻¹
- Missing or inconsistent CRS metadata

---

## 8. Final Notes

Following this guide ensures:
- reproducible simulations,
- physical consistency,
- full compatibility with LPERFECT and Hi-WeFAI workflows.

For advanced workflows, domain and rain preparation can be automated in HPC preprocessing pipelines.

---

## 9. Italy-tailored input data sources

This section lists **commonly used Italian / Europe-wide datasets** that work well with LPERFECT, and practical notes on preparing them.

### 9.1 DEM sources (Italy / Europe)

Typical options (choose one and be consistent end-to-end):

- **TINITALY DEM** (national-scale DEM for Italy, often used in research)
- **EU-DEM / Copernicus DEM** (Europe-wide DEM products; useful for broader domains)
- **Regional DEMs** (e.g., Campania, Lazio, etc., often provided by regional geoportals)

**Recommendation:** use a *projected CRS* for Italy (e.g., UTM zones 32N/33N depending on domain) to keep units in meters and cell areas meaningful.

### 9.2 Land cover for CN derivation

- **CORINE Land Cover (CLC 2018 / 2020)** is commonly used for Italy-wide CN workflows.
- You will map CLC classes + **Hydrologic Soil Group (HSG)** to a CN value using a CN lookup table (TR-55 style).

### 9.3 Soil / Hydrologic Soil Group (HSG)

- A practical soil input for Europe is **HYSOGs 250m**, which can be reclassified into HSG classes (A/B/C/D).

---

## 10. Command-line preprocessing recipes (GDAL / CDO / NCO)

Below are **operational command-line recipes**. They are written as “building blocks” you can compose in pipelines.
All commands assume a Bash-like shell on Linux/macOS. Replace paths, EPSG codes, and resolution as needed.

### 10.1 Step A — Choose grid, CRS, and resolution

Pick a target CRS and resolution that matches your operational needs (e.g., 10–100 m for local domains, 100–250 m for regional domains).

Example target CRS for Italy (choose the appropriate zone):
- UTM 32N: `EPSG:32632`
- UTM 33N: `EPSG:32633`

### 10.2 Step B — Prepare DEM (warp → align → fill sinks)

#### B1) Reproject/resample the raw DEM onto your target grid (GDAL)

```bash
gdalwarp -t_srs EPSG:32633 -tr 100 100 -r bilinear -of GTiff \
  -dstnodata -9999 \
  raw_dem.tif dem_utm33_100m.tif
```

#### B2) Ensure consistent extent and alignment (optional but recommended)

If you have a reference raster you want to match exactly (e.g., a model grid), use `-te` (extent) and `-tap` (align pixels):

```bash
gdalwarp -t_srs EPSG:32633 -tr 100 100 -tap -te xmin ymin xmax ymax \
  -r bilinear -of GTiff -dstnodata -9999 \
  raw_dem.tif dem_aligned.tif
```

#### B3) Fill sinks / depressions (hydrological conditioning)

Options:
- **GRASS GIS**: `r.fill.dir` or `r.hydrodem`
- **TauDEM**
- **RichDEM**

Example (RichDEM, if available):

```bash
rd_fill_depressions dem_aligned.tif dem_filled.tif
```

(If you prefer GRASS, a typical approach is to import the raster into a GRASS location matching your CRS and run `r.hydrodem`.)

### 10.3 Step C — Compute D8 from DEM

Tool choices:
- GRASS GIS `r.watershed` (supports D8 outputs)
- TauDEM `d8flowdir`
- RichDEM `rd_flowdir_d8`

Example (TauDEM):

```bash
mpiexec -n 8 d8flowdir -p d8_esri.tif -fel dem_filled.tif
```

**Important:** ensure the produced D8 encoding matches LPERFECT’s `model.encoding` (ESRI or clockwise).

### 10.4 Step D — Prepare CN from CORINE (CLC) + HSG (HYSOGs)

This is a common Italy workflow:

1) Reproject/resample **CLC** to the DEM grid (nearest neighbor!).
2) Reproject/resample **HSG** raster to the DEM grid (nearest neighbor!).
3) Reclassify (CLC,HSG) pairs into CN via a lookup table.

#### D1) Warp CLC onto the DEM grid (GDAL, preserve classes)

```bash
gdalwarp -t_srs EPSG:32633 -tr 100 100 -r near -of GTiff \
  -dstnodata 0 \
  clc_italy.tif clc_utm33_100m.tif
```

#### D2) Warp HYSOGs/HSG onto the DEM grid (GDAL, preserve classes)

```bash
gdalwarp -t_srs EPSG:32633 -tr 100 100 -r near -of GTiff \
  -dstnodata 0 \
  HYSOGs250m.tif hsg_utm33_100m.tif
```

#### D3) Build CN raster by reclassification

If you already have a Python script that maps `(CLC, HSG) → CN` (recommended for maintainability), the command-line step is simply:

```bash
python cn_from_clc_hsg.py \
  --dem dem_filled.tif \
  --clc clc_utm33_100m.tif \
  --hsg hsg_utm33_100m.tif \
  --out cn_utm33_100m.tif
```

If you want a pure-GDAL approach, you can reclassify using `gdal_calc.py`, but it becomes unwieldy for full CN tables; Python is preferred.

### 10.5 Step E — Optional channel mask generation

Options:
- derive from **flow accumulation threshold**,
- ingest an existing river network raster and rasterize it to the grid.

#### E1) Flow accumulation threshold (example approach)

Using TauDEM (example):

```bash
mpiexec -n 8 aread8 -p d8_esri.tif -ad8 flow_acc.tif
```

Then threshold to channel mask (GDAL calc):

```bash
gdal_calc.py -A flow_acc.tif --outfile=channel_mask.tif \
  --calc="(A>={{THRESH}}).astype(uint8)" --NoDataValue=0
```

### 10.6 Step F — Convert rasters to NetCDF and assemble `domain.nc`

You can use GDAL to create NetCDF from rasters, but assembling a single **CF-friendly** domain file is usually cleaner with **NCO** / **xarray**. Below are practical patterns.

#### F1) Convert each GeoTIFF to NetCDF (GDAL)

```bash
gdal_translate -of netCDF dem_filled.tif dem.nc
gdal_translate -of netCDF d8_esri.tif d8.nc
gdal_translate -of netCDF cn_utm33_100m.tif cn.nc
gdal_translate -of netCDF channel_mask.tif channel_mask.nc
```

> Depending on GDAL version, variable names may be `Band1` or similar. You will likely rename variables afterward.

#### F2) Rename variables to match `domain.ncml` (NCO)

```bash
ncrename -v Band1,dem dem.nc
ncrename -v Band1,d8  d8.nc
ncrename -v Band1,cn  cn.nc
ncrename -v Band1,channel_mask channel_mask.nc
```

#### F3) Merge variables into a single `domain.nc` (NCO)

```bash
ncks -A d8.nc domain.nc  # or start from dem.nc
ncks -A cn.nc domain.nc
ncks -A channel_mask.nc domain.nc
```

A common clean approach is:
1) `cp dem.nc domain.nc`
2) append others with `ncks -A` as above.

#### F4) Ensure dimensions are `latitude,longitude` and create consistent coords (CDO/NCO)

If your NetCDF came with different dim names, rename dimensions/vars with NCO.
(Exact commands vary by file; inspect with `ncdump -h` first.)

```bash
ncdump -h domain.nc | head -n 60
```

Then apply `ncml/domain.ncml` as a compliance checklist: make sure the variables and attributes match.

### 10.7 Step G — Prepare rainfall inputs

#### G1) Convert gridded rainfall GeoTIFF time slices to a CF time-series NetCDF

If you have rainfall frames as GeoTIFFs (one per timestep), you can stack them using CDO:

```bash
# Example: rain_000.tif rain_001.tif ... each already aligned to the domain grid
for f in rain_*.tif; do
  gdal_translate -of netCDF "$f" "${f%.tif}.nc"
done

# Merge into a single time-series (CDO). You may need to set time stamps first.
cdo mergetime rain_*.nc rain_time_dependent.nc
```

Then set correct units and rename to `rain_rate`:

```bash
ncrename -v Band1,rain_rate rain_time_dependent.nc
ncatted -a units,rain_rate,o,c,"mm h-1" rain_time_dependent.nc
```

Finally ensure `time(time)` exists and is compliant (`hours since 1900-01-01 00:00:0.0`).

#### G2) Create a static rainfall file

For constant/design rainfall on the grid:

```bash
# Create a constant field from a reference raster (here: dem) via gdal_calc
gdal_calc.py -A dem_filled.tif --outfile=rain_static.tif \
  --calc="A*0+{{CONST_MMPH}}" --NoDataValue=-9999

gdal_translate -of netCDF rain_static.tif rain_static.nc
ncrename -v Band1,rain_rate rain_static.nc
ncatted -a units,rain_rate,o,c,"mm h-1" rain_static.nc
```

Validate against `ncml/rain_static.ncml`.

---

## 11. Validation against NcML

Use `ncdump -h` and compare to the NcML specs:

- `ncml/domain.ncml`
- `ncml/rain_time_dependent.ncml`
- `ncml/rain_static.ncml`

Typical checks:

```bash
ncdump -h domain.nc | sed -n '1,120p'
ncdump -h rain_time_dependent.nc | sed -n '1,160p'
ncdump -h rain_static.nc | sed -n '1,120p'
```

Check:
- variable names match (`dem`, `d8`, `cn`, `rain_rate`),
- dimensions are correct (`latitude,longitude` and optionally `time`),
- units are correct (`m`, `mm h-1`, etc.),
- CRS / grid_mapping is present (recommended).

---

## 12. Practical “Italy” defaults

These choices work well for many Italian regional deployments:

- **CRS:** UTM 33N (`EPSG:32633`) for central/southern Italy (adjust if needed).
- **Resolution:** 100 m (regional), 25–50 m (local), 250 m (very large domains / fast ensembles).
- **CN derivation:** CORINE + HYSOGs/HSG, with TR-55 style lookup table.
- **Channel mask:** flow-accumulation threshold (calibrated to basin size and resolution).

---
