# `utils/output_plot.py` — Plot LPERFECT flood depth over DEM hillshade

This script creates a map figure by combining:
- **Hillshade** from a DEM (domain NetCDF), and
- **Flood depth overlay** from LPERFECT output NetCDF.

It includes the “do it all” improvements:
- FillValue/NaN handling
- Threshold masking
- Percentile-based `vmax` for better contrast
- Optional log scaling
- Optional regridding (grid alignment)
- Batch frame generation for all time steps
- Optional GeoJSON/Shapefile overlay (if installed)
- Pedagogical line-by-line comments in the code
- CLI + logging for reproducible runs

---

## Install

Required:
```bash
pip install matplotlib numpy xarray netCDF4
```

Optional overlays:
```bash
pip install geopandas shapely
```

---

## Basic usage

Interactive display:
```bash
python utils/output_plot.py --flood data/flood_depth.nc --domain data/domain.nc --time-index 0
```

Save PNG:
```bash
python utils/output_plot.py --flood data/flood_depth.nc --domain data/domain.nc --time-index 0 --out outputs/flood.png
```

---

## Better-looking maps (recommended)

Mask shallow flooding (e.g., 5 cm) and use percentile contrast:
```bash
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --threshold 0.05 \
  --vmax-percentile 99.5 \
  --out outputs/flood_thr5cm.png
```

Log scale for wide depth ranges:
```bash
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --log-scale \
  --out outputs/flood_log.png
```

---

## Grid alignment (scientific correctness)

Default is `--regrid flood_to_dem` (flood interpolated onto DEM grid), which is usually safest for visual alignment.

```bash
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --regrid flood_to_dem \
  --out outputs/flood_rg.png
```

Options:
- `none`
- `flood_to_dem` (default)
- `dem_to_flood`

---

## Batch mode (one PNG per time step)

```bash
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --all-times \
  --out-dir outputs/frames \
  --threshold 0.05
```

Turn frames into a video:
```bash
ffmpeg -framerate 5 -i outputs/frames/flood_depth_t%03d.png -pix_fmt yuv420p outputs/flood.mp4
```

---

## Overlay vector data (boundaries/roads/basins)

```bash
pip install geopandas shapely
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --overlay data/boundaries.geojson \
  --out outputs/flood_with_bounds.png
```

---

## Key CLI options

- Inputs: `--flood`, `--domain`
- Time: `--time-index`, `--all-times --out-dir`
- Variables: `--flood-var`, `--dem-var`, `--lat-name`, `--lon-name`
- Alignment: `--regrid`
- Styling: `--threshold`, `--vmin`, `--vmax`, `--vmax-percentile`, `--log-scale`, `--alpha`, `--cmap-flood`
- Hillshade: `--azdeg`, `--altdeg`, `--vert-exag`
- Output: `--out`, `--dpi`
- Debug: `--log-level DEBUG`

---

## Notes

- The script assumes **rectilinear grids** (1D lat, 1D lon). If you have curvilinear 2D lon/lat, you’ll need a different plotting approach.
- If you see misalignment, try `--regrid flood_to_dem` and check that both datasets use the same coordinate reference (usually lon/lat in EPSG:4326).

