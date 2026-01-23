# `utils/output_plot.py` — Plot LPERFECT flood outputs over DEM hillshade

This script creates a map figure by combining:
- **Hillshade** from a DEM (domain NetCDF), and
- **Flood output overlays** from the LPERFECT NetCDF (depth, risk, or inundation masks).

It includes the “do it all” improvements:
- FillValue/NaN handling
- Threshold masking
- Percentile-based `vmax` for better contrast
- Optional log scaling
- Dry pixels are transparent (any `flood_depth <= 0` is masked automatically)
- Variable-aware defaults for flood depth vs. risk/mask layers
- Optional regridding (grid alignment)
- Batch frame generation for all time steps
- Optional GeoJSON/Shapefile overlay (if installed), with centroid labels
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

## Choose a variable to plot

Supported variables (use `--plot-var`):
- `flood_depth`
- `risk_index`
- `inundation_mask`
- `flood_depth_max`
- `inundation_mask_max`

Example (risk index with the built-in green→yellow→orange→red scale):
```bash
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --plot-var risk_index \
  --out outputs/risk_index.png
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
  --plot-var flood_depth \
  --log-scale \
  --out outputs/flood_log.png
```

Clip to a bounding box (e.g., Regione Campania) to speed up rendering and focus the map:
```bash
python utils/output_plot.py \
  --flood data/flood_depth.nc \
  --domain data/domain.nc \
  --bbox 13.7 39.9 15.9 41.6 \
  --out outputs/flood_campania.png
```
Order is `min_lon min_lat max_lon max_lat` (EPSG:4326).

Common Campania-area boxes you can pass to `--bbox`:
- **Regione Campania:** `13.7 39.9 15.9 41.6` (broad regional view)
- **Naples Province (Metropolitan City):** `13.8 40.6 14.9 41.2` (includes islands)
- **City of Naples:** `14.14 40.78 14.30 40.88` (urban core)
- **Ischia Island:** `13.86 40.69 13.99 40.76` (island detail)

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
  --overlay-label-field name \
  --out outputs/flood_with_bounds.png
```

- Use `--overlay-notext data/boundaries.geojson` for the same overlay without labels.
- If `--overlay-label-field` is omitted, the script tries common name columns (`name`, `nome`, etc.) or the first string-like column.
- Labels are placed at feature centroids with a light outline for readability.
- Overlays are clipped to the plotting domain bounds (including any `--bbox` selection).

---

## Key CLI options

- Inputs: `--flood`, `--domain`
- Time: `--time-index`, `--all-times --out-dir`
- Variables: `--plot-var` (alias `--flood-var`), `--dem-var`, `--lat-name`, `--lon-name`
- Alignment: `--regrid`
- Styling: `--threshold`, `--vmin`, `--vmax`, `--vmax-percentile`, `--log-scale`, `--alpha`, `--cmap-flood`
- Subset: `--bbox min_lon min_lat max_lon max_lat`
- Overlays: `--overlay`, `--overlay-notext`, `--overlay-label-field`, `--overlay-label-size`, `--vector-alpha`, `--vector-linewidth`
- Hillshade: `--azdeg`, `--altdeg`, `--vert-exag`
- Output: `--out`, `--dpi`
- Debug: `--log-level DEBUG`

---

## Notes

- The script assumes **rectilinear grids** (1D lat, 1D lon). If you have curvilinear 2D lon/lat, you’ll need a different plotting approach.
- If you see misalignment, try `--regrid flood_to_dem` and check that both datasets use the same coordinate reference (usually lon/lat in EPSG:4326).
