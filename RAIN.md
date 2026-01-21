# Rainfall inputs for LPERFECT

LPERFECT consumes rainfall forcing files that follow the CF-1.10 schema in `cdl/rain_time_dependent.cdl`. Use the following utilities and guides to prepare rainfall inputs from different data sources.

## Conversion utilities

- **Weather radar GeoTIFFs → rain rate NetCDF:** [`docs/wr_to_rain.md`](docs/wr_to_rain.md)
- **WRF accumulations → rain rate NetCDF:** [`docs/meteouniparthenope-wrf_to_rainrate.md`](docs/meteouniparthenope-wrf_to_rainrate.md)
- **Weather station CSVs → rain rate NetCDF:** [`utils/ws_to_rainrate.md`](utils/ws_to_rainrate.md)

## Data preparation overview

For general guidance on rainfall forcing datasets (time-dependent, static, and multi-source bundles), see [`docs/data.md`](docs/data.md).
