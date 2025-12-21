# LPERFECT
**Lagrangian Parallel Environmental Runoff and Flood Evaluation for Computational Terrain**
## Data preparation

A detailed guide to preparing **domain** and **rainfall forcing** NetCDF inputs (with GDAL/CDO/NCO examples and Italy-oriented datasets) is provided here:

- [`data.md`](../data.md) (in this repository root)


LPERFECT is a lightweight **Lagrangian runoff + flood routing model** designed for fast flood screening on gridded domains.
It is engineered to be **restartable**, **MPI-parallel** (row-slab + particle migration), and **production-ready** from a software-architecture point of view (modular package layout, clear separation of concerns).

---

## 1) Context in Hi-WeFAI

LPERFECT fits naturally in the **Hi-WeFAI** workflow, where blended precipitation forcing (radar nowcasting, stations, NWP)
feeds a fast impact layer producing flood depth and risk maps.

Project website: https://www.hiwefai-project.org

---

## 2) Project structure

```
LPERFECT/
  main.py                  # Thin entry point (CLI + orchestration)
  config.json              # Example configuration
  requirements.txt         # Python dependencies
  lperfect/                # Python package
    __init__.py
    cli.py                 # CLI parsing
    config.py              # default config + merging
    logging_utils.py       # logging setup
    time_utils.py          # time helpers (UTC, ISO parsing)
    domain.py              # domain NetCDF reader + broadcast
    d8.py                  # D8 downstream lookup
    runoff.py              # SCS-CN runoff
    rain.py                # rain inputs + blending (rank0 read)
    particles.py           # particle data structure + pack/unpack
    mpi_utils.py           # slab decomposition + migration + scatter/gather
    hydraulics.py          # particle spawning/advection + volume grid
    risk.py                # flow accumulation + risk index
    io_netcdf.py           # output + restart NetCDF I/O (rank0 only)
    simulation.py          # main simulation driver
```

All Python files include **pedantic, educational comments**.

---

## 3) Installation

```bash
pip install -r requirements.txt
```

MPI runs require an MPI runtime (OpenMPI / MPICH) and `mpi4py`.

---

## 4) Run

### Serial
```bash
python main.py --config config.json
```

### MPI
```bash
mpirun -np 8 python main.py --config config.json
```

Useful overrides:
```bash
python main.py --config config.json --out-nc flood.nc
python main.py --config config.json --restart-in restart_state.nc
python main.py --config config.json --restart-out restart_state.nc
```

---

## 5) Input NetCDFs

### 5.1 Domain NetCDF (`domain.domain_nc`)
Required 2D variables (dims `y,x`):
- `dem` (float): DEM elevation (m). Finite values define the **active domain**.
- `d8` (int): D8 directions (`model.encoding` = `esri` or `cw0_7`).
- `cn` (float): Curve Number (0–100). Values outside become 0.

Optional 2D variable:
- `channel_mask` (0/1 or bool): Cells treated as channel (faster travel time).

Required coordinates (names configurable via `domain.varmap`):
- `x(x)` and `y(y)` with CF-style `units` if possible.

Optional CF grid mapping:
- If `dem` has attribute `grid_mapping`, and the referenced variable exists, LPERFECT preserves it in outputs.

### 5.2 Rain NetCDFs
Each rain source can be:
- 2D `(y,x)` static
- 3D `(time,y,x)`

Modes:
- `intensity_mmph` (mm/h)
- `depth_mm_per_step` (mm per timestep)

Time selection:
- `nearest`: pick the time slice closest to simulation timestamp (`model.start_time` + elapsed)
- `step`: pick `time[k]` by index

---

## 6) Model description

### 6.1 Runoff generation (SCS Curve Number)
LPERFECT uses cumulative SCS-CN:
- accumulate precipitation `P_cum` (mm)
- convert to cumulative runoff `Q_cum` (mm)
- incremental runoff is `dQ = Q_cum(t) - Q_cum(t-1)`

### 6.2 Lagrangian routing
Incremental runoff is converted into **particles** (small fixed water volumes).  
Particles hop along D8 with travel-time gating.

### 6.3 Parallelization
MPI mode uses:
- **slab decomposition** (contiguous row blocks per rank)
- **particle migration** via `Alltoallv` after advection
- **rank0-only I/O** (domain read, rain read, restart+output writes)

---

## 7) Restart usage

Start fresh (and write restart every `restart.every` steps):
```bash
python main.py --config config.json
```

Resume:
```bash
python main.py --config config.json --restart-in restart_state.nc
```

MPI restarts:
- rank0 loads restart
- cumulative fields are scattered as slabs
- particles are redistributed by row ownership

---

## 8) Outputs

Results NetCDF contains:
- `flood_depth(y,x)` in meters
- `risk_index(y,x)` unitless (if enabled)

Example plot:
```python
import xarray as xr
import matplotlib.pyplot as plt

ds = xr.open_dataset("flood_depth.nc")
ds["flood_depth"].plot()
plt.title("Flood depth (m)")
plt.show()

ds["risk_index"].plot()
plt.title("Risk index")
plt.show()
```

---
