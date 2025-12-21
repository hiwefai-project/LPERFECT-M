# LPERFECT – Model Description
**Lagrangian Parallel Environmental Runoff and Flood Evaluation for Computational Terrain**

---

## 1. Overview

LPERFECT is a **Lagrangian, particle-based hydrological model** designed to estimate surface runoff, flood propagation, and hydrogeological risk over gridded terrains.
It is optimized for **parallel execution (MPI)** and **operational workflows**, and is fully integrated in the **Hi-WeFAI** (*High-Performance Weather and Flood AI Systems*) project.

Unlike full hydrodynamic solvers, LPERFECT focuses on **computational efficiency, robustness, and scalability**, making it suitable for:
- flood nowcasting,
- ensemble simulations,
- AI-driven rainfall impact assessment,
- HPC and cloud deployments.

---

## 2. Conceptual Architecture

![LPERFECT conceptual pipeline](figures_svg/lperfect_pipeline.svg)

The LPERFECT processing pipeline is organized as follows:

1. **Rainfall ingestion** from multiple, heterogeneous sources (radar, stations, NWP), blended in a time-aware fashion.
2. **Runoff generation** using the cumulative SCS Curve Number (CN) method.
3. **Particle spawning**, converting incremental runoff into discrete water parcels.
4. **Lagrangian routing** over a D8 flow network with configurable travel times.
5. **MPI-based particle migration** across domain partitions.
6. **Flood depth reconstruction** and **hydrogeological risk assessment**.
7. **CF-compliant NetCDF outputs** for downstream analysis.

---

## 3. Lagrangian Representation of Surface Water

Water is represented by **discrete particles**, each carrying a fixed reference volume.

Key properties:
- particles are created locally from runoff excess,
- particles move independently along D8-defined paths,
- mass conservation is ensured by construction,
- no global CFL constraint is required.

![Particle routing](figures_svg/geographical_example.svg)

The figure above illustrates the full chain **DEM → D8 → particles → flood depth** on a simplified geographical example.

---

## 4. Runoff Generation Model

LPERFECT adopts the **SCS Curve Number (CN)** method in cumulative form.

### Governing equations

Potential maximum retention:

\[
S = \frac{25400}{CN} - 254
\]

Initial abstraction:

\[
I_a = \alpha S
\]

Cumulative runoff:

\[
Q =
\begin{cases}
0 & P \le I_a \\
\frac{(P - I_a)^2}{P - I_a + S} & P > I_a
\end{cases}
\]

Where:
- \(P\) is cumulative precipitation (mm),
- \(Q\) is cumulative runoff (mm),
- \(\alpha\) is the initial abstraction ratio.

Incremental runoff is obtained as:

\[
\Delta Q = Q(t) - Q(t-1)
\]

---

## 5. D8-Based Lagrangian Routing

Each grid cell routes particles to **one downstream neighbor**, defined by a D8 flow-direction raster.

Supported encodings:
- **ESRI** (1, 2, 4, 8, 16, 32, 64, 128)
- **Clockwise** (0–7)

### Travel time control

Particles advance only when their internal clock \(\tau \le 0\).
After each hop:
- hillslope cell: \(\tau = t_{hill}\),
- channel cell: \(\tau = t_{channel}\).

This mechanism enables sub-timestep stability and channel acceleration without explicit hydrodynamics.

---

## 6. Parallelization Strategy (MPI)

![MPI slab decomposition](figures_svg/mpi_slabs.svg)

LPERFECT uses **row-slab domain decomposition**:

- each MPI rank owns a contiguous block of rows,
- particles belong to the rank owning their current row,
- particles crossing slab boundaries are migrated using `MPI_Alltoallv`.

This strategy ensures:
- minimal communication overhead,
- excellent weak scalability,
- natural compatibility with the Lagrangian formulation.

---

## 7. Restartable State Model

The full simulation state is checkpointed in **NetCDF format**, enabling restart and workflow chaining.

![Restart workflow](figures_svg/lperfect_architecture.svg)

Stored state includes:
- cumulative precipitation and runoff fields,
- particle positions, volumes, and timers,
- elapsed simulation time,
- mass balance diagnostics.

---

## 8. Flood Depth Reconstruction

Flood depth is reconstructed by aggregating particle volumes:

\[
h(i,j) = \frac{\sum V_p(i,j)}{A(i,j)}
\]

Where:
- \(V_p\) is the particle volume,
- \(A(i,j)\) is the grid-cell area.

The result is a spatially distributed flood-depth field.

---

## 9. Hydrogeological Risk Index

LPERFECT computes a **dimensionless hydrogeological risk index** by combining:

1. **Direct hydrological forcing** (cumulative runoff),
2. **Morphological control** (flow accumulation).

### Definition

\[
R = \beta \, \hat{Q} + (1 - \beta) \, \hat{A}
\]

Where:
- \(\hat{Q}\) is normalized runoff,
- \(\hat{A}\) is normalized flow accumulation,
- \(\beta\) is a configurable balance parameter.

Robust normalization is performed using percentile thresholds.

---

## 10. Input and Output Data Model

### Inputs
- Domain NetCDF: DEM, D8, CN, optional channel mask
- Rainfall NetCDFs: time-dependent, multi-source

### Outputs
- `flood_depth(y,x)` [m]
- `risk_index(y,x)` [-]

All outputs are **CF-1.10 compliant**, ensuring interoperability with standard geoscientific tools.

---

## 11. Example Execution

```bash
mpirun -np 8 python main.py --config config.json
```

```python
import xarray as xr
ds = xr.open_dataset("flood_depth.nc")
ds.flood_depth.plot()
```

---

## 12. Scope and Limitations

LPERFECT is:
- fast, scalable, and mass conservative,
- suitable for ensemble and nowcasting applications,
- designed for HPC and cloud environments.

LPERFECT is **not**:
- a full shallow-water or Navier–Stokes solver,
- a replacement for detailed hydraulic models.

---

## 13. Outlook

Planned extensions include:
- distributed NetCDF I/O,
- GPU acceleration,
- alternative infiltration models,
- coupling with hydraulic solvers,
- uncertainty propagation and probabilistic risk metrics.

---

**LPERFECT** provides a computationally efficient bridge between **rainfall intelligence** and **flood-risk awareness**, supporting the goals of the Hi-WeFAI project.
