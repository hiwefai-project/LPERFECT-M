# Model Description
**Lagrangian Parallel Environmental Runoff and Flood Evaluation for Computational Terrain**

---

## 1. Overview

LPERFECT is a **Lagrangian, particle-based hydrological model** designed to estimate surface runoff, flood propagation, and hydrogeological risk over gridded terrains.  
The model is optimized for **parallel execution using MPI**, supports **restartable simulations**, and is conceived for **operational and research workflows** within the **Hi-WeFAI** (*High-Performance Weather and Flood AI Systems*) project.

LPERFECT is intentionally positioned between very simple runoff-index models and full two-dimensional hydrodynamic solvers. Its goal is to provide a **computationally efficient, physically interpretable flood-impact model** that can be executed rapidly, repeatedly, and at scale.

Typical application domains include:
- flood nowcasting and early warning,
- ensemble and scenario-based simulations,
- AI-driven rainfall impact assessment,
- HPC and cloud-based operational chains.

---

## 2. Conceptual Architecture

LPERFECT follows a modular processing pipeline:

1. **Rainfall ingestion**, possibly from multiple heterogeneous sources, blended in a time-aware manner.
2. **Runoff generation** based on the cumulative SCS Curve Number (CN) method.
3. **Particle spawning**, where incremental runoff is converted into discrete Lagrangian water parcels.
4. **Lagrangian routing** over a D8 flow-direction network with configurable travel times.
5. **Parallel particle migration** across MPI subdomains.
6. **Flood-depth reconstruction** from particle volumes.
7. **Hydrogeological risk assessment**.
8. **CF-compliant NetCDF outputs** and optional restart checkpoints.

All spatial information is handled exclusively via **NetCDF datasets**, ensuring consistency and interoperability.

---

## 3. Lagrangian Representation of Surface Water

Surface water is represented using **discrete particles**, each carrying a fixed reference volume.  
Particles are independent entities characterized by:
- grid position (row, column),
- water volume,
- an internal travel-time counter.

The Lagrangian formulation has several advantages:
- mass conservation is guaranteed by construction,
- no global Courant–Friedrichs–Lewy (CFL) constraint is required,
- particle motion is naturally parallelizable,
- routing logic remains simple and robust.

Particles are generated locally from runoff excess and move independently along predefined flow paths.

---

## 4. Runoff Generation Model

LPERFECT adopts the **SCS Curve Number (CN)** method in cumulative form.

The potential maximum retention is defined as:

\[
S = \frac{25400}{CN} - 254
\]

The initial abstraction is:

\[
I_a = \alpha S
\]

where \( \alpha \) is a configurable initial abstraction ratio.

Cumulative runoff is computed as:

\[
Q =
\begin{cases}
0 & P \le I_a \\
\frac{(P - I_a)^2}{P - I_a + S} & P > I_a
\end{cases}
\]

where:
- \(P\) is cumulative precipitation,
- \(Q\) is cumulative runoff.

Incremental runoff at each timestep is obtained by differencing cumulative runoff values.

---

## 5. D8-Based Lagrangian Routing

Routing is performed over a **D8 flow-direction network**, where each grid cell drains to a single downstream neighbor.

Two common encodings are supported:
- ESRI encoding (1, 2, 4, 8, 16, 32, 64, 128),
- clockwise encoding (0–7).

Each particle advances downstream only when its internal travel-time counter reaches zero. After each hop:
- a hillslope travel time is assigned for non-channel cells,
- a shorter channel travel time is assigned for channel cells,
- travel times can be fixed (user-supplied scalars) or automatically derived from cell area and representative hillslope/channel velocities, bounded by configurable minimum/maximum values.

This mechanism allows sub-timestep control of motion and captures faster transport along drainage networks without explicit hydraulic equations.

---

## 6. Parallelization Strategy

LPERFECT uses **row-slab domain decomposition** for MPI parallelization.

Key principles:
- each MPI rank owns a contiguous block of grid rows,
- particles belong to the rank owning their current row,
- after advection, particles crossing subdomain boundaries are migrated to the correct rank using collective MPI communication.

This strategy:
- minimizes communication overhead,
- avoids global synchronization during routing,
- scales efficiently for large domains and large particle counts.

All NetCDF input and output operations are performed by **rank 0**, while computational work is distributed across ranks.

---

## 7. Restartable State Model

LPERFECT supports **fully restartable simulations** using NetCDF checkpoint files.

The restart state includes:
- cumulative precipitation and runoff fields,
- particle positions, volumes, and travel-time counters,
- elapsed simulation time,
- mass-balance diagnostics.

On restart:
- rank 0 reads the checkpoint,
- state fields are scattered to MPI ranks,
- particles are redistributed according to spatial ownership.

This design supports fault tolerance, long-running workflows, and ensemble execution.

---

## 8. Flood-Depth Reconstruction

Flood depth is reconstructed by aggregating particle volumes on the grid.

For each grid cell:

\[
h(i,j) = \frac{\sum V_p(i,j)}{A(i,j)}
\]

where:
- \(V_p(i,j)\) is the volume of particles in the cell,
- \(A(i,j)\) is the cell area.

The resulting flood-depth field is expressed in meters.

---

## 9. Hydrogeological Risk Index

LPERFECT computes a **dimensionless hydrogeological risk index** combining:
1. direct hydrological forcing (cumulative runoff),
2. morphological control (flow accumulation).

The risk index is defined as:

\[
R = \beta \, \hat{Q} + (1 - \beta) \, (\hat{Q} \cdot \hat{A})
\]

where:
- \( \hat{Q} \) is normalized runoff,
- \( \hat{A} \) is normalized flow accumulation,
- \( \beta \) is a configurable balance parameter.

Normalization is performed using robust percentile thresholds to reduce sensitivity to outliers.

---

## 10. Input and Output Data Model

### Inputs
- Domain NetCDF containing DEM, D8, CN, and optional channel mask.
- One or more rainfall NetCDF datasets following `cdl/rain_time_dependent.cdl`.

### Outputs
- flood depth field, expressed in meters,
- hydrogeological risk index, dimensionless,
- optional GeoJSON with outflow hit points when enabled, listing the sea/lake impact cells and particle counts per save interval,
- a logged simulation quality report summarizing mass-balance and hydrological checks.

All outputs follow **CF-1.10 conventions**, ensuring compatibility with standard geoscientific tools.

---

## 11. Scope and Limitations

LPERFECT is:
- fast, scalable, and mass-conservative,
- suitable for ensemble and nowcasting applications,
- designed for HPC and cloud environments.

LPERFECT is **not**:
- a full shallow-water or Navier–Stokes solver,
- a substitute for detailed hydraulic modeling in urban-scale studies.

---

## 12. Outlook

Planned developments include:
- distributed NetCDF I/O,
- GPU acceleration,
- alternative infiltration models,
- coupling with hydraulic solvers,
- uncertainty propagation and probabilistic risk assessment.

---

**LPERFECT** provides a computationally efficient bridge between rainfall intelligence and flood-risk awareness, supporting the objectives of the Hi-WeFAI project.
