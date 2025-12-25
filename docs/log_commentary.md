# Commented run log (2025-12-25)

This document annotates the provided `python main.py` run log to clarify what each block reports during the LPERFECT simulation.

## Startup and configuration

- `Domain 'domain_1' loaded (serial): (15600, 16800)`: the model read the domain grid in serial mode; the grid dimensions are 15,600 x 16,800 cells.
- `Starting simulation for domain 'domain_1' (1/1)`: only one domain is being simulated in this run.
- `Compute device: cpu`: the simulation uses the CPU (no GPU acceleration).
- `Travel time mode=auto ...`: the model auto-computed median travel times (hillslope ≈162 s, channel ≈54 s) for routing; used to set stable time stepping.
- `LPERFECT start: T_s=3600.0 dt_s=5.0 steps=720 ranks=1`: the run covers 1 hour (3,600 s) with a 5 s step, totaling 720 iterations on a single MPI rank.
- `Rain source 'rain': using rain rate at ...`: the forcing system selects the initial rainfall time slice (index 0) before stepping.

## Mass balance logging cadence

The `mass_balance` lines report conservation every 10 steps (50 s). Fields:
- `rain_m3`: cumulative rainfall volume injected.
- `runoff_m3`: water converted to runoff (leaving hillslopes).
- `outflow_m3`: volume that exited the domain boundary.
- `system_m3`: volume stored in the modeled system (surface/river cells).
- `err_m3`: numerical mass-balance error; values near 0 indicate stable computation.
- `hops`: particle hops processed; grows with flow routing complexity.

## Early-phase behavior (steps 10–120)

- Runoff stays zero through step 30, showing initial infiltration/storage delay.
- Runoff/outflow begin rising by step 40, with small storage accumulation (≈1.7 m³).
- By step 120 (10 minutes), runoff reaches ~62 m³ while outflow is still small; storage is ~1.75×10³ m³ and errors remain near machine precision.

## Rain forcing updates

Rainfall time slices change at:
- 14:37:23 (index=1)
- 14:42:14 (index=2)
- 14:47:13 (index=3)
- 14:52:19 (index=4)
- 14:57:25 (index=5)

Each change corresponds to a new rainfall rate applied in subsequent steps; notice slight drops in `rain_m3` increments after some updates.

## Mid-run hydrologic response (steps 130–420)

- Post first rain update, storage continues growing, indicating catchment filling while outflow increases gradually.
- Runoff growth accelerates with each rain period; e.g., by step 300 storage surpasses 11,000 m³ and outflow is ~2.8×10¹ m³.
- After the third rain update (index=3), runoff climbs more steeply, reflecting heavier rainfall; outflow and storage both trend upward without instability.

## Late run and peak values (steps 430–720)

- Storage exceeds 37,000 m³ by step 440 and continues increasing as rainfall persists.
- Following the final rain updates (indexes 4 and 5), cumulative rain increments decrease slightly, but the system keeps filling; storage reaches ~1.62×10⁵ m³ by step 720.
- Mass-balance error remains extremely small (|err_m3| ≤ O(10⁻¹⁰)), indicating stable numerical integration.

## Completion and quality report

- The simulation ends at step 720 (1 hour) and prints a quality summary.
- Volume balance confirms total rain ≈1.446×10⁸ m³; runoff and outflow are orders of magnitude smaller, as expected for this setup.
- Final mass error is −5.821×10⁻¹¹ m³ (−0.000% of runoff), with low mean error—evidence of near-perfect conservation.
- Particle statistics show ~27.9 million particles spawned, ~5.85 million exiting, ~22 million still active, and ~130 million total hops—consistent with dense routing activity.
- Output GeoJSON `data/20251223Z1400_outflow.geojson` contains 24,605 features representing the modeled outflow.
- `LPERFECT finished across 1 domain(s).`: confirms clean shutdown across all ranks.
