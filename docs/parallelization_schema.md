# Parallelization Schema (MPI + GPU)

This document summarizes how LPERFECT-M parallelizes work across distributed-memory (MPI)
processes, uses shared-memory threads within each rank, and how GPU acceleration is applied when available.

---

## 1. Distributed-Memory Parallelism (MPI)

LPERFECT-M offers two MPI schemas:

1. **Slab (default)**: load-balanced row-slab decomposition (described below).
2. **Particles**: domain fields are replicated; particles are evenly distributed across ranks (and rebalanced automatically).

If MPI is available but explicitly disabled, only rank 0 runs while other launched ranks exit immediately.
When `model.runoff_only_risk=true` (or `--runoff-only-risk`), particle transport is skipped and particle-specific schema behavior is effectively inactive.

### 1.1 Domain Ownership (Slab schema)

- The global grid is split by **rows** into `N` slabs, weighted by the number of active cells per row.
- Rank `r` owns rows `[row_start_r, row_end_r]` and all variables restricted to that slab.
- Particles are **owned by the rank that owns their current row**.
- Guardrails: `compute.mpi.min_rows_per_rank` avoids tiny slabs; set `compute.mpi.decomposition` to
  `balanced` (default) or `even` to force uniform slabs.

### 1.2 Particle Migration (Slab schema)

After each routing step:

1. Particles are advanced along the D8 flow direction.
2. Any particle that crosses a slab boundary is placed into a **send buffer**.
3. Ranks exchange boundary-crossing particles with their neighbors. When no particles leave a slab, the exchange is skipped entirely to avoid empty collectives.
4. Received particles are merged into the local particle list.

This approach keeps communication localized and avoids global all-to-all exchanges; a neighbor-only path is used whenever all migrants target adjacent slabs.

### 1.3 I/O Strategy

- **Rank 0 only** (default) or **all ranks** (configurable):
  - When `compute.parallelization.io="rank0"`, rank 0 handles domain/rain I/O, gathers outputs, and writes a single NetCDF.
  - When `compute.parallelization.io="all"`, each rank writes its own artifacts with `_rankXXXX` suffixes while still keeping hydrological totals synchronized.
- In slab mode, input fields (restart, domain, and per-step rainfall) are read on rank 0 and **scattered** to other ranks,
  avoiding whole-grid broadcasts every timestep. Output fields are gathered to rank 0 and written to disk (or written per-rank if requested).
- Restart files follow the same pattern: centralized by default, optional per-rank outputs when configured.

### 1.4 Independent switches

- Toggle MPI on/off with `compute.mpi.enabled = true|false|null` or the CLI `--mpi-mode enabled|disabled|auto`.
- Select decomposition with `compute.mpi.decomposition = balanced|even` or `--mpi-decomposition` (slab schema only).
- Control slab granularity via `compute.mpi.min_rows_per_rank` or `--mpi-min-rows`. Ranks are automatically pruned when the active rows cannot meet this constraint, keeping only the leading ranks needed for the domain.
- Periodic/automatic rebalancing: `compute.mpi.balance.every_steps` or `every_sim_s` forces a resplit that accounts for current particle density; `compute.mpi.balance.auto` rebalances when the max/min particle ratio exceeds `imbalance_threshold`.
- Choose schema via `compute.parallelization.schema = slab|particles` (CLI: `--parallel-schema`). Set I/O mode with `compute.parallelization.io = rank0|all` (CLI: `--parallel-io`).

### 1.5 Particle-Only Schema

- Particles are evenly distributed across ranks (and rebalanced automatically when a rain source advances to a new time slice, plus any configured periodic rebalance).
- Domain fields (`dem`, `d8`, `cn`, masks) remain replicated on all ranks for hydrological consistency.
- Rainfall is read on rank 0, broadcast to all ranks, and new particles are spawned on rank 0 before **per-step even redistribution** that scatters only the fresh particles (existing particles stay in place). This preserves balance without a full gather/rebroadcast cycle.
- Advection, shared-memory threading, and GPU usage remain per rank; no slab-based migration is needed because ownership is particle-only.

---

## 2. GPU Acceleration (Per Rank)

LPERFECT-M optionally uses **CuPy** to accelerate array-based computations on GPUs.
GPU acceleration is **local to each MPI rank** and does not change the MPI
communication pattern.

### 2.1 What Runs on the GPU

When `compute.device = "gpu"` and CuPy is available:

- The SCS-CN runoff calculations operate on GPU arrays.
- Intermediate cumulative precipitation and runoff fields are stored on the GPU
  during runoff calculation.

All MPI communication and particle routing remain CPU-based.

### 2.2 Fallback Behavior

If the GPU device is requested but CuPy is unavailable:

- LPERFECT-M **logs a warning** and falls back to CPU execution.
- The simulation continues with NumPy arrays.

---

## 3. Shared-Memory Parallelism (Per Rank)

LPERFECT-M can accelerate CPU-bound sections inside each rank using **threaded shared-memory**
parallelism. This targets particle advection and slab accumulation, which complements GPU
acceleration (runoff) without requiring a GPU.

- Enable via `compute.shared_memory.enabled = true`.
- Configure worker threads with `compute.shared_memory.workers` (defaults to CPU cores).
- Guardrails: `min_particles_per_worker` avoids overhead on small local particle counts, and
  `chunk_size` tunes per-task batch size for cache efficiency.
- Results are **deterministic**: chunks are concatenated in submission order to preserve
  reproducibility while reducing contention.

Shared-memory threading is independent of GPU usage: you can run CPU-only (threads),
GPU-only, or combine GPU runoff with threaded particle routing on the host.

---

## 4. Hybrid MPI + GPU Execution

For hybrid runs:

- Each MPI rank can target a GPU on its node.
- GPU acceleration is **independent per rank**.
- Users are responsible for setting appropriate GPU affinity (e.g., with
  `CUDA_VISIBLE_DEVICES` or MPI runtime options).

Example:

```bash
mpirun -np 4 python main.py --config config.json --device gpu
```

---

## 5. Summary

- **MPI** distributes the grid by row slabs and migrates particles across boundaries.
- **Shared-memory threads** speed up particle advection and slab accumulation within each rank.
- **GPU** acceleration (via CuPy) speeds up runoff calculations within each rank.
- **I/O** is centralized on rank 0 for simplicity and CF-compliant NetCDF outputs.
- Each layer (MPI, threads, GPU) is independently switchable; you can run GPU-only, threads-only, MPI-only, or any combination.

---

## 6. Parallelization evaluation metrics

The table below captures the metrics we track when validating the combined MPI + threading + GPU implementation. They come from a 1-hour synthetic rainfall experiment (5 s time step) using row-slab decomposition; values illustrate typical behavior on a dual-socket CPU node.

| Metric | What it measures | Example value | How to interpret |
| --- | --- | --- | --- |
| Wall-clock per step (CPU-only) | Mean elapsed time per routing step on 1 rank | **96 ms/step** | Baseline for strong-scaling comparisons. |
| Strong-scaling speedup (4 ranks) | Wall-clock(1 rank) / Wall-clock(4 ranks) | **3.4×** | Near-ideal scaling for slab-balanced domains (≈85 % efficiency). |
| Particle throughput | Particle hops processed per second | **≈1.9×10⁶ hops/s** | Tracks how well threading accelerates local advection/aggregation. |
| Migration ratio | Fraction of particles exchanged across ranks per step | **5–8 %** | Indicates communication pressure; higher ratios warrant larger slabs or topology-aware placement. |
| GPU runoff speedup | Runtime improvement of runoff kernels vs CPU | **2.7×** | Applies only when `compute.device="gpu"`; particle routing stays on CPU. |
| Memory per rank | Peak resident set size during the run | **≈3.2 GB** | Helps size jobs for many-rank deployments; scales roughly with local slab height. |

Measurement notes:
- Wall-clock timings come from the simulation log; mean values exclude startup/I/O. When using GPUs, only the runoff section benefits, so whole-step speedup will be lower than kernel-only speedup.
- Particle throughput is derived from the `hops` counter divided by wall-clock runtime; it captures the combined effect of threading and vectorization within each rank.
- Migration ratio is computed as `(migrated / total active particles) * 100` per step and helps identify decomposition changes before scaling to many nodes.
- Slab-mode steps with zero migrants now skip Alltoall entirely, and particle-mode runs scatter only the newly spawned particles each step; both behaviors reduce synchronization noise when rain is sparse.
- If `compute.shared_memory.enabled=true`, ensure `workers` matches the cores allocated by the launcher (e.g., `SLURM_CPUS_PER_TASK`) to reproduce the speedups above.

### Performance-tuned compute block (CPU, 4 MPI ranks)

For large national domains (e.g., 15 600 × 16 800) on a CPU-only node with four MPI ranks, use this `compute` section to minimize communication while keeping ranks busy:

```json
"compute": {
  "device": "cpu",
  "mpi": {
    "enabled": true,
    "decomposition": "balanced",
    "min_rows_per_rank": 512,
    "balance": {
      "every_steps": 120,
      "every_sim_s": 0,
      "auto": true,
      "imbalance_threshold": 1.5
    }
  },
  "shared_memory": {
    "enabled": true,
    "workers": null,
    "min_particles_per_worker": 2000,
    "chunk_size": 16384
  }
}
```

- `balanced` slabs plus `min_rows_per_rank` keep ranks from owning tiny slices, which reduces migration pressure and amortizes I/O per slab.
- `balance.auto` with a 1.5 threshold re-splits when particle load skews; `every_steps=120` forces a check every 10 minutes (with `dt_s=5 s`) for long events.
- Neighbor-only migration kicks in automatically when particles cross only adjacent slabs, shrinking MPI traffic without changing results.
- Shared-memory threading is enabled per rank and defaults `workers` to hardware cores; the lower `min_particles_per_worker` and smaller `chunk_size` improve throughput when particle counts ramp up gradually.

### How LPERFECT now computes and exports metrics

LPERFECT can emit a GPT-friendly JSON report with per-step and summary metrics that work across:

- **Shared memory:** per-step throughput and wall times reflect thread usage and chunking inside each rank.
- **Distributed memory:** migration ratios and per-step wall times use the MPI max across ranks to capture the slowest slab.
- **GPU + CPU hybrids:** runoff timings isolate the GPU-heavy CN computation, while advection/migration timings stay CPU-side.

Enable the report via either `--parallel-metrics` (CLI) or `metrics.parallelization.enabled=true` (JSON). Optional fields:

- `metrics.parallelization.output`: path to write the JSON (also logged on rank 0).
- `metrics.parallelization.max_samples`: cap on per-step samples; large runs are evenly down-sampled.

For AI summaries that focus on hydrology and compute settings (rather than step-by-step timing), enable `metrics.assistant` or the `--ai-metrics` CLI flag to emit a compact JSON status report (documented in `docs/ai.md`).

The JSON contains:

- `scenario`: ranks, shared-memory settings, device, domain shape, timestep.
- `summary`: wall-clock totals and quantiles, particle throughput, migration ratios, segment timings (runoff/advect/migration), spawned/outflow particle counts.
- `per_step_samples`: sampled steps with wall time, segment breakdown, hops, throughput, and migration ratios.

This structure is stable and designed for downstream agents to ingest directly for automated performance tuning across CPU, GPU, and MPI configurations.
