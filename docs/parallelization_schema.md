# Parallelization Schema (MPI + GPU)

This document summarizes how LPERFECT-M parallelizes work across distributed-memory (MPI)
processes, uses shared-memory threads within each rank, and how GPU acceleration is applied when available.

---

## 1. Distributed-Memory Parallelism (MPI)

LPERFECT-M uses **row-slab domain decomposition** for MPI execution.
Each rank owns a contiguous block of grid rows and performs all particle updates that
occur within its slab.

### 1.1 Domain Ownership

- The global grid is split by **rows** into `N` slabs.
- Rank `r` owns rows `[row_start_r, row_end_r]` and all variables restricted to that slab.
- Particles are **owned by the rank that owns their current row**.

### 1.2 Particle Migration

After each routing step:

1. Particles are advanced along the D8 flow direction.
2. Any particle that crosses a slab boundary is placed into a **send buffer**.
3. Ranks exchange boundary-crossing particles with their neighbors.
4. Received particles are merged into the local particle list.

This approach keeps communication localized and avoids global all-to-all exchanges.

### 1.3 I/O Strategy

- **Rank 0** performs all NetCDF input and output.
- Input fields are read on rank 0 and **scattered** to other ranks.
- Output fields are gathered to rank 0 and written to disk.
- Restart files follow the same pattern (rank 0 writes, rank 0 reads and redistributes).

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
