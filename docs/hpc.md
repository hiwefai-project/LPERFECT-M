# Running LPERFECT-M on Slurm HPC systems

This guide explains how to launch **LPERFECT-M** on a high-performance computing (HPC) cluster managed by **Slurm**. It covers CPU-only, shared-memory, MPI-distributed, GPU-accelerated, and hybrid configurations, with ready-to-use `sbatch` templates and tuning tips.

---

## Prerequisites

- **Python environment** with `mpi4py` (for MPI jobs) and **CuPy** (for GPU jobs).
- **MPI runtime** (OpenMPI/MPICH) available via modules or environment.
- **CUDA** or ROCm drivers and runtime on GPU nodes.
- Access to a **shared filesystem** that is visible from all nodes for the domain/rain NetCDF inputs and outputs.
- A **JSON configuration** (e.g., `config.json`) tailored to your domain and rain sources. See `docs/configure.md` for full options.

Recommended module skeleton (adjust names to your cluster):

```bash
module load python/3.11
module load openmpi/4.1        # or mpich/4.x
module load cuda/12.2          # only for GPU nodes
```

---

## Slurm job-script template

Use `sbatch` to submit a script similar to:

```bash
#!/bin/bash
#SBATCH --job-name=lperfect
#SBATCH --time=02:00:00
#SBATCH --partition=compute          # adjust to your cluster
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1          # increase for MPI
#SBATCH --cpus-per-task=1            # increase for shared-memory threading
#SBATCH --gres=gpu:0                 # set to >0 for GPU runs
#SBATCH --output=logs/%x-%j.out

set -euo pipefail

module purge
module load python/3.11
module load openmpi/4.1
# module load cuda/12.2              # enable on GPU partitions

CONFIG=/path/to/config.json          # must be on shared storage
OUTDIR=/path/to/outputs
mkdir -p "${OUTDIR}"

# Example invocation (choose one of the modes below)
srun --mpi=pmix_v3 python main.py --config "${CONFIG}" --out-nc "${OUTDIR}/flood_depth.nc"
```

Adjust the `#SBATCH` directives and the final `srun` line according to the modes below.

---

## Execution modes

LPERFECT-M supports three kinds of parallelism:

1. **Shared memory (threads)** within each rank (`compute.shared_memory.enabled=true`).
2. **Distributed memory (MPI)** across ranks (`mpirun`/`srun` with `--ntasks` > 1).
3. **GPU acceleration** per rank (`compute.device="gpu"` or `--device gpu`), independent of MPI/threads.

You can enable them independently or in combination. Use the sections below as recipes.

### 1) Baseline CPU (single process)

- Slurm: `--ntasks=1`, `--cpus-per-task=1`, `--gres=gpu:0`.
- Command:

```bash
srun python main.py --config "${CONFIG}" --out-nc "${OUTDIR}/flood_depth.nc"
```

### 2) CPU + shared memory (threads only)

- Slurm: `--ntasks=1`, set `--cpus-per-task` to the desired thread count (e.g., 8).
- Config: enable threading in JSON, for example:

```json
{
  "compute": {
    "shared_memory": {
      "enabled": true,
      "workers": 8,
      "min_particles_per_worker": 20000
    }
  }
}
```

- Command:

```bash
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun python main.py --config "${CONFIG}"
```

### 3) MPI (distributed memory, CPU-only)

- Slurm (single node example): `--nodes=1`, `--ntasks-per-node=4`, `--cpus-per-task=1`.
- Slurm (multi-node example): `--nodes=2`, `--ntasks-per-node=4` (8 ranks total).
- Command:

```bash
srun --mpi=pmix_v3 python main.py --config "${CONFIG}"
```

### 4) MPI + shared memory (hybrid CPU)

- Slurm: choose ranks and threads, e.g., `--nodes=1`, `--ntasks-per-node=4`, `--cpus-per-task=8`.
- Config: set `compute.shared_memory.enabled=true` and `workers` to match `--cpus-per-task`.
- Tips:
  - Keep `workers` ≤ `SLURM_CPUS_PER_TASK` to avoid oversubscription.
  - Particle routing threads are confined to each rank; no cross-rank threading occurs.

### 5) Single-GPU (no MPI)

- Slurm: `--ntasks=1`, `--cpus-per-task=1-4`, `--gres=gpu:1`.
- Command:

```bash
srun --gres=gpu:1 python main.py --config "${CONFIG}" --device gpu
```

- Notes:
  - CuPy is required; if unavailable, the model logs a warning and falls back to CPU.
  - Keep inputs/outputs on shared storage unless you copy them to local SSD and back.

### 6) MPI + multi-GPU (one GPU per rank)

- Slurm example for 2 nodes, 2 GPUs each:
  - `--nodes=2`
  - `--ntasks-per-node=2`
  - `--gpus-per-task=1` (or `--gres=gpu:2` + `--ntasks-per-node=2`)
  - `--cpus-per-task=2` (adjust for host-side work)
- Command:

```bash
srun --mpi=pmix_v3 --gpu-bind=single:1 python main.py --config "${CONFIG}" --device gpu
```

- Notes:
  - Each MPI rank uses only the GPU assigned by Slurm (`--gpu-bind` or the launcher’s default binding).
  - MPI communication remains on CPU; GPU is used for CuPy-enabled runoff kernels.

### 7) MPI + GPU + shared memory (full hybrid)

- Slurm: combine GPU binding with CPU threads, e.g.,
  - `--nodes=1`
  - `--ntasks-per-node=2`
  - `--gpus-per-task=1`
  - `--cpus-per-task=8`
- Config:
  - `compute.device="gpu"`
  - `compute.shared_memory.enabled=true`
  - `compute.shared_memory.workers=8`
- Command:

```bash
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
srun --mpi=pmix_v3 --gpu-bind=single:1 python main.py --config "${CONFIG}" --device gpu
```

---

## Data, restarts, and outputs on HPC

- **Inputs**: Place domain and rain NetCDFs on shared storage. Rank 0 reads inputs and scatters to other ranks; ensure low-latency access for rank 0 (e.g., node-local copy if beneficial).
- **Outputs**: Rank 0 writes NetCDF outputs and restarts. Point `--out-nc` and `restart.out` to a shared path or to node-local storage with a final copy back.
- **Restarts**: To resume a job that was preempted or stopped, submit a new job with `--restart-in /path/to/restart_state.nc`. Restart files are written by rank 0 according to `restart.every`.

---

## Performance and troubleshooting tips

- Match `--cpus-per-task` to `compute.shared_memory.workers` to avoid oversubscription.
- For GPU jobs, keep `CUDA_VISIBLE_DEVICES` managed by Slurm (`--gpu-bind` or `--gpus-per-task`). Avoid manual assignment unless debugging.
- Use **local scratch** for temporary I/O if your cluster provides fast node-local SSDs; copy results back before job end.
- Prefer `srun` over `mpirun` inside Slurm allocations to inherit task affinity and GPU binding cleanly.
- Verify CuPy picks up the correct CUDA runtime by running a short smoke test on a GPU node:

```bash
srun --gres=gpu:1 python - <<'PY'
import cupy
print("CuPy device:", cupy.cuda.runtime.getDevice())
print("CuPy version:", cupy.__version__)
PY
```

- Check logs (`--output` path) for:
  - Warnings about missing GPUs → falls back to CPU.
  - Particle migration counts → indicates MPI communication is occurring.
  - Restart writes → confirm long jobs are checkpointing as expected.
