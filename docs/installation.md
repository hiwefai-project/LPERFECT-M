# Installation

This guide walks through installing **LPERFECT** with either a Python virtual environment + `pip` or with **Conda**. Each section is tailored to common environments so you can choose the path that matches your workflow.

- For **development**, see the PyCharm and Visual Studio Code sections.
- For **production**, see the workstation/Linux server and HPC sections (modulefile + Slurm, or Snap + Slurm).
- GPU acceleration is optional; CPU-only installs work everywhere.

## Prerequisites (all methods)

- A supported Python version (see `requirements.txt`).
- NetCDF tooling (GDAL/CDO/NCO) if you plan to prepare inputs locally.
- For MPI runs: an MPI runtime such as OpenMPI or MPICH plus `mpi4py` (installed via the commands below).
- For GPU runs: CUDA drivers and a compatible CuPy build (`cupy-cuda11x` or `cupy-cuda12x`).

---

## Method A — Python virtual environment + pip

### 1) Development in PyCharm (CPU or GPU)

1. **Create the environment**
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt
   # Optional GPU extras
   pip install "cupy-cuda12x"  # pick the CUDA build that matches your drivers
   # Optional MPI runtime if not provided by the OS
   # On Debian/Ubuntu: sudo apt-get install libopenmpi-dev
   pip install mpi4py
   ```
2. **Add the interpreter in PyCharm**
   - `File` → `Settings` → `Project: LPERFECT` → `Python Interpreter` → `Add Interpreter` → `Existing environment`.
   - Browse to `.venv/bin/python`.
3. **Set run configuration**
   - `Run` → `Edit Configurations...`.
   - Script path: `main.py`.
   - Parameters (example): `--config config.json --device cpu`.
4. **GPU or MPI runs in PyCharm (optional)**
   - GPU: add `--device gpu` to the parameters after installing CuPy.
   - MPI: configure a **Compound** run with `mpirun` or use a PyCharm external tool: `mpirun -np 4 ${PYTHON} main.py --config config.json`.

### 2) Development in Visual Studio Code (CPU or GPU)

1. **Create the environment** (same as PyCharm):
   ```bash
   python -m venv .venv
   source .venv/bin/activate
   pip install --upgrade pip
   pip install -r requirements.txt
   pip install "cupy-cuda12x"  # optional GPU
   pip install mpi4py           # if you have an MPI runtime installed
   ```
2. **Select interpreter**
   - Open the Command Palette → `Python: Select Interpreter` → pick `.venv`.
3. **Configure run/debug**
   - Create `.vscode/launch.json` with a `program` pointing to `main.py` and args such as `--config config.json`.
4. **GPU or MPI runs**
   - GPU: add `--device gpu` in `args`.
   - MPI: run from the integrated terminal: `mpirun -np 4 .venv/bin/python main.py --config config.json`.

### 3) Production on a workstation/Linux server (multi-core, with or without GPU)

1. **Install system dependencies**
   ```bash
   sudo apt-get update
   sudo apt-get install -y python3-venv build-essential libopenmpi-dev  # adjust for your distro
   ```
2. **Create an isolated environment**
   ```bash
   python -m venv /opt/lperfect-env
   source /opt/lperfect-env/bin/activate
   pip install --upgrade pip
   pip install -r /path/to/LPERFECT/requirements.txt
   ```
3. **Enable multi-core / GPU**
   - Multi-core per rank is controlled in the config (`compute.shared_memory.*`).
   - GPU: `pip install "cupy-cuda12x"` and run with `--device gpu`.
   - MPI across cores or nodes: `mpirun -np 16 python main.py --config config.json --mpi-decomposition balanced`.
4. **Service-style runs**
   - For recurring jobs, wrap the activation + run command in a systemd unit or cron job.

### 4) Production on an HPC cluster (modulefile + Slurm)

This assumes your site provides compiler/MPI/Python modules.

1. **Load modules**
   ```bash
   module load python/3.10
   module load openmpi/4.1
   module load cudatoolkit/12.2   # only if using GPUs
   ```
2. **Create a per-project environment**
   ```bash
   python -m venv $HOME/lperfect-env
   source $HOME/lperfect-env/bin/activate
   pip install --upgrade pip
   pip install -r $HOME/LPERFECT/requirements.txt
   pip install mpi4py
   pip install "cupy-cuda12x"  # optional GPU
   ```
3. **Submit with Slurm**
   - Example CPU job script:
     ```bash
     #!/bin/bash
     #SBATCH -J lperfect_cpu
     #SBATCH -N 1
     #SBATCH -n 32
     #SBATCH -t 02:00:00
     #SBATCH -p compute
     module load python/3.10 openmpi/4.1
     source $HOME/lperfect-env/bin/activate
     srun -n 32 python main.py --config config.json --mpi-decomposition balanced
     ```
   - Example GPU job script (1 node, 4 GPUs):
     ```bash
     #!/bin/bash
     #SBATCH -J lperfect_gpu
     #SBATCH -N 1
     #SBATCH --gres=gpu:4
     #SBATCH -n 4
     #SBATCH -t 02:00:00
     #SBATCH -p gpu
     module load python/3.10 openmpi/4.1 cudatoolkit/12.2
     source $HOME/lperfect-env/bin/activate
     srun -n 4 python main.py --config config.json --device gpu --mpi-decomposition balanced
     ```

### 5) Production on an HPC cluster (Snap + Slurm)

Some clusters provide user-level **Snap** for isolated software stacks.

1. **Install Snap (if needed)**
   - Follow site guidance or install to your home directory. Ensure `snap` is on your `PATH` and permitted on compute nodes.
2. **Create a Snap-based environment**
   ```bash
   snap install --devmode --dangerous python-venv.snap   # or site-provided snap
   python -m venv $HOME/lperfect-snap-env
   source $HOME/lperfect-snap-env/bin/activate
   pip install --upgrade pip
   pip install -r $HOME/LPERFECT/requirements.txt
   pip install mpi4py
   pip install "cupy-cuda12x"  # optional GPU if the snap exposes CUDA
   ```
3. **Submit with Slurm**
   - Ensure the `snap` mount and environment activation are available inside the job.
   - Example:
     ```bash
     #!/bin/bash
     #SBATCH -J lperfect_snap
     #SBATCH -N 1
     #SBATCH -n 32
     #SBATCH -t 02:00:00
     module load slurm  # keep site defaults
     source $HOME/lperfect-snap-env/bin/activate
     srun -n 32 python main.py --config config.json --mpi-decomposition balanced
     ```

---

## Method B — Conda

Use Conda when you prefer managed binaries or cannot install system packages.

### 1) Development in PyCharm

1. **Create environment**
   ```bash
   conda create -n lperfect-dev python=3.10
   conda activate lperfect-dev
   conda install mpi4py
   pip install -r requirements.txt
   # Optional GPU
   conda install -c conda-forge cupy
   ```
2. **Select interpreter** in PyCharm: add the Conda env path (e.g., `~/miniconda3/envs/lperfect-dev/bin/python`).
3. **Run configurations**: same parameters as the pip workflow; MPI via `mpirun -np 4 python main.py ...`.

### 2) Development in Visual Studio Code

1. **Create environment** (same as above) and activate it in the terminal.
2. **Interpreter**: `Python: Select Interpreter` → pick `lperfect-dev`.
3. **Run/Debug**: set `program` to `main.py`; add arguments in `launch.json`.
4. **MPI**: run via terminal `mpirun -np 4 python main.py --config config.json`.

### 3) Production on a workstation/Linux server

1. **Create environment**
   ```bash
   conda create -y -n lperfect-prod python=3.10 mpi4py
   conda activate lperfect-prod
   pip install -r /path/to/LPERFECT/requirements.txt
   conda install -c conda-forge cupy  # optional GPU
   ```
2. **Run**
   ```bash
   mpirun -np 16 python main.py --config config.json --mpi-decomposition balanced
   # or GPU per rank
   mpirun -np 4 python main.py --config config.json --device gpu
   ```

### 4) Production on an HPC cluster (modulefile + Slurm)

1. **Load minimal modules** (if needed for Conda to work with site MPI)
   ```bash
   module load miniconda
   module load openmpi/4.1
   module load cudatoolkit/12.2   # optional
   ```
2. **Create environment**
   ```bash
   conda create -y -n lperfect-hpc python=3.10 mpi4py
   conda activate lperfect-hpc
   pip install -r $HOME/LPERFECT/requirements.txt
   conda install -c conda-forge cupy  # optional GPU
   ```
3. **Submit with Slurm**
   ```bash
   #!/bin/bash
   #SBATCH -J lperfect_conda
   #SBATCH -N 1
   #SBATCH -n 32
   #SBATCH -t 02:00:00
   module load miniconda openmpi/4.1
   source $(conda info --base)/etc/profile.d/conda.sh
   conda activate lperfect-hpc
   srun -n 32 python main.py --config config.json --mpi-decomposition balanced
   ```

### 5) Production on an HPC cluster (Snap + Slurm)

If Snap is preferred but Conda is allowed inside it, combine both:

1. **Set up Snap and Conda** (site instructions may provide a Snap with Conda baked in).
2. **Create environment** inside the Snap context:
   ```bash
   conda create -y -n lperfect-snap python=3.10 mpi4py
   conda activate lperfect-snap
   pip install -r $HOME/LPERFECT/requirements.txt
   conda install -c conda-forge cupy  # optional GPU
   ```
3. **Submit with Slurm**
   ```bash
   #!/bin/bash
   #SBATCH -J lperfect_snap_conda
   #SBATCH -N 1
   #SBATCH -n 32
   #SBATCH -t 02:00:00
   source $(conda info --base)/etc/profile.d/conda.sh
   conda activate lperfect-snap
   srun -n 32 python main.py --config config.json --mpi-decomposition balanced
   ```

---

## Quick reference

- **GPU support**: install a matching CuPy build and use `--device gpu`.
- **MPI**: install `mpi4py` and run with `mpirun -np <ranks>`.
- **Shared-memory threads**: configure `compute.shared_memory` in `config.json` (per rank, even under MPI).
- **Isolation**: prefer per-project envs (virtualenv or Conda) to keep dependencies reproducible across workstations and clusters.
