# -*- coding: utf-8 -*-
"""Core simulation driver for LPERFECT."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import typing primitives.
from typing import Any, Dict, Iterable, List, Optional  # import typing import Any, Dict, Iterable, List, Optional
from collections import Counter, defaultdict  # import collections import Counter, defaultdict
from time import perf_counter  # import time import perf_counter
# Import dataclasses.
from dataclasses import dataclass  # import dataclasses import dataclass

# Import logging.
import logging  # import logging

# Import numpy.
import numpy as np  # import numpy as np

# Import datetime helpers.
from datetime import timedelta  # import datetime import timedelta
from pathlib import Path  # import pathlib import Path
import os  # import os
import json  # import json

# Import local modules.
from .time_utils import datetime_to_hours_since_1900, parse_iso8601_to_datetime64, parse_iso8601_to_utc_datetime  # import .time_utils import datetime_to_hours_since_1900, parse_iso8601_to_datetime64, parse_iso8601_to_utc_datetime
from .d8 import build_downstream_index  # import .d8 import build_downstream_index
from .rain import build_rain_sources, blended_rain_step_mm_rank0  # import .rain import build_rain_sources, blended_rain_step_mm_rank0
from .runoff import precompute_scs_cn_params, scs_cn_cumulative_runoff_mm  # import .runoff import precompute_scs_cn_params, scs_cn_cumulative_runoff_mm
from .compute_backend import gpu_available, normalize_device  # import .compute_backend import gpu_available, normalize_device
from .hydraulics import spawn_particles_from_runoff_slab, advect_particles_one_step, local_volgrid_from_particles_slab  # import .hydraulics import spawn_particles_from_runoff_slab, advect_particles_one_step, local_volgrid_from_particles_slab
from .particles import Particles, empty_particles, concat_particles  # import .particles import Particles, empty_particles, concat_particles
from .risk import compute_flow_accum_area_m2, compute_risk_index  # import .risk import compute_flow_accum_area_m2, compute_risk_index
from .io_netcdf import write_results_netcdf_rank0, save_restart_netcdf_rank0, load_restart_netcdf_rank0  # import .io_netcdf import write_results_netcdf_rank0, save_restart_netcdf_rank0, load_restart_netcdf_rank0
from .shared_memory import SharedMemoryConfig  # import .shared_memory import SharedMemoryConfig

# Import MPI helpers.
from .mpi_utils import (  # import .mpi_utils import (
    HAVE_MPI,  # execute statement
    MPI,  # execute statement
    MPIConfig,  # execute statement
    PartitionPlan,  # execute statement
    build_partition_from_active,  # execute statement
    build_partition_from_weights,  # execute statement
    scatter_field_partition,  # execute statement
    gather_field_partition_to_rank0,  # execute statement
    gather_particles_to_rank0,  # execute statement
    scatter_particles_from_rank0,  # execute statement
    migrate_particles_partition,  # execute statement
    scatter_new_particles_even,  # execute statement
    rebalance_particles_even,  # execute statement
)  # execute statement

# Import Domain class.
from .domain import Domain  # import .domain import Domain

# Create a logger for this module.
logger = logging.getLogger("lperfect")  # set logger


def run_simulation(
    comm: Any,
    rank: int,
    size: int,
    cfg: Dict[str, Any],
    dom: Domain,
    domain_label: str | None = None,
    mpi_world_size: int | None = None,
) -> None:  # define function run_simulation
    """Run the full LPERFECT simulation (serial or MPI)."""  # execute statement

    # ------------------------------
    # Model parameters
    # ------------------------------
    mcfg = cfg["model"]  # set mcfg
    T_s = float(mcfg["T_s"])                                  # Total simulation time.  # set T_s
    dt_s = float(mcfg["dt_s"])                                # Timestep.  # set dt_s
    encoding = str(mcfg["encoding"])                           # D8 encoding.  # set encoding
    ia_ratio = float(mcfg["ia_ratio"])                         # Initial abstraction ratio.  # set ia_ratio
    particle_vol_m3 = float(mcfg["particle_vol_m3"])           # Target particle volume.  # set particle_vol_m3
    travel_time_mode = str(mcfg.get("travel_time_mode", "fixed")).lower().strip()  # Travel time mode.  # set travel_time_mode
    travel_time_s_cfg = float(mcfg["travel_time_s"])            # Hillslope hop time (fixed).  # set travel_time_s_cfg
    travel_time_channel_s_cfg = float(mcfg["travel_time_channel_s"])  # Channel hop time (fixed).  # set travel_time_channel_s_cfg
    travel_time_auto_cfg = mcfg.get("travel_time_auto", {})     # Auto travel-time config.  # set travel_time_auto_cfg
    outflow_sink = bool(mcfg["outflow_sink"])                  # Drop particles leaving domain.  # set outflow_sink
    log_every = int(mcfg.get("log_every", 10))                 # Diagnostics frequency.  # set log_every

    # Resolve compute device.
    compute_cfg = cfg.get("compute", {})  # set compute_cfg
    device = normalize_device(compute_cfg.get("device", "cpu"))  # set device
    parallel_cfg = compute_cfg.get("parallelization", {}) if isinstance(compute_cfg, dict) else {}  # set parallel_cfg
    parallel_schema = str(parallel_cfg.get("schema", "slab") or "slab").lower().strip()  # set parallel_schema
    if parallel_schema not in {"slab", "particles"}:  # check condition parallel_schema not supported:
        raise ValueError("compute.parallelization.schema must be 'slab' or 'particles'.")  # raise ValueError
    io_mode = str(parallel_cfg.get("io", "rank0") or "rank0").lower().strip()  # set io_mode
    if io_mode not in {"rank0", "all"}:  # check condition io_mode invalid:
        raise ValueError("compute.parallelization.io must be 'rank0' or 'all'.")  # raise ValueError
    io_rank0_only = io_mode == "rank0"  # set io_rank0_only
    particle_parallel = parallel_schema == "particles"  # set particle_parallel
    world_size = int(mpi_world_size) if mpi_world_size is not None else int(size)
    mpi_cfg = MPIConfig.from_dict(compute_cfg.get("mpi", {}), world_size=world_size)
    mpi_active = mpi_cfg.enabled and size > 1 and comm is not None
    if device == "gpu" and not gpu_available():  # check condition device == "gpu" and not gpu_available():
        if rank == 0:  # check condition rank == 0:
            logger.warning("GPU requested but CuPy not available; falling back to CPU.")  # execute statement
        device = "cpu"  # set device
    # Domain shape.
    nrows, ncols = dom.dem.shape  # set nrows, ncols

    # Optionally prune ranks when slabs would be empty or violate min_rows_per_rank.
    active_rows = int(np.count_nonzero(np.any(dom.active_mask, axis=1)))
    target_size = size
    if mpi_active and comm is not None:
        max_ranks_by_rows = max(1, active_rows // max(1, mpi_cfg.min_rows_per_rank))
        target_size = max(1, min(size, max_ranks_by_rows))
        if target_size < size:
            color = 0 if rank < target_size else MPI.UNDEFINED
            new_comm = comm.Split(color=color, key=rank)
            if new_comm == MPI.COMM_NULL:
                return
            comm = new_comm
            rank = comm.Get_rank()
            size = comm.Get_size()
    mpi_active = mpi_cfg.enabled and size > 1 and comm is not None

    if rank == 0:  # check condition rank == 0:
        logger.info("Compute device: %s", device)  # execute statement
        logger.info(
            "Distributed memory: %s (active ranks=%d, world_size=%d, decomposition=%s, min_rows_per_rank=%d, active_rows=%d)",
            "enabled" if mpi_active else "disabled",
            size if mpi_active else 1,
            world_size,
            mpi_cfg.decomposition,
            mpi_cfg.min_rows_per_rank,
            active_rows,
        )
        if mpi_active:
            logger.info(
                "Parallelization schema=%s, io_mode=%s (rank0_only=%s)",
                parallel_schema,
                io_mode,
                str(io_rank0_only),
            )
    shared_cfg = SharedMemoryConfig.from_dict(compute_cfg.get("shared_memory", {}))  # set shared_cfg
    if rank == 0 and shared_cfg.enabled and shared_cfg.workers > 1:  # check condition rank == 0 and shared_cfg.enabled and shared_cfg.workers > 1:
        logger.info(
            "Shared-memory parallelism: %s with %d workers (chunk=%d, min_particles=%d)",
            "auto-enabled" if shared_cfg.auto_enabled else "enabled",
            shared_cfg.workers,
            shared_cfg.chunk_size,
            shared_cfg.min_particles_per_worker,
        )  # execute statement

    metrics_root = cfg.get("metrics", {}) if isinstance(cfg.get("metrics", {}), dict) else {}  # set metrics_root
    metrics_cfg = metrics_root.get("parallelization", {}) if isinstance(metrics_root, dict) else {}  # set metrics_cfg
    metrics_enabled = bool(metrics_cfg.get("enabled", False))  # set metrics_enabled
    metrics_output = metrics_cfg.get("output", None)  # set metrics_output
    metrics_max_samples = int(metrics_cfg.get("max_samples", 256) or 0)  # set metrics_max_samples
    metrics_format = str(metrics_cfg.get("format", "detailed") or "detailed")  # set metrics_format
    if metrics_format not in {"detailed", "compact"}:  # check condition metrics_format not in {"detailed", "compact"}:
        raise ValueError("metrics.parallelization.format must be 'detailed' or 'compact'.")  # raise ValueError
    metrics_indent = 2 if metrics_format == "detailed" else None  # set metrics_indent
    metrics_separators = None if metrics_format == "detailed" else (",", ":")  # set metrics_separators
    assistant_cfg = metrics_root.get("assistant", {}) if isinstance(metrics_root, dict) else {}  # set assistant_cfg
    assistant_metrics_enabled = bool(assistant_cfg.get("enabled", False))  # set assistant_metrics_enabled
    assistant_metrics_output = assistant_cfg.get("output", None)  # set assistant_metrics_output
    assistant_metrics_format = str(assistant_cfg.get("format", "detailed") or "detailed")  # set assistant_metrics_format
    if assistant_metrics_format not in {"detailed", "compact"}:  # check condition assistant_metrics_format not valid:
        raise ValueError("metrics.assistant.format must be 'detailed' or 'compact'.")  # raise ValueError
    assistant_indent = 2 if assistant_metrics_format == "detailed" else None  # set assistant_indent
    assistant_separators = None if assistant_metrics_format == "detailed" else (",", ":")  # set assistant_separators
    balance_cfg = compute_cfg.get("mpi", {}).get("balance", {}) if isinstance(compute_cfg.get("mpi", {}), dict) else {}
    balance_every_steps = int(balance_cfg.get("every_steps", 0) or 0)
    balance_every_sim_s = float(balance_cfg.get("every_sim_s", 0.0) or 0.0)
    balance_auto = bool(balance_cfg.get("auto", False))
    balance_threshold = float(balance_cfg.get("imbalance_threshold", 2.0) or 2.0)
    next_balance_time: float | None = balance_every_sim_s if balance_every_sim_s > 0.0 else None

    # Start time for time-aware rain selection and output timestamps.
    start_time = parse_iso8601_to_datetime64(mcfg.get("start_time", None))  # set start_time
    start_datetime = parse_iso8601_to_utc_datetime(mcfg.get("start_time", None))  # set start_datetime

    # Precompute downstream lookup (replicated).
    valid, ds_r, ds_c = build_downstream_index(dom.d8, encoding)  # set valid, ds_r, ds_c
    # Static active-cell weights (reused when recomputing slab partitions).
    active_row_weights = np.sum(dom.active_mask.astype(np.int32), axis=1).astype(np.int64)

    # Build a rank-aware partition plan (contiguous row slabs, load-balanced by active cells).
    if particle_parallel and mpi_active:
        partition = PartitionPlan(counts=np.array([nrows], dtype=np.int32), starts=np.array([0], dtype=np.int32))
        r0, r1 = 0, nrows
    else:
        partition = build_partition_from_active(
            active_mask=dom.active_mask,
            size=max(1, size),
            strategy=mpi_cfg.decomposition,
            min_rows_per_rank=mpi_cfg.min_rows_per_rank,
        )
        r0, r1 = partition.bounds(rank if size > 1 else 0)

    def _rebalance_particles(reason: str, step_idx: int, step_time_s: float) -> None:
        """Gather state, rebuild slabs with particle-aware weights, and resplit."""
        nonlocal partition, P_slab, Q_slab, particles, r0, r1, CN_slab, next_balance_time
        if not mpi_active or comm is None:
            return

        if particle_parallel:
            counts_before = comm.allgather(int(particles.r.size))
            if rank == 0:
                logger.info(
                    "Particle-parallel rebalance (%s) step=%d time=%.1fs counts=%s",
                    reason,
                    step_idx + 1,
                    step_time_s,
                    counts_before,
                )
            particles = rebalance_particles_even(comm, particles)
            counts_after = comm.allgather(int(particles.r.size))
            if rank == 0:
                logger.info(
                    "Particle-parallel rebalance completed: %s -> %s",
                    counts_before,
                    counts_after,
                )
            return

        # Log current distribution.
        local_count = int(particles.r.size)
        counts_by_rank = comm.allgather(local_count)
        if rank == 0:
            logger.info(
                "MPI load-balance (%s) at step=%d time=%.1fs: particle_counts=%s",
                reason,
                step_idx + 1,
                step_time_s,
                counts_by_rank,
            )

        # Gather slabs and particles to rank0.
        local_particle_rows = (
            np.bincount(particles.r.astype(np.int32), minlength=nrows) if particles.r.size > 0 else np.zeros(nrows, dtype=np.int64)
        )  # local per-row particle counts
        particle_row_counts = np.zeros_like(local_particle_rows)  # buffer for the global row histogram
        comm.Allreduce(local_particle_rows, particle_row_counts, op=MPI.SUM)  # sum histograms across ranks
        row_weights = active_row_weights + particle_row_counts  # combine static active rows with live particle density
        row_weights = np.where(row_weights <= 0, 1, row_weights)  # enforce positive weights for every row

        # Build a new partition weighted by active cells and current particle density.
        if rank == 0:
            new_plan = build_partition_from_weights(
                row_weights=row_weights,
                size=max(1, size),
                strategy=mpi_cfg.decomposition,
                min_rows_per_rank=mpi_cfg.min_rows_per_rank,
            )
            counts = new_plan.counts
            starts = new_plan.starts
        else:
            counts = None
            starts = None

        current_plan = partition
        counts = comm.bcast(counts, root=0)  # broadcast candidate counts
        starts = comm.bcast(starts, root=0)  # broadcast candidate starts
        counts_arr = np.asarray(counts, dtype=np.int32)
        starts_arr = np.asarray(starts, dtype=np.int32)
        if np.array_equal(counts_arr, current_plan.counts) and np.array_equal(starts_arr, current_plan.starts):
            if rank == 0:
                logger.debug("MPI load-balance skipped: partition unchanged (reason=%s)", reason)
            return
        new_partition = PartitionPlan(counts=counts_arr, starts=starts_arr)

        # Gather slabs and particles to rank0 only after we know a redistribution is required.
        P_full = gather_field_partition_to_rank0(comm, current_plan, P_slab, nrows, ncols)
        Q_full = gather_field_partition_to_rank0(comm, current_plan, Q_slab, nrows, ncols)
        particles_all = gather_particles_to_rank0(comm, particles)

        # Scatter to the refreshed layout.
        P_slab = scatter_field_partition(comm, new_partition, P_full, nrows, ncols, np.float32)
        Q_slab = scatter_field_partition(comm, new_partition, Q_full, nrows, ncols, np.float32)
        particles = scatter_particles_from_rank0(comm, new_partition, particles_all)
        partition = new_partition
        r0, r1 = partition.bounds(rank if size > 1 else 0)
        CN_slab = dom.cn[r0:r1, :]
        new_counts_by_rank = comm.allgather(int(particles.r.size))
        if rank == 0:
            logger.info(
                "MPI load-balance completed: %s -> %s (reason=%s step=%d time=%.1fs)",
                counts_by_rank,
                new_counts_by_rank,
                reason,
                step_idx + 1,
                step_time_s,
            )

    def _maybe_rebalance(step_idx: int, step_time_s: float) -> None:
        """Trigger periodic or automatic particle-aware rebalancing."""
        nonlocal next_balance_time
        if not mpi_active or comm is None:
            return

        trigger = False
        reason = ""

        if balance_every_steps > 0 and ((step_idx + 1) % balance_every_steps == 0):
            trigger = True
            reason = "interval_steps"

        if balance_every_sim_s > 0.0 and next_balance_time is not None and (step_time_s + 1e-9) >= next_balance_time:
            trigger = True
            reason = reason or "interval_time"
            while next_balance_time is not None and (step_time_s + 1e-9) >= next_balance_time:
                next_balance_time += balance_every_sim_s

        if balance_auto:
            local_count = int(particles.r.size)
            counts_by_rank = comm.allgather(local_count)
            max_c = max(counts_by_rank) if counts_by_rank else 0
            min_c = min(counts_by_rank) if counts_by_rank else 0
            if max_c > 0 and (min_c == 0 or (max_c / max(1, min_c)) >= balance_threshold):
                trigger = True
                reason = reason or "auto_imbalance"

        if trigger:
            _rebalance_particles(reason or "manual", step_idx, step_time_s)

    def _compute_hop_distances(cell_area_m2: float | np.ndarray) -> np.ndarray:
        """Approximate hop distances (m) using cell area and downstream orientation."""
        nrows_loc, ncols_loc = valid.shape
        rr, cc = np.indices((nrows_loc, ncols_loc))
        dr = ds_r - rr
        dc = ds_c - cc
        diag = (np.abs(dr) == 1) & (np.abs(dc) == 1) & valid
        orth = valid & ~diag

        if np.isscalar(cell_area_m2):
            base = float(np.sqrt(cell_area_m2))
            dist = np.zeros_like(valid, dtype=np.float32)
            dist[orth] = base
            dist[diag] = base * np.sqrt(2.0)
            return dist

        base = np.sqrt(cell_area_m2.astype(np.float32))
        dist = np.zeros_like(base, dtype=np.float32)
        dist[orth] = base[orth]
        dist[diag] = base[diag] * np.sqrt(2.0)
        return dist

    def _compute_auto_travel_times() -> tuple[np.ndarray, np.ndarray]:
        """Derive hop travel times from velocities and geometry."""
        vel_hill = float(travel_time_auto_cfg.get("hillslope_velocity_ms", 0.5))
        vel_ch = float(travel_time_auto_cfg.get("channel_velocity_ms", 1.5))
        min_s = float(travel_time_auto_cfg.get("min_s", 0.25))
        max_s = float(travel_time_auto_cfg.get("max_s", 3600.0))

        if vel_hill <= 0.0 or vel_ch <= 0.0:
            raise ValueError("travel_time_auto velocities must be positive.")

        dist = _compute_hop_distances(dom.cell_area_m2)
        hill = np.clip(dist / vel_hill, min_s, max_s).astype(np.float32)
        ch = np.clip(dist / vel_ch, min_s, max_s).astype(np.float32)
        return hill, ch

    if travel_time_mode == "auto":
        travel_time_s, travel_time_channel_s = _compute_auto_travel_times()  # set travel_time_s
        if rank == 0:
            logger.info(
                "Travel time mode=auto: hillslope median=%.3fs channel median=%.3fs",
                float(np.median(travel_time_s[valid])),
                float(np.median(travel_time_channel_s[valid])),
            )
    elif travel_time_mode == "fixed":
        travel_time_s = travel_time_s_cfg  # set travel_time_s
        travel_time_channel_s = travel_time_channel_s_cfg  # set travel_time_channel_s
    else:
        raise ValueError(f"Unknown travel_time_mode '{travel_time_mode}'. Use 'fixed' or 'auto'.")

    # Rain sources configuration (replicated config, rank0 reads files).
    rain_sources = build_rain_sources(cfg)  # set rain_sources

    # Output configuration.
    output_cfg = cfg.get("output", {})  # set output_cfg
    out_nc = output_cfg.get("out_netcdf", None)  # set out_nc
    save_every_s = float(output_cfg.get("save_every_s", 0.0) or 0.0)  # set save_every_s
    rotate_every_s = float(output_cfg.get("rotate_every_s", 0.0) or 0.0)  # set rotate_every_s
    outflow_geojson = output_cfg.get("outflow_geojson", None)  # set outflow_geojson
    record_outflow_points = bool(outflow_geojson)  # set record_outflow_points
    if save_every_s < 0.0 or rotate_every_s < 0.0:  # check condition save_every_s < 0.0 or rotate_every_s < 0.0:
        raise ValueError("output.save_every_s and output.rotate_every_s must be non-negative.")  # raise ValueError("output.save_every_s and output.rotate_every_s must be non-negative.")
    if save_every_s > 0.0 and rotate_every_s > 0.0:  # check condition save_every_s > 0.0 and rotate_every_s > 0.0:
        if rank == 0:
            logger.warning(
                "Both output.save_every_s=%.0f and output.rotate_every_s=%.0f are set; "
                "rotating outputs take precedence.",
                save_every_s,
                rotate_every_s,
            )
        save_every_s = 0.0
    rotate_enabled = rotate_every_s > 0.0  # set rotate_enabled
    output_interval_s = rotate_every_s if rotate_enabled else save_every_s  # set output_interval_s
    if output_interval_s > 0.0 and not out_nc:  # check condition output_interval_s > 0.0 and not out_nc:
        raise ValueError("output.out_netcdf must be set to enable periodic outputs.")  # raise ValueError("output.out_netcdf must be set to enable periodic outputs.")
    output_file_initialized = False  # set output_file_initialized
    rotation_index = 0  # set rotation_index
    next_output_s: float | None = None  # set next_output_s
    inundation_threshold_m = float(output_cfg.get("inundation_threshold_m", 0.01))  # set inundation_threshold_m

    risk_cfg = cfg.get("risk", {})  # set risk_cfg
    do_risk = bool(risk_cfg.get("enabled", True))  # set do_risk
    risk_accum = None  # set risk_accum
    flood_depth_max = np.where(dom.active_mask, 0.0, np.nan).astype(np.float32)  # set flood_depth_max
    inundation_mask_max = np.where(dom.active_mask, 0, 0).astype(np.int8)  # set inundation_mask_max

    # ------------------------------
    # Initialize state (possibly from restart)
    # ------------------------------
    rst_cfg = cfg.get("restart", {})  # set rst_cfg
    restart_in = rst_cfg.get("in", None)  # set restart_in

    if size == 1:  # check condition size == 1:
        # Serial state uses full arrays.
        if restart_in:  # check condition restart_in:
            r = load_restart_netcdf_rank0(restart_in)  # set r
            P_slab = r["P_cum_mm"].astype(np.float32)  # set P_slab
            Q_slab = r["Q_cum_mm"].astype(np.float32)  # set Q_slab
            particles = Particles(r=r["r"], c=r["c"], vol=r["vol"], tau=r["tau"])  # set particles
            cum_rain = r["cum_rain_vol_m3"]  # set cum_rain
            cum_runoff = r["cum_runoff_vol_m3"]  # set cum_runoff
            cum_outflow = r["cum_outflow_vol_m3"]  # set cum_outflow
            elapsed_s0 = r["elapsed_s"]  # set elapsed_s0
        else:  # fallback branch
            P_slab = np.zeros((nrows, ncols), dtype=np.float32)  # set P_slab
            Q_slab = np.zeros((nrows, ncols), dtype=np.float32)  # set Q_slab
            particles = empty_particles()  # set particles
            cum_rain = cum_runoff = cum_outflow = 0.0  # set cum_rain
            elapsed_s0 = 0.0  # set elapsed_s0
        r0, r1 = 0, nrows  # set r0, r1
    else:  # fallback branch
        # MPI: rank0 loads restart and scatters state.
        if restart_in and rank == 0:  # check condition restart_in and rank == 0:
            r = load_restart_netcdf_rank0(restart_in)  # set r

            # Optional strict checks.
            if bool(rst_cfg.get("strict_grid_check", True)):  # check condition bool(rst_cfg.get("strict_grid_check", True)):
                if r["P_cum_mm"].shape != (nrows, ncols):  # check condition r["P_cum_mm"].shape != (nrows, ncols):
                    raise ValueError("Restart grid mismatch for P_cum_mm")  # raise ValueError("Restart grid mismatch for P_cum_mm")
                if r["Q_cum_mm"].shape != (nrows, ncols):  # check condition r["Q_cum_mm"].shape != (nrows, ncols):
                    raise ValueError("Restart grid mismatch for Q_cum_mm")                  # raise ValueError("Restart grid mismatch for Q_cum_mm")

            P_full = r["P_cum_mm"].astype(np.float32)  # set P_full
            Q_full = r["Q_cum_mm"].astype(np.float32)  # set Q_full
            particles_all = Particles(r=r["r"], c=r["c"], vol=r["vol"], tau=r["tau"])  # set particles_all
            cum_rain = r["cum_rain_vol_m3"]  # set cum_rain
            cum_runoff = r["cum_runoff_vol_m3"]  # set cum_runoff
            cum_outflow = r["cum_outflow_vol_m3"]  # set cum_outflow
            elapsed_s0 = r["elapsed_s"]  # set elapsed_s0
        else:  # fallback branch
            P_full = None  # set P_full
            Q_full = None  # set Q_full
            particles_all = None  # set particles_all
            cum_rain = cum_runoff = cum_outflow = 0.0  # set cum_rain
            elapsed_s0 = 0.0  # set elapsed_s0

        # Broadcast scalar diagnostics and elapsed time.
        scal = np.array([cum_rain, cum_runoff, cum_outflow, elapsed_s0], dtype=np.float64) if rank == 0 else np.zeros(4, dtype=np.float64)  # set scal
        comm.Bcast(scal, root=0)  # execute statement
        cum_rain, cum_runoff, cum_outflow, elapsed_s0 = map(float, scal.tolist())  # set cum_rain, cum_runoff, cum_outflow, elapsed_s0

        if particle_parallel:
            P_slab = comm.bcast(P_full if P_full is not None else np.zeros((nrows, ncols), dtype=np.float32), root=0).astype(np.float32)
            Q_slab = comm.bcast(Q_full if Q_full is not None else np.zeros((nrows, ncols), dtype=np.float32), root=0).astype(np.float32)
            particles_all = particles_all if particles_all is not None else empty_particles()  # normalize None to empty on rank0
            particles_all = particles_all if rank == 0 else empty_particles()  # keep non-root empty
            particles = rebalance_particles_even(comm, particles_all)  # distribute evenly
            r0, r1 = 0, nrows  # ensure full-range indexing
        else:
            # Scatter slab fields.
            P_slab = scatter_field_partition(comm, partition, P_full, nrows, ncols, np.float32)  # set P_slab
            Q_slab = scatter_field_partition(comm, partition, Q_full, nrows, ncols, np.float32)  # set Q_slab

            # Scatter particles by row ownership.
            particles = scatter_particles_from_rank0(comm, partition, particles_all)  # set particles

    if rank == 0:
        logger.debug("State arrays initialized with dtype P=%s Q=%s", P_slab.dtype, Q_slab.dtype)

    CN_slab = dom.cn[r0:r1, :]
    cn_params = precompute_scs_cn_params(CN_slab, ia_ratio=ia_ratio, device=device)

    # ------------------------------
    # Output helpers and bookkeeping
    # ------------------------------
    if rotate_enabled and rotate_every_s > 0.0:  # check condition rotate_enabled and rotate_every_s > 0.0:
        rotation_index = int(np.floor(elapsed_s0 / rotate_every_s))  # set rotation_index

    if output_interval_s > 0.0:  # check condition output_interval_s > 0.0:
        next_output_s = output_interval_s * (int(elapsed_s0 // output_interval_s) + 1)  # set next_output_s

    output_written_any = False  # set output_written_any
    output_written_on_final_step = False  # set output_written_on_final_step

    def _output_time_hours(elapsed_s: float) -> float:
        """Return simulation timestamp in CF time units."""
        return datetime_to_hours_since_1900(start_datetime + timedelta(seconds=float(elapsed_s)))

    def _rotated_out_path(idx: int) -> str:
        """Build rotated output file path."""
        base = Path(str(out_nc))
        suffix = base.suffix or ".nc"
        return str(base.with_name(f"{base.stem}_{idx:04d}{suffix}"))

    def _tag_rank(path: str) -> str:
        """Return path suffixed with rank when io_mode='all'."""
        base = Path(path)
        suffix = base.suffix or ".nc"
        return str(base.with_name(f"{base.stem}_rank{rank:04d}{suffix}"))

    def _prepare_output_target() -> tuple[str, str]:
        """Return (path, mode) for the next output write."""
        nonlocal output_file_initialized, rotation_index
        if rotate_enabled:
            idx = rotation_index
            path = _rotated_out_path(idx)
            rotation_index += 1
            if rank == 0:
                logger.info("Rotating output NetCDF to %s (index=%04d)", path, idx)
                if os.path.exists(path):
                    os.remove(path)
            return (path if io_rank0_only or size == 1 else _tag_rank(path)), "w"

        if not out_nc:
            raise ValueError("Output NetCDF path not set.")

        if not output_file_initialized:
            if rank == 0 and os.path.exists(out_nc):
                os.remove(out_nc)
            output_file_initialized = True
            path = str(out_nc)
        else:
            path = str(out_nc)
            return path if io_rank0_only or size == 1 else _tag_rank(path), "a"

        return path if io_rank0_only or size == 1 else _tag_rank(path), "w"

    def _flood_from_volume(volgrid: np.ndarray) -> np.ndarray:
        """Convert cell volumes to depth, masking inactive cells."""
        if np.isscalar(dom.cell_area_m2):
            flood_depth = volgrid / float(dom.cell_area_m2)
        else:
            flood_depth = np.divide(
                volgrid, dom.cell_area_m2, out=np.zeros_like(volgrid), where=(dom.cell_area_m2 > 0)
            )
        return np.where(dom.active_mask, flood_depth, np.nan).astype(np.float32)

    def _compute_risk_field(runoff_mm: np.ndarray, flood_depth: np.ndarray) -> np.ndarray:
        """Compute risk index if enabled, otherwise fill with NaNs."""
        nonlocal risk_accum
        if not do_risk:
            return np.full_like(flood_depth, np.nan, dtype=np.float32)
        if risk_accum is None:
            risk_accum = compute_flow_accum_area_m2(dom.d8, encoding, dom.cell_area_m2, dom.active_mask).astype(np.float32)
            if rank == 0:
                logger.debug("Computed flow accumulation field (dtype=%s)", risk_accum.dtype)
        return compute_risk_index(
            runoff_cum_mm=runoff_mm.astype(np.float32, copy=False),
            flow_accum_m2=risk_accum,
            active_mask=dom.active_mask,
            balance=float(risk_cfg.get("balance", 0.55)),
            p_low=float(risk_cfg.get("p_low", 5.0)),
            p_high=float(risk_cfg.get("p_high", 95.0)),
        ).astype(np.float32)

    def _gather_outputs_for_rank0() -> tuple[np.ndarray | None, np.ndarray | None]:
        """Gather flood depth and risk to rank0 (or return locally in serial)."""
        if size == 1:
            volgrid = local_volgrid_from_particles_slab(particles, r0=0, r1=nrows, ncols=ncols, shared_cfg=shared_cfg)
            flood_depth = _flood_from_volume(volgrid)
            risk_field = _compute_risk_field(Q_slab, flood_depth)
            return flood_depth, risk_field

        if particle_parallel:
            vol_local = local_volgrid_from_particles_slab(particles, r0=0, r1=nrows, ncols=ncols, shared_cfg=shared_cfg)
            if io_rank0_only:
                vol_full = comm.reduce(vol_local, op=MPI.SUM, root=0)
                if rank == 0:
                    flood_depth = _flood_from_volume(vol_full)
                    risk_field = _compute_risk_field(Q_slab, flood_depth)
                    return flood_depth, risk_field
                return None, None
            vol_full = comm.allreduce(vol_local, op=MPI.SUM)
            flood_depth = _flood_from_volume(vol_full)
            risk_field = _compute_risk_field(Q_slab, flood_depth)
            return flood_depth, risk_field

        volgrid_slab = local_volgrid_from_particles_slab(particles, r0=r0, r1=r1, ncols=ncols, shared_cfg=shared_cfg)
        vol_full = gather_field_partition_to_rank0(comm, partition, volgrid_slab, nrows, ncols)
        Q_full = gather_field_partition_to_rank0(comm, partition, Q_slab, nrows, ncols)
        if io_rank0_only:
            if rank == 0:
                flood_depth = _flood_from_volume(vol_full)
                risk_field = _compute_risk_field(Q_full, flood_depth)
                return flood_depth, risk_field
            return None, None

        if rank == 0:
            flood_depth = _flood_from_volume(vol_full)
            risk_field = _compute_risk_field(Q_full, flood_depth)
        else:
            flood_depth = None
            risk_field = None
        flood_depth = comm.bcast(flood_depth, root=0)
        risk_field = comm.bcast(risk_field, root=0)
        return flood_depth, risk_field

    def _write_output(elapsed_s: float) -> bool:
        """Compute and write outputs for the given elapsed time."""
        if not out_nc:
            return False
        out_path, mode = _prepare_output_target()
        time_hours = _output_time_hours(elapsed_s)
        flood_depth, risk_field = _gather_outputs_for_rank0()
        should_write = (io_rank0_only and rank == 0) or (not io_rank0_only)
        if should_write and flood_depth is not None and risk_field is not None:
            inundation_mask = np.where(
                dom.active_mask & np.isfinite(flood_depth) & (flood_depth >= inundation_threshold_m), 1, 0
            ).astype(np.int8, copy=False)
            nonlocal flood_depth_max, inundation_mask_max
            # Update the running maximum in place to avoid large temporary arrays.
            np.fmax(flood_depth_max, flood_depth, out=flood_depth_max)
            # Preserve inactive cells as NaN after the in-place update.
            flood_depth_max[~dom.active_mask] = np.nan
            # Track any cell that has ever been inundated without allocating a fresh array.
            np.bitwise_or(inundation_mask_max, inundation_mask, out=inundation_mask_max)

            logger.info(
                "Adding output time slice at t=%.3fs (%.3f h) to %s",
                float(elapsed_s),
                float(time_hours),
                out_path,
            )
            write_results_netcdf_rank0(
                out_path,
                cfg,
                dom,
                flood_depth,
                risk_field,
                inundation_mask,
                flood_depth_max,
                inundation_mask_max,
                time_hours=time_hours,
                mode=mode,
            )
        return True

    # ------------------------------
    # Quality and outflow tracking
    # ------------------------------
    system_vol = float(particles.vol.sum())
    if size > 1:
        system_vol = float(comm.allreduce(system_vol, op=MPI.SUM))
    mass_error_series: list[float] = []
    total_hops = 0
    total_spawned_particles = 0
    total_outflow_particles = 0
    outflow_interval_counts: dict[int, Counter[tuple[int, int]]] = defaultdict(Counter)

    def _interval_index(elapsed_s: float) -> int:
        """Return zero-based interval index for a given elapsed time."""
        if output_interval_s <= 0.0:
            return 0
        return int(np.floor(elapsed_s / output_interval_s))

    def _build_quality_report(
        final_elapsed_s: float,
        system_vol_final: float,
        mass_error_series: list[float],
        cum_rain: float,
        cum_runoff: float,
        cum_outflow: float,
        total_hops: int,
        total_spawned_particles: int,
        total_outflow_particles: int,
        active_particles: int,
        steps_run: int,
    ) -> Dict[str, Any]:
        """Assemble numerical + hydrological diagnostics."""
        expected = cum_runoff - cum_outflow
        final_mass_error = system_vol_final - expected
        series = np.asarray(mass_error_series if mass_error_series else [final_mass_error], dtype=np.float64)
        mass_abs_max = float(np.max(np.abs(series)))
        mass_abs_mean = float(np.mean(series))

        tol_base = max(cum_runoff, cum_rain, 1.0)
        numerical_ok = abs(final_mass_error) <= (1e-3 * tol_base + 1e-6)

        runoff_to_rain = float(cum_runoff / cum_rain) if cum_rain > 0.0 else np.nan
        hydro_ok = (cum_rain <= 0.0 and cum_runoff <= 1e-6) or (cum_runoff <= cum_rain + 1e-3 * tol_base)

        return {
            "runtime": {"elapsed_s": float(final_elapsed_s), "steps": int(steps_run)},
            "volume_balance": {
                "rain_m3": float(cum_rain),
                "runoff_m3": float(cum_runoff),
                "outflow_m3": float(cum_outflow),
                "stored_m3": float(system_vol_final),
                "expected_storage_m3": float(expected),
                "final_mass_error_m3": float(final_mass_error),
                "mass_error_abs_max_m3": mass_abs_max,
                "mass_error_abs_mean_m3": mass_abs_mean,
                "mass_error_pct_of_runoff": float(final_mass_error / cum_runoff * 100.0) if cum_runoff > 0.0 else np.nan,
            },
            "hydrology": {
                "runoff_to_rain_ratio": runoff_to_rain,
                "hydrology_ok": bool(hydro_ok),
                "notes": "runoff should not exceed rainfall volume",
            },
            "particles": {
                "spawned_total": int(total_spawned_particles),
                "outflow_total": int(total_outflow_particles),
                "active_particles": int(active_particles),
                "total_hops": int(total_hops),
            },
            "status": {"numerical_ok": bool(numerical_ok), "hydrological_ok": bool(hydro_ok)},
        }

    def _log_quality_report(report: Dict[str, Any]) -> None:
        """Emit a human-readable quality summary."""
        logger.info("Simulation quality report:")
        bal = report["volume_balance"]
        hyd = report["hydrology"]
        part = report["particles"]
        logger.info(
            "  Volume balance (m3): rain=%.3e runoff=%.3e outflow=%.3e stored=%.3e expected=%.3e",
            bal["rain_m3"],
            bal["runoff_m3"],
            bal["outflow_m3"],
            bal["stored_m3"],
            bal["expected_storage_m3"],
        )
        logger.info(
            "  Mass error: final=%.3e (%.3f%% of runoff) max|err|=%.3e mean|err|=%.3e",
            bal["final_mass_error_m3"],
            bal["mass_error_pct_of_runoff"],
            bal["mass_error_abs_max_m3"],
            bal["mass_error_abs_mean_m3"],
        )
        logger.info(
            "  Hydrology: runoff/rain=%.3f ok=%s (%s)",
            hyd["runoff_to_rain_ratio"],
            hyd["hydrology_ok"],
            hyd["notes"],
        )
        logger.info(
            "  Particles: spawned=%d outflow=%d active=%d total_hops=%d",
            part["spawned_total"],
            part["outflow_total"],
            part["active_particles"],
            part["total_hops"],
        )

    def _write_outflow_geojson(
        interval_counts: dict[int, Counter[tuple[int, int]]],
        final_elapsed_s: float,
    ) -> None:
        """Write GeoJSON with per-interval outflow particle counts."""
        if not outflow_geojson or rank != 0:
            return

        if output_interval_s > 0.0:
            start_interval = _interval_index(elapsed_s0)
            end_interval = _interval_index(final_elapsed_s)
            interval_base = min(start_interval, min(interval_counts.keys(), default=start_interval))
            n_intervals = max(1, end_interval - interval_base + 1)
        else:
            interval_base = 0
            n_intervals = 1

        per_cell: dict[tuple[int, int], list[int]] = {}
        for idx, counter in interval_counts.items():
            for (r, c), cnt in counter.items():
                arr = per_cell.setdefault((r, c), [0] * n_intervals)
                if idx - interval_base >= n_intervals:
                    arr.extend([0] * (idx - interval_base - n_intervals + 1))
                arr[idx - interval_base] += int(cnt)

        features = []
        for (r, c), counts in per_cell.items():
            if len(counts) < n_intervals:
                counts = counts + [0] * (n_intervals - len(counts))
            x = float(dom.x_vals[c])
            y = float(dom.y_vals[r])
            features.append(
                {
                    "type": "Feature",
                    "geometry": {"type": "Point", "coordinates": [x, y]},
                    "properties": {
                        "row": int(r),
                        "col": int(c),
                        "particles_per_interval": counts,
                        "total_particles": int(sum(counts)),
                        "interval_index_offset": int(interval_base),
                        "save_interval_s": float(output_interval_s),
                    },
                }
            )

        Path(outflow_geojson).parent.mkdir(parents=True, exist_ok=True)
        collection = {
            "type": "FeatureCollection",
            "features": features,
            "properties": {
                "description": "LPERFECT outflow impact points",
                "interval_index_offset": int(interval_base),
                "save_interval_s": float(output_interval_s),
                "elapsed_start_s": float(elapsed_s0),
                "elapsed_end_s": float(final_elapsed_s),
            },
        }
        with open(outflow_geojson, "w", encoding="utf-8") as f:
            json.dump(collection, f)
        logger.info("Outflow GeoJSON written: %s (features=%d)", outflow_geojson, len(features))
    # ------------------------------
    # Main time loop
    # ------------------------------
    run_wall_start = perf_counter()
    step_wall_samples: list[float] = []
    step_migration_ratio: list[float] = []
    step_hops_throughput: list[float] = []
    sampled_steps: list[dict[str, Any]] = []
    if elapsed_s0 > T_s + 1e-9:  # check condition elapsed_s0 > T_s + 1e-9:
        raise ValueError(f"Restart elapsed_s={elapsed_s0} exceeds total T_s={T_s}.")  # raise ValueError(f"Restart elapsed_s={elapsed_s0} exceeds total T_s={T_s}.")

    remaining_s = max(0.0, T_s - elapsed_s0)  # set remaining_s
    steps = int(np.ceil(remaining_s / dt_s))  # set steps

    if rank == 0:  # check condition rank == 0:
        logger.info("LPERFECT start: T_s=%.1f dt_s=%.1f steps=%d ranks=%d", T_s, dt_s, steps, size)  # execute statement

    sample_stride = max(1, int(np.ceil(steps / metrics_max_samples))) if (metrics_enabled and metrics_max_samples > 0) else 1

    for k in range(steps):  # loop over k in range(steps):
        step_t0 = perf_counter()
        runoff_t = 0.0
        advection_t = 0.0
        migration_t = 0.0
        active_particles_before_migration = int(particles.r.size)
        migrated_particles_global = 0
        # Elapsed time at end of step.
        step_time_s = elapsed_s0 + min((k + 1) * dt_s, remaining_s)  # set step_time_s

        # Absolute timestamp (for rain selection) if configured.
        sim_time = (start_time + np.timedelta64(int(round(step_time_s)), "s")) if start_time is not None else None  # set sim_time

        # Rank0 reads and blends rainfall, then scatters slabs.
        if rank == 0:  # check condition rank == 0:
            rain_full, rain_updated = blended_rain_step_mm_rank0(  # set rain_step_mm and rain_updated
                sources=rain_sources,  # set sources
                shape=(nrows, ncols),  # set shape
                dt_s=dt_s,  # set dt_s
                step_idx=k,  # set step_idx
                sim_time=sim_time,  # set sim_time
            )  # execute statement
            # Force inactive cells to zero rainfall.
            rain_full = np.where(dom.active_mask, rain_full, 0.0).astype(np.float32)  # set rain_step_mm
        else:  # fallback branch
            rain_full = None  # set rain_full
            rain_updated = False  # set rain_updated when not rank0

        if size > 1:  # check condition size > 1:
            if particle_parallel:
                rain_slab_mm = comm.bcast(rain_full, root=0).astype(np.float32, copy=False)  # set rain_slab_mm
            else:
                rain_slab_mm = scatter_field_partition(comm, partition, rain_full, nrows, ncols, np.float32)  # set rain_slab_mm
            rain_updated = bool(comm.bcast(rain_updated, root=0))  # broadcast rain_updated to all ranks
        else:
            rain_slab_mm = rain_full.astype(np.float32, copy=False)  # type: ignore[union-attr]

        # Compute rain volume injected this step (for diagnostics).
        if particle_parallel and size > 1 and rank != 0:
            rain_vol_step_local = 0.0  # avoid double-counting shared rain fields
        elif np.isscalar(dom.cell_area_m2):  # check condition np.isscalar(dom.cell_area_m2):
            rain_vol_step_local = float((rain_slab_mm / 1000.0).sum() * float(dom.cell_area_m2))  # set rain_vol_step_local
        else:  # fallback branch
            rain_vol_step_local = float(((rain_slab_mm / 1000.0) * dom.cell_area_m2[r0:r1, :]).sum())  # set rain_vol_step_local

        # Update cumulative precipitation.
        P_slab = P_slab + rain_slab_mm  # set P_slab

        # Compute cumulative runoff with CN and incremental runoff for this step.
        runoff_t0 = perf_counter()
        Q_cum_slab_full = scs_cn_cumulative_runoff_mm(
            P_slab,
            CN_slab,
            ia_ratio=ia_ratio,
            device=device,
            params=cn_params,
        )  # set Q_cum_slab
        Q_cum_slab = Q_cum_slab_full.astype(np.float32)  # set Q_cum_slab
        dQ_mm = np.maximum(Q_cum_slab - Q_slab, 0.0)  # set dQ_mm
        Q_slab = Q_cum_slab  # set Q_slab
        runoff_t = perf_counter() - runoff_t0

        # Convert incremental runoff to meters.
        runoff_depth_m = dQ_mm / 1000.0  # set runoff_depth_m

        # Spawn new particles from incremental runoff.
        if particle_parallel and size > 1 and rank != 0:
            newp = empty_particles()
            spawned_vol_local = 0.0
        else:
            newp, spawned_vol_local = spawn_particles_from_runoff_slab(  # set newp, spawned_vol_local
                runoff_depth_m_slab=runoff_depth_m,  # set runoff_depth_m_slab
                r0=r0,  # set r0
                cell_area_m2=dom.cell_area_m2,  # set cell_area_m2
                particle_vol_m3=particle_vol_m3,  # set particle_vol_m3
                active_mask_global=dom.active_mask,  # set active_mask_global
            )  # execute statement
        if particle_parallel and mpi_active:
            distributed_newp = scatter_new_particles_even(comm, newp)  # spread spawn evenly without touching existing layout
            particles = concat_particles(particles, distributed_newp)  # set particles
            new_particles_local = int(distributed_newp.r.size)
        else:
            particles = concat_particles(particles, newp)  # set particles
            new_particles_local = int(newp.r.size)
        new_particles = comm.allreduce(new_particles_local, op=MPI.SUM) if size > 1 else new_particles_local
        total_spawned_particles += new_particles

        # Optional load balancing of particles across ranks before advection.
        _maybe_rebalance(step_idx=k, step_time_s=step_time_s)

        # Advect particles.
        advect_t0 = perf_counter()
        particles, outflow_vol_local, nhops_local, outflow_points_local, outflow_particles_local = advect_particles_one_step(  # set particles, outflow_vol_local, nhops_local, outflow_points_local, outflow_particles_local
            particles=particles,  # set particles
            valid=valid,  # set valid
            ds_r=ds_r,  # set ds_r
            ds_c=ds_c,  # set ds_c
            dt_s=dt_s,  # set dt_s
            travel_time_s=travel_time_s,  # set travel_time_s
            travel_time_channel_s=travel_time_channel_s,  # set travel_time_channel_s
            channel_mask=dom.channel_mask,  # set channel_mask
            outflow_sink=outflow_sink,  # set outflow_sink
            shared_cfg=shared_cfg,  # set shared_cfg
            track_outflow_points=record_outflow_points,  # set track_outflow_points
            return_outflow_particles=True,  # set return_outflow_particles
        )  # execute statement
        advection_t = perf_counter() - advect_t0

        # Migrate particles between slabs (MPI).
        if size > 1 and not particle_parallel:  # check condition size > 1:
            if metrics_enabled:
                active_particles_before_migration = int(comm.allreduce(active_particles_before_migration, op=MPI.SUM))
            mig_t0 = perf_counter()
            particles, migrated_local = migrate_particles_partition(comm, particles, plan=partition)  # migrate and count
            migration_t = perf_counter() - mig_t0
            if metrics_enabled:
                migrated_particles_global = int(comm.allreduce(migrated_local, op=MPI.SUM))

        # Compute current system volume (sum of particle volumes).
        system_vol_local = float(particles.vol.sum())  # set system_vol_local
        outflow_particle_count_step = int(sum(outflow_points_local.values())) if record_outflow_points else 0

        # Reduce diagnostics to global totals.
        if size > 1:  # check condition size > 1:
            rain_vol_step = comm.allreduce(rain_vol_step_local, op=MPI.SUM)  # set rain_vol_step
            runoff_vol_step = comm.allreduce(spawned_vol_local, op=MPI.SUM)  # set runoff_vol_step
            outflow_vol_step = comm.allreduce(outflow_vol_local, op=MPI.SUM)  # set outflow_vol_step
            system_vol = comm.allreduce(system_vol_local, op=MPI.SUM)  # set system_vol
            nhops = comm.allreduce(nhops_local, op=MPI.SUM)  # set nhops
            outflow_particle_count_global = comm.allreduce(outflow_particle_count_step, op=MPI.SUM)
        else:  # fallback branch
            rain_vol_step = rain_vol_step_local  # set rain_vol_step
            runoff_vol_step = spawned_vol_local  # set runoff_vol_step
            outflow_vol_step = outflow_vol_local  # set outflow_vol_step
            system_vol = system_vol_local  # set system_vol
            nhops = nhops_local  # set nhops
            outflow_particle_count_global = outflow_particle_count_step

        # Update cumulative totals.
        cum_rain += rain_vol_step  # execute statement
        cum_runoff += runoff_vol_step  # execute statement
        cum_outflow += outflow_vol_step  # execute statement
        total_outflow_particles += int(outflow_particle_count_global)
        total_hops += int(nhops)

        # Mass balance error: stored volume minus (generated - outflow).
        expected = cum_runoff - cum_outflow  # set expected
        mass_error = system_vol - expected  # set mass_error
        mass_error_series.append(float(mass_error))

        if record_outflow_points and outflow_points_local:
            interval_idx = _interval_index(step_time_s)
            outflow_interval_counts[interval_idx].update(outflow_points_local)

        # Capture the wall-clock duration for this step (max across ranks when MPI is active).
        step_duration_local = perf_counter() - step_t0  # set step_duration_local
        step_duration = step_duration_local  # set step_duration
        if size > 1:  # check condition size > 1:
            step_duration = float(comm.allreduce(step_duration_local, op=MPI.MAX))  # set step_duration

        # Track how many particles are currently active across all ranks.
        active_particles_local_now = int(particles.r.size)  # set active_particles_local_now
        active_particles_global = active_particles_local_now  # set active_particles_global
        if size > 1:  # check condition size > 1:
            active_particles_global = int(comm.allreduce(active_particles_local_now, op=MPI.SUM))  # set active_particles_global

        particles_per_second = float(active_particles_global) / max(step_duration, 1e-9)  # set particles_per_second

        if metrics_enabled:
            step_wall_samples.append(step_duration)
            if size > 1:
                migration_ratio_pct = (
                    float(migrated_particles_global) / float(active_particles_before_migration + 1e-9) * 100.0
                    if active_particles_before_migration > 0
                    else 0.0
                )
            else:
                migration_ratio_pct = 0.0
            step_migration_ratio.append(migration_ratio_pct)
            throughput = float(nhops) / max(step_duration, 1e-9)
            step_hops_throughput.append(throughput)

            if metrics_max_samples <= 0 or (k % sample_stride == 0):
                sampled_steps.append(
                    {
                        "step": int(k + 1),
                        "sim_time_s": float(step_time_s),
                        "wall_clock_s": float(step_duration),
                        "segments": {
                            "runoff_spawn_s": float(runoff_t),
                            "advection_s": float(advection_t),
                            "migration_s": float(migration_t),
                        },
                        "particles": {
                            "active_before_migration": int(active_particles_before_migration),
                            "migrated": int(migrated_particles_global),
                            "outflow": int(outflow_particle_count_global),
                        },
                        "hops": int(nhops),
                        "throughput_hops_per_s": float(throughput),
                        "migration_ratio_pct": float(migration_ratio_pct),
                    }
                )

        # Periodic diagnostics (rank0 only).
        if rank == 0 and log_every > 0 and ((k + 1) % log_every == 0 or (k + 1) == steps):  # check condition rank == 0 and log_every > 0 and ((k + 1) % log_every == 0 or (k + 1) == steps):
            logger.info(  # execute statement
                "mass_balance step=%d time_s=%.1f rain_m3=%.3e runoff_m3=%.3e outflow_m3=%.3e system_m3=%.3e err_m3=%.3e hops=%d active_particles=%d particles_per_s=%.1f",  # execute statement
                (k + 1),
                step_time_s,
                rain_vol_step,
                runoff_vol_step,
                outflow_vol_step,
                system_vol,
                mass_error,
                nhops,
                active_particles_global,
                particles_per_second,
            )  # execute statement

        # Restart writing if enabled.
        rst_out = rst_cfg.get("out", None)  # set rst_out
        rst_every = int(rst_cfg.get("every", 0))  # set rst_every
        if rst_out and rst_every > 0 and (((k + 1) % rst_every == 0) or ((k + 1) == steps)):  # check condition rst_out and rst_every > 0 and (((k + 1) % rst_every == 0) or ((k + 1) == steps)):
            if size == 1:  # check condition size == 1:
                save_restart_netcdf_rank0(  # execute statement
                    out_path=str(rst_out),  # set out_path
                    cfg=cfg,  # set cfg
                    dom=dom,  # set dom
                    elapsed_s=float(step_time_s),  # set elapsed_s
                    cum_rain_vol_m3=float(cum_rain),  # set cum_rain_vol_m3
                    cum_runoff_vol_m3=float(cum_runoff),  # set cum_runoff_vol_m3
                    cum_outflow_vol_m3=float(cum_outflow),  # set cum_outflow_vol_m3
                    P_cum_mm_full=P_slab,  # set P_cum_mm_full
                    Q_cum_mm_full=Q_slab,  # set Q_cum_mm_full
                    particles_all=particles,  # set particles_all
                )  # execute statement
            else:  # fallback branch
                if particle_parallel:
                    particles_all = gather_particles_to_rank0(comm, particles)
                    if io_rank0_only and rank == 0:
                        save_restart_netcdf_rank0(
                            out_path=str(rst_out),
                            cfg=cfg,
                            dom=dom,
                            elapsed_s=float(step_time_s),
                            cum_rain_vol_m3=float(cum_rain),
                            cum_runoff_vol_m3=float(cum_runoff),
                            cum_outflow_vol_m3=float(cum_outflow),
                            P_cum_mm_full=P_slab,
                            Q_cum_mm_full=Q_slab,
                            particles_all=particles_all,
                        )
                    elif not io_rank0_only:
                        save_restart_netcdf_rank0(
                            out_path=_tag_rank(str(rst_out)),
                            cfg=cfg,
                            dom=dom,
                            elapsed_s=float(step_time_s),
                            cum_rain_vol_m3=float(cum_rain),
                            cum_runoff_vol_m3=float(cum_runoff),
                            cum_outflow_vol_m3=float(cum_outflow),
                            P_cum_mm_full=P_slab,
                            Q_cum_mm_full=Q_slab,
                            particles_all=particles,
                        )
                else:
                    P_full = gather_field_partition_to_rank0(comm, partition, P_slab, nrows, ncols)  # set P_full
                    Q_full = gather_field_partition_to_rank0(comm, partition, Q_slab, nrows, ncols)  # set Q_full
                    particles_all = gather_particles_to_rank0(comm, particles)  # set particles_all
                    if rank == 0:  # check condition rank == 0:
                        save_restart_netcdf_rank0(  # execute statement
                            out_path=str(rst_out),  # set out_path
                            cfg=cfg,  # set cfg
                            dom=dom,  # set dom
                            elapsed_s=float(step_time_s),  # set elapsed_s
                        cum_rain_vol_m3=float(cum_rain),  # set cum_rain_vol_m3
                        cum_runoff_vol_m3=float(cum_runoff),  # set cum_runoff_vol_m3
                        cum_outflow_vol_m3=float(cum_outflow),  # set cum_outflow_vol_m3
                        P_cum_mm_full=P_full,  # type: ignore[arg-type]  # set P_cum_mm_full
                        Q_cum_mm_full=Q_full,  # type: ignore[arg-type]  # set Q_cum_mm_full
                        particles_all=particles_all,  # set particles_all
                    )  # execute statement

        # Periodic output writing if requested.
        is_last_step = (k + 1) == steps
        outputs_written_this_step = False
        if output_interval_s > 0.0 and next_output_s is not None:
            while step_time_s + 1e-9 >= next_output_s:
                wrote = _write_output(step_time_s)
                outputs_written_this_step = outputs_written_this_step or wrote
                output_written_any = output_written_any or wrote
                if is_last_step:
                    output_written_on_final_step = output_written_on_final_step or wrote
                next_output_s += output_interval_s

        if is_last_step and out_nc:
            if output_interval_s <= 0.0 or not outputs_written_this_step:
                if _write_output(step_time_s):
                    output_written_any = True
                    output_written_on_final_step = True

    # ------------------------------
    # Final output guard (in case no write happened yet)
    # ------------------------------
    final_elapsed_s = elapsed_s0 + remaining_s  # set final_elapsed_s
    if not out_nc:  # check condition not out_nc:
        if rank == 0:  # check condition rank == 0:
            logger.info("No output.out_netcdf configured; skipping final output.")  # execute statement
    else:
        if steps == 0 and not output_written_any:  # check condition steps == 0 and not output_written_any:
            _write_output(final_elapsed_s)  # execute statement
        if not output_written_on_final_step:  # check condition not output_written_on_final_step:
            _write_output(final_elapsed_s)  # execute statement

    # ------------------------------
    # Final diagnostics and optional GeoJSON export
    # ------------------------------
    final_system_vol_local = float(particles.vol.sum())
    active_particles_local = int(particles.r.size)
    if size > 1:
        final_system_vol = float(comm.allreduce(final_system_vol_local, op=MPI.SUM))
        active_particles = int(comm.allreduce(active_particles_local, op=MPI.SUM))
    else:
        final_system_vol = final_system_vol_local
        active_particles = active_particles_local

    if not mass_error_series:
        expected = cum_runoff - cum_outflow
        mass_error_series.append(float(final_system_vol - expected))

    merged_interval_counts: dict[int, Counter[tuple[int, int]]] = defaultdict(Counter)
    if record_outflow_points:
        if size > 1:
            gathered = comm.gather(outflow_interval_counts, root=0)
            if rank == 0:
                for counts in gathered:
                    for idx, counter in counts.items():
                        merged_interval_counts[idx].update(counter)
        else:
            merged_interval_counts = outflow_interval_counts

    total_wall_local = perf_counter() - run_wall_start  # set total_wall_local
    total_wall = total_wall_local  # set total_wall
    if (metrics_enabled or assistant_metrics_enabled) and size > 1:  # check condition (metrics_enabled or assistant_metrics_enabled) and size > 1:
        total_wall = float(comm.allreduce(total_wall_local, op=MPI.MAX))  # set total_wall

    if metrics_enabled:
        wall_arr = np.asarray(step_wall_samples, dtype=np.float64) if step_wall_samples else np.asarray([0.0])
        mig_arr = np.asarray(step_migration_ratio, dtype=np.float64) if step_migration_ratio else np.asarray([0.0])

        seg_runoff = float(np.sum([sample["segments"]["runoff_spawn_s"] for sample in sampled_steps]))
        seg_advect = float(np.sum([sample["segments"]["advection_s"] for sample in sampled_steps]))
        seg_mig = float(np.sum([sample["segments"]["migration_s"] for sample in sampled_steps]))
        if size > 1:
            seg_runoff = float(comm.allreduce(seg_runoff, op=MPI.MAX))
            seg_advect = float(comm.allreduce(seg_advect, op=MPI.MAX))
            seg_mig = float(comm.allreduce(seg_mig, op=MPI.MAX))
    if rank == 0:
        report = _build_quality_report(
            final_elapsed_s=final_elapsed_s,
            system_vol_final=final_system_vol,
            mass_error_series=mass_error_series,
            cum_rain=cum_rain,
            cum_runoff=cum_runoff,
            cum_outflow=cum_outflow,
            total_hops=total_hops,
            total_spawned_particles=total_spawned_particles,
            total_outflow_particles=total_outflow_particles,
            active_particles=active_particles,
            steps_run=steps,
        )
        _log_quality_report(report)
        if record_outflow_points:
            _write_outflow_geojson(merged_interval_counts, final_elapsed_s)
        if metrics_enabled:
            metrics_report = {
                "scenario": {
                    "domain": domain_label,
                    "shape_rc": [int(nrows), int(ncols)],
                    "ranks": int(size),
                    "mpi": {
                        "enabled": bool(mpi_active),
                        "size": int(size),
                        "world_size": int(world_size),
                        "decomposition": mpi_cfg.decomposition,
                        "min_rows_per_rank": int(mpi_cfg.min_rows_per_rank),
                    "partition": partition.to_dict(),
                },
                "parallelization": {
                    "schema": parallel_schema,
                    "io_mode": io_mode,
                },
                "shared_memory": {
                    "enabled": bool(shared_cfg.enabled),
                    "auto_enabled": bool(shared_cfg.auto_enabled),
                    "workers": int(shared_cfg.workers),
                    "chunk_size": int(shared_cfg.chunk_size),
                        "min_particles_per_worker": int(shared_cfg.min_particles_per_worker),
                    },
                    "device": device,
                    "dt_s": float(dt_s),
                    "T_s": float(T_s),
                    "steps": int(steps),
                },
                "summary": {
                    "wall_clock_s": float(total_wall),
                    "wall_clock_per_step_mean_s": float(np.mean(wall_arr)),
                    "wall_clock_per_step_p50_s": float(np.percentile(wall_arr, 50)),
                    "wall_clock_per_step_p95_s": float(np.percentile(wall_arr, 95)),
                    "particle_throughput_hops_per_s": float(total_hops / max(total_wall, 1e-9)),
                    "migration_ratio_pct_mean": float(np.mean(mig_arr)),
                    "migration_ratio_pct_p95": float(np.percentile(mig_arr, 95)),
                    "gpu_runoff_time_s": float(seg_runoff),
                    "advect_time_s": float(seg_advect),
                    "migration_time_s": float(seg_mig),
                    "spawned_particles": int(total_spawned_particles),
                    "outflow_particles": int(total_outflow_particles),
                },
                "per_step_samples": sampled_steps,
                "notes": [
                    "wall_clock metrics use max across ranks when MPI is active",
                    "throughput uses total hops divided by wall time, independent of simulation time",
                    "migration_ratio_pct reflects particles that changed MPI rank per step",
                    "segments are approximate and exclude synchronization overhead outside each block",
                ],
            }
            logger.info(
                "Parallelization metrics (%s JSON): %s",
                metrics_format,
                json.dumps(metrics_report, indent=metrics_indent, separators=metrics_separators),
            )
            if metrics_output:
                Path(metrics_output).parent.mkdir(parents=True, exist_ok=True)
                with open(metrics_output, "w", encoding="utf-8") as f:
                    json.dump(metrics_report, f, indent=metrics_indent, separators=metrics_separators)
                logger.info("Parallelization metrics written to %s using %s format", metrics_output, metrics_format)
        if assistant_metrics_enabled:
            active_cells = int(np.count_nonzero(dom.active_mask))  # set active_cells
            total_cells = int(dom.active_mask.size)  # set total_cells
            if np.isscalar(dom.cell_area_m2):  # check condition np.isscalar(dom.cell_area_m2):
                cell_area_summary = {  # set cell_area_summary
                    "mode": "scalar",  # execute statement
                    "value_m2": float(dom.cell_area_m2),  # execute statement
                }  # execute statement
            else:  # fallback branch
                area = np.asarray(dom.cell_area_m2, dtype=np.float64)  # set area
                active_area = area[dom.active_mask] if area.shape == dom.active_mask.shape else area  # set active_area
                active_area = active_area[np.isfinite(active_area)]  # set active_area
                cell_area_summary = {  # set cell_area_summary
                    "mode": "grid",  # execute statement
                    "min_m2": float(np.nanmin(active_area)) if active_area.size else float("nan"),  # execute statement
                    "mean_m2": float(np.nanmean(active_area)) if active_area.size else float("nan"),  # execute statement
                    "max_m2": float(np.nanmax(active_area)) if active_area.size else float("nan"),  # execute statement
                }  # execute statement
            travel_time_summary: dict[str, Any]  # set travel_time_summary type
            if np.isscalar(travel_time_s):  # check condition np.isscalar(travel_time_s):
                travel_time_summary = {  # set travel_time_summary
                    "mode": travel_time_mode,  # execute statement
                    "hillslope_s": float(travel_time_s),  # execute statement
                    "channel_s": float(travel_time_channel_s),  # execute statement
                }  # execute statement
            else:  # fallback branch
                hillslope_vals = travel_time_s[valid]  # set hillslope_vals
                channel_vals = travel_time_channel_s[valid]  # set channel_vals
                hillslope_safe = hillslope_vals if hillslope_vals.size else np.asarray([np.nan], dtype=np.float32)  # set hillslope_safe
                channel_safe = channel_vals if channel_vals.size else np.asarray([np.nan], dtype=np.float32)  # set channel_safe
                travel_time_summary = {  # set travel_time_summary
                    "mode": travel_time_mode,  # execute statement
                    "hillslope_median_s": float(np.nanmedian(hillslope_safe)),  # execute statement
                    "channel_median_s": float(np.nanmedian(channel_safe)),  # execute statement
                    "hillslope_bounds_s": [float(np.nanmin(hillslope_safe)), float(np.nanmax(hillslope_safe))],  # execute statement
                    "channel_bounds_s": [float(np.nanmin(channel_safe)), float(np.nanmax(channel_safe))],  # execute statement
                }  # execute statement
            hyd_bal = report["volume_balance"]  # set hyd_bal
            hyd_status = report["hydrology"]  # set hyd_status
            assistant_report = {  # set assistant_report
                "domain": {  # execute statement
                    "name": domain_label,  # execute statement
                    "shape_rc": [int(nrows), int(ncols)],  # execute statement
                    "active_cells": active_cells,  # execute statement
                    "total_cells": total_cells,  # execute statement
                    "cell_area_m2": cell_area_summary,  # execute statement
                },  # execute statement
                "time": {  # execute statement
                    "start_time": str(mcfg.get("start_time", "")),  # execute statement
                    "dt_s": float(dt_s),  # execute statement
                    "simulation_window_s": float(T_s),  # execute statement
                    "steps_executed": int(steps),  # execute statement
                    "elapsed_wall_clock_s": float(total_wall),  # execute statement
                },  # execute statement
                "compute": {  # execute statement
                    "device": device,  # execute statement
                    "mpi_enabled": bool(mpi_active),  # execute statement
                    "mpi_ranks": int(size),  # execute statement
                    "mpi_world_size": int(world_size),  # execute statement
                    "shared_memory_workers": int(shared_cfg.workers),  # execute statement
                    "travel_time": travel_time_summary,  # execute statement
                    "particle_volume_m3": float(particle_vol_m3),  # execute statement
                },  # execute statement
                "hydrology": {  # execute statement
                    "initial_abstraction_ratio": float(ia_ratio),  # execute statement
                    "outflow_sink_enabled": bool(outflow_sink),  # execute statement
                    "rain_m3": hyd_bal["rain_m3"],  # execute statement
                    "runoff_m3": hyd_bal["runoff_m3"],  # execute statement
                    "outflow_m3": hyd_bal["outflow_m3"],  # execute statement
                    "stored_m3": hyd_bal["stored_m3"],  # execute statement
                    "runoff_to_rain_ratio": hyd_status["runoff_to_rain_ratio"],  # execute statement
                    "mass_error_m3": hyd_bal["final_mass_error_m3"],  # execute statement
                    "mass_error_pct_of_runoff": hyd_bal["mass_error_pct_of_runoff"],  # execute statement
                    "hydrology_ok": bool(hyd_status["hydrology_ok"]),  # execute statement
                },  # execute statement
                "status": report["status"],  # execute statement
                "notes": [  # execute statement
                    "GPT-friendly summary for AI-assisted review of computational and hydrological behavior.",  # execute statement
                    "Travel-time statistics are derived from active cells when automatic mode is used.",  # execute statement
                ],  # execute statement
            }  # execute statement
            logger.info(
                "AI assistant metrics (%s JSON): %s",
                assistant_metrics_format,
                json.dumps(assistant_report, indent=assistant_indent, separators=assistant_separators),
            )
            if assistant_metrics_output:
                Path(assistant_metrics_output).parent.mkdir(parents=True, exist_ok=True)
                with open(assistant_metrics_output, "w", encoding="utf-8") as f:
                    json.dump(assistant_report, f, indent=assistant_indent, separators=assistant_separators)
                logger.info("AI assistant metrics written to %s using %s format", assistant_metrics_output, assistant_metrics_format)
