# -*- coding: utf-8 -*-
"""Core simulation driver for LPERFECT."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import typing primitives.
from typing import Any, Dict  # import typing import Any, Dict

# Import logging.
import logging  # import logging

# Import numpy.
import numpy as np  # import numpy as np

# Import local modules.
from .time_utils import parse_iso8601_to_datetime64  # import .time_utils import parse_iso8601_to_datetime64
from .d8 import build_downstream_index  # import .d8 import build_downstream_index
from .rain import build_rain_sources, blended_rain_step_mm_rank0  # import .rain import build_rain_sources, blended_rain_step_mm_rank0
from .runoff import scs_cn_cumulative_runoff_mm  # import .runoff import scs_cn_cumulative_runoff_mm
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
    slab_bounds,  # execute statement
    scatter_field_slab,  # execute statement
    gather_field_slab_to_rank0,  # execute statement
    gather_particles_to_rank0,  # execute statement
    scatter_particles_from_rank0,  # execute statement
    migrate_particles_slab,  # execute statement
)  # execute statement

# Import Domain class.
from .domain import Domain  # import .domain import Domain

# Create a logger for this module.
logger = logging.getLogger("lperfect")  # set logger


def run_simulation(comm: Any, rank: int, size: int, cfg: Dict[str, Any], dom: Domain) -> None:  # define function run_simulation
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
    if device == "gpu" and not gpu_available():  # check condition device == "gpu" and not gpu_available():
        if rank == 0:  # check condition rank == 0:
            logger.warning("GPU requested but CuPy not available; falling back to CPU.")  # execute statement
        device = "cpu"  # set device
    if rank == 0:  # check condition rank == 0:
        logger.info("Compute device: %s", device)  # execute statement
    shared_cfg = SharedMemoryConfig.from_dict(compute_cfg.get("shared_memory", {}))  # set shared_cfg
    if rank == 0 and shared_cfg.enabled and shared_cfg.workers > 1:  # check condition rank == 0 and shared_cfg.enabled and shared_cfg.workers > 1:
        logger.info(
            "Shared-memory parallelism: enabled with %d workers (chunk=%d, min_particles=%d)",
            shared_cfg.workers,
            shared_cfg.chunk_size,
            shared_cfg.min_particles_per_worker,
        )  # execute statement

    # Start time for time-aware rain selection.
    start_time = parse_iso8601_to_datetime64(mcfg.get("start_time", None))  # set start_time

    # Domain shape.
    nrows, ncols = dom.dem.shape  # set nrows, ncols

    # Precompute downstream lookup (replicated).
    valid, ds_r, ds_c = build_downstream_index(dom.d8, encoding)  # set valid, ds_r, ds_c

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
            dist = np.zeros_like(valid, dtype=np.float64)
            dist[orth] = base
            dist[diag] = base * np.sqrt(2.0)
            return dist

        base = np.sqrt(cell_area_m2.astype(np.float64))
        dist = np.zeros_like(base, dtype=np.float64)
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
        hill = np.clip(dist / vel_hill, min_s, max_s)
        ch = np.clip(dist / vel_ch, min_s, max_s)
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

    # ------------------------------
    # Initialize state (possibly from restart)
    # ------------------------------
    rst_cfg = cfg.get("restart", {})  # set rst_cfg
    restart_in = rst_cfg.get("in", None)  # set restart_in

    if size == 1:  # check condition size == 1:
        # Serial state uses full arrays.
        if restart_in:  # check condition restart_in:
            r = load_restart_netcdf_rank0(restart_in)  # set r
            P_slab = r["P_cum_mm"]  # set P_slab
            Q_slab = r["Q_cum_mm"]  # set Q_slab
            particles = Particles(r=r["r"], c=r["c"], vol=r["vol"], tau=r["tau"])  # set particles
            cum_rain = r["cum_rain_vol_m3"]  # set cum_rain
            cum_runoff = r["cum_runoff_vol_m3"]  # set cum_runoff
            cum_outflow = r["cum_outflow_vol_m3"]  # set cum_outflow
            elapsed_s0 = r["elapsed_s"]  # set elapsed_s0
        else:  # fallback branch
            P_slab = np.zeros((nrows, ncols), dtype=np.float64)  # set P_slab
            Q_slab = np.zeros((nrows, ncols), dtype=np.float64)  # set Q_slab
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

            P_full = r["P_cum_mm"]  # set P_full
            Q_full = r["Q_cum_mm"]  # set Q_full
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

        # Scatter slab fields.
        P_slab = scatter_field_slab(comm, P_full, nrows, ncols, np.float64)  # set P_slab
        Q_slab = scatter_field_slab(comm, Q_full, nrows, ncols, np.float64)  # set Q_slab

        # Scatter particles by row ownership.
        particles = scatter_particles_from_rank0(comm, particles_all, nrows=nrows)  # set particles

        # Local slab bounds.
        r0, r1 = slab_bounds(nrows, size, rank)  # set r0, r1

    # ------------------------------
    # Main time loop
    # ------------------------------
    if elapsed_s0 > T_s + 1e-9:  # check condition elapsed_s0 > T_s + 1e-9:
        raise ValueError(f"Restart elapsed_s={elapsed_s0} exceeds total T_s={T_s}.")  # raise ValueError(f"Restart elapsed_s={elapsed_s0} exceeds total T_s={T_s}.")

    remaining_s = max(0.0, T_s - elapsed_s0)  # set remaining_s
    steps = int(np.ceil(remaining_s / dt_s))  # set steps

    if rank == 0:  # check condition rank == 0:
        logger.info("LPERFECT start: T_s=%.1f dt_s=%.1f steps=%d ranks=%d", T_s, dt_s, steps, size)  # execute statement

    for k in range(steps):  # loop over k in range(steps):
        # Elapsed time at end of step.
        step_time_s = elapsed_s0 + min((k + 1) * dt_s, remaining_s)  # set step_time_s

        # Absolute timestamp (for rain selection) if configured.
        sim_time = (start_time + np.timedelta64(int(round(step_time_s)), "s")) if start_time is not None else None  # set sim_time

        # Rank0 reads and blends rainfall, then broadcasts.
        if rank == 0:  # check condition rank == 0:
            rain_step_mm = blended_rain_step_mm_rank0(  # set rain_step_mm
                sources=rain_sources,  # set sources
                shape=(nrows, ncols),  # set shape
                dt_s=dt_s,  # set dt_s
                step_idx=k,  # set step_idx
                sim_time=sim_time,  # set sim_time
            )  # execute statement
            # Force inactive cells to zero rainfall.
            rain_step_mm = np.where(dom.active_mask, rain_step_mm, 0.0).astype(np.float32)  # set rain_step_mm
        else:  # fallback branch
            rain_step_mm = np.empty((nrows, ncols), dtype=np.float32)  # set rain_step_mm

        if size > 1:  # check condition size > 1:
            comm.Bcast(rain_step_mm, root=0)  # execute statement

        # Slice rainfall to local slab.
        rain_slab_mm = rain_step_mm[r0:r1, :].astype(np.float64)  # set rain_slab_mm

        # Compute rain volume injected this step (for diagnostics).
        if np.isscalar(dom.cell_area_m2):  # check condition np.isscalar(dom.cell_area_m2):
            rain_vol_step_local = float((rain_slab_mm / 1000.0).sum() * float(dom.cell_area_m2))  # set rain_vol_step_local
        else:  # fallback branch
            rain_vol_step_local = float(((rain_slab_mm / 1000.0) * dom.cell_area_m2[r0:r1, :]).sum())  # set rain_vol_step_local

        # Update cumulative precipitation.
        P_slab = P_slab + rain_slab_mm  # set P_slab

        # Compute cumulative runoff with CN and incremental runoff for this step.
        CN_slab = dom.cn[r0:r1, :]  # set CN_slab
        Q_cum_slab = scs_cn_cumulative_runoff_mm(P_slab, CN_slab, ia_ratio=ia_ratio, device=device)  # set Q_cum_slab
        dQ_mm = np.maximum(Q_cum_slab - Q_slab, 0.0)  # set dQ_mm
        Q_slab = Q_cum_slab  # set Q_slab

        # Convert incremental runoff to meters.
        runoff_depth_m = dQ_mm / 1000.0  # set runoff_depth_m

        # Spawn new particles from incremental runoff.
        newp, spawned_vol_local = spawn_particles_from_runoff_slab(  # set newp, spawned_vol_local
            runoff_depth_m_slab=runoff_depth_m,  # set runoff_depth_m_slab
            r0=r0,  # set r0
            cell_area_m2=dom.cell_area_m2,  # set cell_area_m2
            particle_vol_m3=particle_vol_m3,  # set particle_vol_m3
            active_mask_global=dom.active_mask,  # set active_mask_global
        )  # execute statement
        particles = concat_particles(particles, newp)  # set particles

        # Advect particles.
        particles, outflow_vol_local, nhops_local = advect_particles_one_step(  # set particles, outflow_vol_local, nhops_local
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
        )  # execute statement

        # Migrate particles between slabs (MPI).
        if size > 1:  # check condition size > 1:
            particles = migrate_particles_slab(comm, particles, nrows=nrows)  # set particles

        # Compute current system volume (sum of particle volumes).
        system_vol_local = float(particles.vol.sum())  # set system_vol_local

        # Reduce diagnostics to global totals.
        if size > 1:  # check condition size > 1:
            rain_vol_step = comm.allreduce(rain_vol_step_local, op=MPI.SUM)  # set rain_vol_step
            runoff_vol_step = comm.allreduce(spawned_vol_local, op=MPI.SUM)  # set runoff_vol_step
            outflow_vol_step = comm.allreduce(outflow_vol_local, op=MPI.SUM)  # set outflow_vol_step
            system_vol = comm.allreduce(system_vol_local, op=MPI.SUM)  # set system_vol
            nhops = comm.allreduce(nhops_local, op=MPI.SUM)  # set nhops
        else:  # fallback branch
            rain_vol_step = rain_vol_step_local  # set rain_vol_step
            runoff_vol_step = spawned_vol_local  # set runoff_vol_step
            outflow_vol_step = outflow_vol_local  # set outflow_vol_step
            system_vol = system_vol_local  # set system_vol
            nhops = nhops_local  # set nhops

        # Update cumulative totals.
        cum_rain += rain_vol_step  # execute statement
        cum_runoff += runoff_vol_step  # execute statement
        cum_outflow += outflow_vol_step  # execute statement

        # Mass balance error: stored volume minus (generated - outflow).
        expected = cum_runoff - cum_outflow  # set expected
        mass_error = system_vol - expected  # set mass_error

        # Periodic diagnostics (rank0 only).
        if rank == 0 and log_every > 0 and ((k + 1) % log_every == 0 or (k + 1) == steps):  # check condition rank == 0 and log_every > 0 and ((k + 1) % log_every == 0 or (k + 1) == steps):
            logger.info(  # execute statement
                "mass_balance step=%d time_s=%.1f rain_m3=%.3e runoff_m3=%.3e outflow_m3=%.3e system_m3=%.3e err_m3=%.3e hops=%d",  # execute statement
                (k + 1), step_time_s, rain_vol_step, runoff_vol_step, outflow_vol_step, system_vol, mass_error, nhops  # execute statement
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
                P_full = gather_field_slab_to_rank0(comm, P_slab, nrows, ncols)  # set P_full
                Q_full = gather_field_slab_to_rank0(comm, Q_slab, nrows, ncols)  # set Q_full
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

    # ------------------------------
    # Final output gathering and writing (rank0)
    # ------------------------------
    out_nc = cfg.get("output", {}).get("out_netcdf", None)  # set out_nc
    if not out_nc:  # check condition not out_nc:
        if rank == 0:  # check condition rank == 0:
            logger.info("No output.out_netcdf configured; skipping final output.")  # execute statement
        return  # return value

    risk_cfg = cfg.get("risk", {})  # set risk_cfg
    do_risk = bool(risk_cfg.get("enabled", True))  # set do_risk

    if size == 1:  # check condition size == 1:
        volgrid = local_volgrid_from_particles_slab(particles, r0=0, r1=nrows, ncols=ncols, shared_cfg=shared_cfg)  # set volgrid
        if np.isscalar(dom.cell_area_m2):  # check condition np.isscalar(dom.cell_area_m2):
            flood = volgrid / float(dom.cell_area_m2)  # set flood
        else:  # fallback branch
            flood = np.divide(volgrid, dom.cell_area_m2, out=np.zeros_like(volgrid), where=(dom.cell_area_m2 > 0))  # set flood
        flood = np.where(dom.active_mask, flood, np.nan)  # set flood

        if do_risk:  # check condition do_risk:
            acc = compute_flow_accum_area_m2(dom.d8, encoding, dom.cell_area_m2, dom.active_mask)  # set acc
            risk = compute_risk_index(  # set risk
                runoff_cum_mm=Q_slab,  # set runoff_cum_mm
                flow_accum_m2=acc,  # set flow_accum_m2
                active_mask=dom.active_mask,  # set active_mask
                balance=float(risk_cfg.get("balance", 0.55)),  # set balance
                p_low=float(risk_cfg.get("p_low", 5.0)),  # set p_low
                p_high=float(risk_cfg.get("p_high", 95.0)),  # set p_high
            )  # execute statement
        else:  # fallback branch
            risk = np.full_like(flood, np.nan, dtype=np.float64)  # set risk

        write_results_netcdf_rank0(out_nc, cfg, dom, flood, risk)  # execute statement

    else:  # fallback branch
        Q_full = gather_field_slab_to_rank0(comm, Q_slab, nrows, ncols)  # set Q_full
        volgrid_slab = local_volgrid_from_particles_slab(particles, r0=r0, r1=r1, ncols=ncols, shared_cfg=shared_cfg)  # set volgrid_slab
        vol_full = gather_field_slab_to_rank0(comm, volgrid_slab, nrows, ncols)  # set vol_full

        if rank == 0:  # check condition rank == 0:
            if np.isscalar(dom.cell_area_m2):  # check condition np.isscalar(dom.cell_area_m2):
                flood = vol_full / float(dom.cell_area_m2)  # type: ignore[operator]  # set flood
            else:  # fallback branch
                flood = np.divide(vol_full, dom.cell_area_m2, out=np.zeros_like(vol_full), where=(dom.cell_area_m2 > 0))  # type: ignore[arg-type]  # set flood
            flood = np.where(dom.active_mask, flood, np.nan)  # set flood

            if do_risk:  # check condition do_risk:
                acc = compute_flow_accum_area_m2(dom.d8, encoding, dom.cell_area_m2, dom.active_mask)  # set acc
                risk = compute_risk_index(  # set risk
                    runoff_cum_mm=Q_full,  # type: ignore[arg-type]  # set runoff_cum_mm
                    flow_accum_m2=acc,  # set flow_accum_m2
                    active_mask=dom.active_mask,  # set active_mask
                    balance=float(risk_cfg.get("balance", 0.55)),  # set balance
                    p_low=float(risk_cfg.get("p_low", 5.0)),  # set p_low
                    p_high=float(risk_cfg.get("p_high", 95.0)),  # set p_high
                )  # execute statement
            else:  # fallback branch
                risk = np.full_like(flood, np.nan, dtype=np.float64)  # set risk

            write_results_netcdf_rank0(out_nc, cfg, dom, flood, risk)  # execute statement
