# -*- coding: utf-8 -*-
"""Core simulation driver for LPERFECT."""

# Import typing primitives.
from typing import Any, Dict

# Import logging.
import logging

# Import numpy.
import numpy as np

# Import local modules.
from .time_utils import parse_iso8601_to_datetime64
from .d8 import build_downstream_index
from .rain import build_rain_sources, blended_rain_step_mm_rank0
from .runoff import scs_cn_cumulative_runoff_mm
from .compute_backend import gpu_available, normalize_device
from .hydraulics import spawn_particles_from_runoff_slab, advect_particles_one_step, local_volgrid_from_particles_slab
from .particles import Particles, empty_particles, concat_particles
from .risk import compute_flow_accum_area_m2, compute_risk_index
from .io_netcdf import write_results_netcdf_rank0, save_restart_netcdf_rank0, load_restart_netcdf_rank0

# Import MPI helpers.
from .mpi_utils import (
    HAVE_MPI,
    MPI,
    slab_bounds,
    scatter_field_slab,
    gather_field_slab_to_rank0,
    gather_particles_to_rank0,
    scatter_particles_from_rank0,
    migrate_particles_slab,
)

# Import Domain class.
from .domain import Domain

# Create a logger for this module.
logger = logging.getLogger("lperfect")


def run_simulation(comm: Any, rank: int, size: int, cfg: Dict[str, Any], dom: Domain) -> None:
    """Run the full LPERFECT simulation (serial or MPI)."""

    # ------------------------------
    # Model parameters
    # ------------------------------
    mcfg = cfg["model"]
    T_s = float(mcfg["T_s"])                                  # Total simulation time.
    dt_s = float(mcfg["dt_s"])                                # Timestep.
    encoding = str(mcfg["encoding"])                           # D8 encoding.
    ia_ratio = float(mcfg["ia_ratio"])                         # Initial abstraction ratio.
    particle_vol_m3 = float(mcfg["particle_vol_m3"])           # Target particle volume.
    travel_time_s = float(mcfg["travel_time_s"])               # Hillslope hop time.
    travel_time_channel_s = float(mcfg["travel_time_channel_s"])# Channel hop time.
    outflow_sink = bool(mcfg["outflow_sink"])                  # Drop particles leaving domain.
    log_every = int(mcfg.get("log_every", 10))                 # Diagnostics frequency.

    # Resolve compute device.
    compute_cfg = cfg.get("compute", {})
    device = normalize_device(compute_cfg.get("device", "cpu"))
    if device == "gpu" and not gpu_available():
        if rank == 0:
            logger.warning("GPU requested but CuPy not available; falling back to CPU.")
        device = "cpu"
    if rank == 0:
        logger.info("Compute device: %s", device)

    # Start time for time-aware rain selection.
    start_time = parse_iso8601_to_datetime64(mcfg.get("start_time", None))

    # Domain shape.
    nrows, ncols = dom.dem.shape

    # Precompute downstream lookup (replicated).
    valid, ds_r, ds_c = build_downstream_index(dom.d8, encoding)

    # Rain sources configuration (replicated config, rank0 reads files).
    rain_sources = build_rain_sources(cfg)

    # ------------------------------
    # Initialize state (possibly from restart)
    # ------------------------------
    rst_cfg = cfg.get("restart", {})
    restart_in = rst_cfg.get("in", None)

    if size == 1:
        # Serial state uses full arrays.
        if restart_in:
            r = load_restart_netcdf_rank0(restart_in)
            P_slab = r["P_cum_mm"]
            Q_slab = r["Q_cum_mm"]
            particles = Particles(r=r["r"], c=r["c"], vol=r["vol"], tau=r["tau"])
            cum_rain = r["cum_rain_vol_m3"]
            cum_runoff = r["cum_runoff_vol_m3"]
            cum_outflow = r["cum_outflow_vol_m3"]
            elapsed_s0 = r["elapsed_s"]
        else:
            P_slab = np.zeros((nrows, ncols), dtype=np.float64)
            Q_slab = np.zeros((nrows, ncols), dtype=np.float64)
            particles = empty_particles()
            cum_rain = cum_runoff = cum_outflow = 0.0
            elapsed_s0 = 0.0
        r0, r1 = 0, nrows
    else:
        # MPI: rank0 loads restart and scatters state.
        if restart_in and rank == 0:
            r = load_restart_netcdf_rank0(restart_in)

            # Optional strict checks.
            if bool(rst_cfg.get("strict_grid_check", True)):
                if r["P_cum_mm"].shape != (nrows, ncols):
                    raise ValueError("Restart grid mismatch for P_cum_mm")
                if r["Q_cum_mm"].shape != (nrows, ncols):
                    raise ValueError("Restart grid mismatch for Q_cum_mm")                

            P_full = r["P_cum_mm"]
            Q_full = r["Q_cum_mm"]
            particles_all = Particles(r=r["r"], c=r["c"], vol=r["vol"], tau=r["tau"])
            cum_rain = r["cum_rain_vol_m3"]
            cum_runoff = r["cum_runoff_vol_m3"]
            cum_outflow = r["cum_outflow_vol_m3"]
            elapsed_s0 = r["elapsed_s"]
        else:
            P_full = None
            Q_full = None
            particles_all = None
            cum_rain = cum_runoff = cum_outflow = 0.0
            elapsed_s0 = 0.0

        # Broadcast scalar diagnostics and elapsed time.
        scal = np.array([cum_rain, cum_runoff, cum_outflow, elapsed_s0], dtype=np.float64) if rank == 0 else np.zeros(4, dtype=np.float64)
        comm.Bcast(scal, root=0)
        cum_rain, cum_runoff, cum_outflow, elapsed_s0 = map(float, scal.tolist())

        # Scatter slab fields.
        P_slab = scatter_field_slab(comm, P_full, nrows, ncols, np.float64)
        Q_slab = scatter_field_slab(comm, Q_full, nrows, ncols, np.float64)

        # Scatter particles by row ownership.
        particles = scatter_particles_from_rank0(comm, particles_all, nrows=nrows)

        # Local slab bounds.
        r0, r1 = slab_bounds(nrows, size, rank)

    # ------------------------------
    # Main time loop
    # ------------------------------
    if elapsed_s0 > T_s + 1e-9:
        raise ValueError(f"Restart elapsed_s={elapsed_s0} exceeds total T_s={T_s}.")

    remaining_s = max(0.0, T_s - elapsed_s0)
    steps = int(np.ceil(remaining_s / dt_s))

    if rank == 0:
        logger.info("LPERFECT start: T_s=%.1f dt_s=%.1f steps=%d ranks=%d", T_s, dt_s, steps, size)

    for k in range(steps):
        # Elapsed time at end of step.
        step_time_s = elapsed_s0 + min((k + 1) * dt_s, remaining_s)

        # Absolute timestamp (for rain selection) if configured.
        sim_time = (start_time + np.timedelta64(int(round(step_time_s)), "s")) if start_time is not None else None

        # Rank0 reads and blends rainfall, then broadcasts.
        if rank == 0:
            rain_step_mm = blended_rain_step_mm_rank0(
                sources=rain_sources,
                shape=(nrows, ncols),
                dt_s=dt_s,
                step_idx=k,
                sim_time=sim_time,
            )
            # Force inactive cells to zero rainfall.
            rain_step_mm = np.where(dom.active_mask, rain_step_mm, 0.0).astype(np.float32)
        else:
            rain_step_mm = np.empty((nrows, ncols), dtype=np.float32)

        if size > 1:
            comm.Bcast(rain_step_mm, root=0)

        # Slice rainfall to local slab.
        rain_slab_mm = rain_step_mm[r0:r1, :].astype(np.float64)

        # Compute rain volume injected this step (for diagnostics).
        if np.isscalar(dom.cell_area_m2):
            rain_vol_step_local = float((rain_slab_mm / 1000.0).sum() * float(dom.cell_area_m2))
        else:
            rain_vol_step_local = float(((rain_slab_mm / 1000.0) * dom.cell_area_m2[r0:r1, :]).sum())

        # Update cumulative precipitation.
        P_slab = P_slab + rain_slab_mm

        # Compute cumulative runoff with CN and incremental runoff for this step.
        CN_slab = dom.cn[r0:r1, :]
        Q_cum_slab = scs_cn_cumulative_runoff_mm(P_slab, CN_slab, ia_ratio=ia_ratio, device=device)
        dQ_mm = np.maximum(Q_cum_slab - Q_slab, 0.0)
        Q_slab = Q_cum_slab

        # Convert incremental runoff to meters.
        runoff_depth_m = dQ_mm / 1000.0

        # Spawn new particles from incremental runoff.
        newp, spawned_vol_local = spawn_particles_from_runoff_slab(
            runoff_depth_m_slab=runoff_depth_m,
            r0=r0,
            cell_area_m2=dom.cell_area_m2,
            particle_vol_m3=particle_vol_m3,
            active_mask_global=dom.active_mask,
        )
        particles = concat_particles(particles, newp)

        # Advect particles.
        particles, outflow_vol_local, nhops_local = advect_particles_one_step(
            particles=particles,
            valid=valid,
            ds_r=ds_r,
            ds_c=ds_c,
            dt_s=dt_s,
            travel_time_s=travel_time_s,
            travel_time_channel_s=travel_time_channel_s,
            channel_mask=dom.channel_mask,
            outflow_sink=outflow_sink,
        )

        # Migrate particles between slabs (MPI).
        if size > 1:
            particles = migrate_particles_slab(comm, particles, nrows=nrows)

        # Compute current system volume (sum of particle volumes).
        system_vol_local = float(particles.vol.sum())

        # Reduce diagnostics to global totals.
        if size > 1:
            rain_vol_step = comm.allreduce(rain_vol_step_local, op=MPI.SUM)
            runoff_vol_step = comm.allreduce(spawned_vol_local, op=MPI.SUM)
            outflow_vol_step = comm.allreduce(outflow_vol_local, op=MPI.SUM)
            system_vol = comm.allreduce(system_vol_local, op=MPI.SUM)
            nhops = comm.allreduce(nhops_local, op=MPI.SUM)
        else:
            rain_vol_step = rain_vol_step_local
            runoff_vol_step = spawned_vol_local
            outflow_vol_step = outflow_vol_local
            system_vol = system_vol_local
            nhops = nhops_local

        # Update cumulative totals.
        cum_rain += rain_vol_step
        cum_runoff += runoff_vol_step
        cum_outflow += outflow_vol_step

        # Mass balance error: stored volume minus (generated - outflow).
        expected = cum_runoff - cum_outflow
        mass_error = system_vol - expected

        # Periodic diagnostics (rank0 only).
        if rank == 0 and log_every > 0 and ((k + 1) % log_every == 0 or (k + 1) == steps):
            logger.info(
                "mass_balance step=%d time_s=%.1f rain_m3=%.3e runoff_m3=%.3e outflow_m3=%.3e system_m3=%.3e err_m3=%.3e hops=%d",
                (k + 1), step_time_s, rain_vol_step, runoff_vol_step, outflow_vol_step, system_vol, mass_error, nhops
            )

        # Restart writing if enabled.
        rst_out = rst_cfg.get("out", None)
        rst_every = int(rst_cfg.get("every", 0))
        if rst_out and rst_every > 0 and (((k + 1) % rst_every == 0) or ((k + 1) == steps)):
            if size == 1:
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
                    particles_all=particles,
                )
            else:
                P_full = gather_field_slab_to_rank0(comm, P_slab, nrows, ncols)
                Q_full = gather_field_slab_to_rank0(comm, Q_slab, nrows, ncols)
                particles_all = gather_particles_to_rank0(comm, particles)
                if rank == 0:
                    save_restart_netcdf_rank0(
                        out_path=str(rst_out),
                        cfg=cfg,
                        dom=dom,
                        elapsed_s=float(step_time_s),
                        cum_rain_vol_m3=float(cum_rain),
                        cum_runoff_vol_m3=float(cum_runoff),
                        cum_outflow_vol_m3=float(cum_outflow),
                        P_cum_mm_full=P_full,  # type: ignore[arg-type]
                        Q_cum_mm_full=Q_full,  # type: ignore[arg-type]
                        particles_all=particles_all,
                    )

    # ------------------------------
    # Final output gathering and writing (rank0)
    # ------------------------------
    out_nc = cfg.get("output", {}).get("out_netcdf", None)
    if not out_nc:
        if rank == 0:
            logger.info("No output.out_netcdf configured; skipping final output.")
        return

    risk_cfg = cfg.get("risk", {})
    do_risk = bool(risk_cfg.get("enabled", True))

    if size == 1:
        volgrid = local_volgrid_from_particles_slab(particles, r0=0, r1=nrows, ncols=ncols)
        if np.isscalar(dom.cell_area_m2):
            flood = volgrid / float(dom.cell_area_m2)
        else:
            flood = np.divide(volgrid, dom.cell_area_m2, out=np.zeros_like(volgrid), where=(dom.cell_area_m2 > 0))
        flood = np.where(dom.active_mask, flood, np.nan)

        if do_risk:
            acc = compute_flow_accum_area_m2(dom.d8, encoding, dom.cell_area_m2, dom.active_mask)
            risk = compute_risk_index(
                runoff_cum_mm=Q_slab,
                flow_accum_m2=acc,
                active_mask=dom.active_mask,
                balance=float(risk_cfg.get("balance", 0.55)),
                p_low=float(risk_cfg.get("p_low", 5.0)),
                p_high=float(risk_cfg.get("p_high", 95.0)),
            )
        else:
            risk = np.full_like(flood, np.nan, dtype=np.float64)

        write_results_netcdf_rank0(out_nc, cfg, dom, flood, risk)

    else:
        Q_full = gather_field_slab_to_rank0(comm, Q_slab, nrows, ncols)
        volgrid_slab = local_volgrid_from_particles_slab(particles, r0=r0, r1=r1, ncols=ncols)
        vol_full = gather_field_slab_to_rank0(comm, volgrid_slab, nrows, ncols)

        if rank == 0:
            if np.isscalar(dom.cell_area_m2):
                flood = vol_full / float(dom.cell_area_m2)  # type: ignore[operator]
            else:
                flood = np.divide(vol_full, dom.cell_area_m2, out=np.zeros_like(vol_full), where=(dom.cell_area_m2 > 0))  # type: ignore[arg-type]
            flood = np.where(dom.active_mask, flood, np.nan)

            if do_risk:
                acc = compute_flow_accum_area_m2(dom.d8, encoding, dom.cell_area_m2, dom.active_mask)
                risk = compute_risk_index(
                    runoff_cum_mm=Q_full,  # type: ignore[arg-type]
                    flow_accum_m2=acc,
                    active_mask=dom.active_mask,
                    balance=float(risk_cfg.get("balance", 0.55)),
                    p_low=float(risk_cfg.get("p_low", 5.0)),
                    p_high=float(risk_cfg.get("p_high", 95.0)),
                )
            else:
                risk = np.full_like(flood, np.nan, dtype=np.float64)

            write_results_netcdf_rank0(out_nc, cfg, dom, flood, risk)
