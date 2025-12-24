# -*- coding: utf-8 -*-
"""Lagrangian routing helpers: particle spawning and advection."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import typing primitives.
from typing import Optional, Tuple  # import typing import Optional, Tuple

# Import numpy.
import numpy as np  # import numpy as np

# Import local particle structures.
from .particles import Particles, empty_particles  # import .particles import Particles, empty_particles
from .shared_memory import SharedMemoryConfig, chunk_bounds, should_parallelize  # import .shared_memory import SharedMemoryConfig, chunk_bounds, should_parallelize


def _travel_time_increment(  # define function _travel_time_increment
    travel_time_s: float | np.ndarray,  # execute statement
    travel_time_channel_s: float | np.ndarray,  # execute statement
    channel_mask: Optional[np.ndarray],  # execute statement
    r_idx: np.ndarray,  # execute statement
    c_idx: np.ndarray,  # execute statement
) -> np.ndarray:  # execute statement
    """Return per-particle travel time increment for destination cells."""  # execute statement
    if channel_mask is None:  # check condition channel_mask is None:
        if np.isscalar(travel_time_s):  # check condition np.isscalar(travel_time_s):
            return np.full_like(r_idx, float(travel_time_s), dtype=np.float64)  # return np.full_like(r_idx, float(travel_time_s), dtype=np.float64)
        return travel_time_s[r_idx, c_idx].astype(np.float64)  # return travel_time_s[r_idx, c_idx].astype(np.float64)

    is_ch = channel_mask[r_idx, c_idx]  # set is_ch

    def _lookup(tt: float | np.ndarray) -> np.ndarray:  # define function _lookup
        if np.isscalar(tt):  # check condition np.isscalar(tt):
            return np.full_like(r_idx, float(tt), dtype=np.float64)  # return np.full_like(r_idx, float(tt), dtype=np.float64)
        return tt[r_idx, c_idx].astype(np.float64)  # return tt[r_idx, c_idx].astype(np.float64)

    return np.where(is_ch, _lookup(travel_time_channel_s), _lookup(travel_time_s))  # return np.where(is_ch, _lookup(travel_time_channel_s), _lookup(travel_time_s))


def cell_area_at(area: float | np.ndarray, rr: np.ndarray, cc: np.ndarray) -> np.ndarray:  # define function cell_area_at
    """Return cell area(s) for indices rr,cc."""  # execute statement
    if np.isscalar(area):  # check condition np.isscalar(area):
        return np.full(rr.shape, float(area), dtype=np.float64)  # return np.full(rr.shape, float(area), dtype=np.float64)
    return area[rr, cc].astype(np.float64)  # return area[rr, cc].astype(np.float64)


def spawn_particles_from_runoff_slab(  # define function spawn_particles_from_runoff_slab
    runoff_depth_m_slab: np.ndarray,  # execute statement
    r0: int,  # execute statement
    cell_area_m2: float | np.ndarray,  # execute statement
    particle_vol_m3: float,  # execute statement
    active_mask_global: np.ndarray,  # execute statement
) -> Tuple[Particles, float]:  # execute statement
    """Convert incremental runoff depth (m) into particles for a local slab."""  # execute statement
    depth = np.maximum(runoff_depth_m_slab, 0.0)  # set depth
    rr_local, cc = np.nonzero(depth > 0.0)  # set rr_local, cc
    if rr_local.size == 0:  # check condition rr_local.size == 0:
        return empty_particles(), 0.0  # return empty_particles(), 0.0

    rr_global = rr_local + r0  # set rr_global

    ok = active_mask_global[rr_global, cc]  # set ok
    rr_local = rr_local[ok]  # set rr_local
    rr_global = rr_global[ok]  # set rr_global
    cc = cc[ok]  # set cc
    if rr_global.size == 0:  # check condition rr_global.size == 0:
        return empty_particles(), 0.0  # return empty_particles(), 0.0

    area = cell_area_at(cell_area_m2, rr_global, cc)  # set area
    vols = depth[rr_local, cc] * area  # set vols
    total_vol = float(vols.sum())  # set total_vol

    n = np.maximum(1, np.round(vols / particle_vol_m3).astype(np.int32))  # set n

    r = np.repeat(rr_global.astype(np.int32), n)  # set r
    c = np.repeat(cc.astype(np.int32), n)  # set c
    vol = np.repeat((vols / n).astype(np.float64), n)  # set vol
    tau = np.zeros_like(vol, dtype=np.float64)  # set tau

    return Particles(r=r, c=c, vol=vol, tau=tau), total_vol  # return Particles(r=r, c=c, vol=vol, tau=tau), total_vol


def advect_particles_one_step(  # define function advect_particles_one_step
    particles: Particles,  # execute statement
    valid: np.ndarray,  # execute statement
    ds_r: np.ndarray,  # execute statement
    ds_c: np.ndarray,  # execute statement
    dt_s: float,  # execute statement
    travel_time_s: float | np.ndarray,  # execute statement
    travel_time_channel_s: float | np.ndarray,  # execute statement
    channel_mask: Optional[np.ndarray],  # execute statement
    outflow_sink: bool,  # execute statement
    shared_cfg: Optional["SharedMemoryConfig"] = None,  # execute statement
) -> Tuple[Particles, float, int]:  # execute statement
    """Advance particles one step, with travel-time gating."""  # execute statement
    n = particles.r.size  # set n
    if n == 0:  # check condition particles.r.size == 0:
        return particles, 0.0, 0  # return particles, 0.0, 0

    # Fast serial path for small workloads.
    if not should_parallelize(n, shared_cfg):  # check condition not should_parallelize(n, shared_cfg):
        particles.tau = particles.tau - dt_s  # execute statement

        can_move = (particles.tau <= 0.0)  # set can_move
        if not np.any(can_move):  # check condition not np.any(can_move):
            return particles, 0.0, 0  # return particles, 0.0, 0

        idx = np.nonzero(can_move)[0]  # set idx
        r0 = particles.r[idx]  # set r0
        c0 = particles.c[idx]  # set c0

        v = valid[r0, c0]  # set v
        nhops = int(np.count_nonzero(v))  # set nhops

        if np.any(v):  # check condition np.any(v):
            rds = ds_r[r0[v], c0[v]]  # set rds
            cds = ds_c[r0[v], c0[v]]  # set cds
            moved = idx[v]  # set moved

            particles.r[moved] = rds  # execute statement
            particles.c[moved] = cds  # execute statement

            tau_inc = _travel_time_increment(travel_time_s, travel_time_channel_s, channel_mask, rds, cds)  # set tau_inc
            particles.tau[moved] = particles.tau[moved] + tau_inc  # execute statement

        outflow_vol = 0.0  # set outflow_vol
        if outflow_sink and np.any(~v):  # check condition outflow_sink and np.any(~v):
            drop = idx[~v]  # set drop
            outflow_vol = float(particles.vol[drop].sum())  # set outflow_vol
            keep = np.ones(particles.r.shape[0], dtype=bool)  # set keep
            keep[drop] = False  # execute statement
            particles = Particles(  # set particles
                r=particles.r[keep],  # set r
                c=particles.c[keep],  # set c
                vol=particles.vol[keep],  # set vol
                tau=particles.tau[keep],  # set tau
            )  # execute statement

        return particles, outflow_vol, nhops  # return particles, outflow_vol, nhops

    # Shared-memory parallel path: work on slices, then concatenate to preserve determinism.
    from concurrent.futures import ThreadPoolExecutor  # import concurrent.futures import ThreadPoolExecutor

    # Pre-compute slices for workers.
    chunks = list(chunk_bounds(n, shared_cfg.chunk_size))  # set chunks
    # Copy once to avoid races.
    r_full = particles.r.copy()  # set r_full
    c_full = particles.c.copy()  # set c_full
    vol_full = particles.vol.copy()  # set vol_full
    tau_full = particles.tau.copy() - dt_s  # set tau_full

    def _process_chunk(start: int, end: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, float, int]:
        r_chunk = r_full[start:end].copy()
        c_chunk = c_full[start:end].copy()
        vol_chunk = vol_full[start:end].copy()
        tau_chunk = tau_full[start:end].copy()

        can_move = tau_chunk <= 0.0
        if not np.any(can_move):
            return r_chunk, c_chunk, vol_chunk, tau_chunk, 0.0, 0

        idx_local = np.nonzero(can_move)[0]
        r0 = r_chunk[idx_local]
        c0 = c_chunk[idx_local]
        v = valid[r0, c0]
        nhops_local = int(np.count_nonzero(v))

        if np.any(v):
            rds = ds_r[r0[v], c0[v]]
            cds = ds_c[r0[v], c0[v]]
            moved_local = idx_local[v]
            r_chunk[moved_local] = rds
            c_chunk[moved_local] = cds
            tau_inc = _travel_time_increment(travel_time_s, travel_time_channel_s, channel_mask, rds, cds)
            tau_chunk[moved_local] = tau_chunk[moved_local] + tau_inc

        outflow_vol_local = 0.0
        if outflow_sink and np.any(~v):
            drop = idx_local[~v]
            outflow_vol_local = float(vol_chunk[drop].sum())
            keep_mask = np.ones(r_chunk.shape[0], dtype=bool)
            keep_mask[drop] = False
            r_chunk = r_chunk[keep_mask]
            c_chunk = c_chunk[keep_mask]
            vol_chunk = vol_chunk[keep_mask]
            tau_chunk = tau_chunk[keep_mask]

        return r_chunk, c_chunk, vol_chunk, tau_chunk, outflow_vol_local, nhops_local

    results = []
    with ThreadPoolExecutor(max_workers=shared_cfg.workers) as ex:
        futures = [ex.submit(_process_chunk, s, e) for s, e in chunks]
        for fut in futures:
            results.append(fut.result())

    r_out = np.concatenate([r for r, _, _, _, _, _ in results]) if results else np.zeros(0, dtype=np.int32)
    c_out = np.concatenate([c for _, c, _, _, _, _ in results]) if results else np.zeros(0, dtype=np.int32)
    vol_out = np.concatenate([v for _, _, v, _, _, _ in results]) if results else np.zeros(0, dtype=np.float64)
    tau_out = np.concatenate([t for _, _, _, t, _, _ in results]) if results else np.zeros(0, dtype=np.float64)

    outflow_vol = float(sum(o for _, _, _, _, o, _ in results))
    nhops = int(sum(h for _, _, _, _, _, h in results))

    return Particles(r=r_out.astype(np.int32), c=c_out.astype(np.int32), vol=vol_out.astype(np.float64), tau=tau_out.astype(np.float64)), outflow_vol, nhops  # return Particles(...), outflow_vol, nhops


def local_volgrid_from_particles_slab(p: Particles, r0: int, r1: int, ncols: int, shared_cfg: Optional["SharedMemoryConfig"] = None) -> np.ndarray:  # define function local_volgrid_from_particles_slab
    """Accumulate particle volume into a local slab grid."""  # execute statement
    slab_h = r1 - r0  # set slab_h
    volgrid = np.zeros((slab_h, ncols), dtype=np.float64)  # set volgrid
    if p.r.size == 0:  # check condition p.r.size == 0:
        return volgrid  # return volgrid
    m = (p.r >= r0) & (p.r < r1)  # set m
    if not np.any(m):  # check condition not np.any(m):
        return volgrid  # return volgrid

    # Optionally parallelize accumulation for large particle counts.
    if should_parallelize(int(np.count_nonzero(m)), shared_cfg):  # check condition should_parallelize(int(np.count_nonzero(m)), shared_cfg):
        from concurrent.futures import ThreadPoolExecutor  # import concurrent.futures import ThreadPoolExecutor

        rr_all = p.r[m] - r0  # set rr_all
        cc_all = p.c[m]  # set cc_all
        vol_all = p.vol[m]  # set vol_all
        chunks = list(chunk_bounds(rr_all.size, shared_cfg.chunk_size))  # set chunks

        def _partial(start: int, end: int) -> np.ndarray:
            vg = np.zeros_like(volgrid)  # set vg
            rr = rr_all[start:end]  # set rr
            cc = cc_all[start:end]  # set cc
            np.add.at(vg, (rr, cc), vol_all[start:end])  # execute statement
            return vg  # return vg

        partials = []
        with ThreadPoolExecutor(max_workers=shared_cfg.workers) as ex:
            futures = [ex.submit(_partial, s, e) for s, e in chunks]
            for fut in futures:
                partials.append(fut.result())
        volgrid = np.sum(partials, axis=0, dtype=np.float64)  # set volgrid
        return volgrid  # return volgrid

    rr = p.r[m] - r0  # set rr
    cc = p.c[m]  # set cc
    np.add.at(volgrid, (rr, cc), p.vol[m])  # execute statement
    return volgrid  # return volgrid
