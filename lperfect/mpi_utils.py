# -*- coding: utf-8 -*-
"""MPI utilities for LPERFECT.

This module provides:
- MPI initialization (optional)
- slab decomposition helpers
- particle migration (Alltoallv)
- scatter/gather helpers for restart and output
"""

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import typing primitives.
from typing import Any, List, Optional, Tuple  # import typing import Any, List, Optional, Tuple

# Import dataclass for structured configs.
from dataclasses import dataclass  # import dataclasses import dataclass

# Import sys for optional early exits when MPI is disabled explicitly.
import sys  # import sys

# Import numpy for counts/displacements arrays.
import numpy as np  # import numpy as np

# Import local particle helpers.
from .particles import Particles, empty_particles, pack_particles_to_float64, unpack_particles_from_float64  # import .particles import Particles, empty_particles, pack_particles_to_float64, unpack_particles_from_float64


# Try importing mpi4py; allow serial fallback.
try:  # start exception handling
    from mpi4py import MPI  # type: ignore  # import mpi4py import MPI  # type: ignore
    HAVE_MPI = True  # set HAVE_MPI
except Exception:  # handle exception Exception:
    MPI = None  # type: ignore  # set MPI
    HAVE_MPI = False  # set HAVE_MPI


@dataclass(frozen=True)
class MPIConfig:
    """User-facing MPI configuration resolved from JSON/CLI."""

    enabled: bool
    decomposition: str
    min_rows_per_rank: int

    @classmethod
    def from_dict(cls, cfg: dict, world_size: int | None = None) -> "MPIConfig":
        """Build MPIConfig with safe defaults."""
        world = int(world_size) if world_size is not None else 1
        enabled_raw = cfg.get("enabled", None)
        enabled = bool(enabled_raw) if enabled_raw is not None else (HAVE_MPI and world > 1)
        # If mpi4py is missing, force-disable even if the user requested it.
        if enabled and not HAVE_MPI:
            enabled = False
        decomposition = str(cfg.get("decomposition", "auto") or "auto").lower()
        min_rows = max(1, int(cfg.get("min_rows_per_rank", 1) or 1))
        return cls(enabled=enabled, decomposition=decomposition, min_rows_per_rank=min_rows)


@dataclass(frozen=True)
class PartitionPlan:
    """Row-slab ownership plan shared by all ranks."""

    counts: np.ndarray
    starts: np.ndarray

    def bounds(self, rank: int) -> tuple[int, int]:
        """Return (r0, r1) bounds for a rank."""
        r0 = int(self.starts[rank])
        r1 = int(self.starts[rank] + self.counts[rank])
        return r0, r1

    def owner_of_rows(self, rows: np.ndarray) -> np.ndarray:
        """Return owning rank for each row index."""
        ends = self.starts + self.counts
        return np.searchsorted(ends, rows, side="right").astype(np.int32)

    def to_dict(self) -> dict[str, Any]:
        """Serialize the plan for metrics/debug logs."""
        return {
            "counts": [int(x) for x in self.counts.tolist()],
            "starts": [int(x) for x in self.starts.tolist()],
        }

    @property
    def size(self) -> int:
        """Return number of ranks participating in this plan."""
        return int(self.counts.size)


def get_comm(force_disabled: bool = False) -> tuple[Any, int, int]:  # define function get_comm
    """Return (comm, rank, size) for MPI or serial, honoring an explicit disable flag."""
    # If MPI is unavailable or explicitly disabled, behave as serial.
    if force_disabled or not HAVE_MPI:
        return None, 0, 1
    comm = MPI.COMM_WORLD
    return comm, comm.Get_rank(), comm.Get_size()


def initialize_mpi(mpi_cfg: MPIConfig) -> tuple[Any, int, int, int, bool]:
    """Return (comm, rank, size, world_size, active) honoring user MPI preferences."""
    if not HAVE_MPI:
        return None, 0, 1, 1, False

    world = MPI.COMM_WORLD
    world_rank = world.Get_rank()
    world_size = world.Get_size()

    # Auto-disable when only one rank is present.
    if world_size == 1:
        return None, 0, 1, 1, False

    # Respect explicit disable requests even if launched under mpirun.
    if not mpi_cfg.enabled:
        if world_rank != 0:
            # Non-root ranks exit quietly so only rank0 proceeds in serial mode.
            MPI.Finalize()
            sys.exit(0)
        return None, 0, 1, world_size, False

    return world, world_rank, world_size, world_size, True


def slab_counts_starts(nrows: int, size: int) -> tuple[np.ndarray, np.ndarray]:  # define function slab_counts_starts
    """Compute slab row counts and starts for each rank."""  # execute statement
    # Start with floor division.
    counts = np.full(size, nrows // size, dtype=np.int32)  # set counts
    # Distribute remainder to the first ranks.
    counts[: (nrows % size)] += 1  # execute statement
    # Compute starts as prefix sums of counts.
    starts = np.zeros(size, dtype=np.int32)  # set starts
    starts[1:] = np.cumsum(counts[:-1])  # execute statement
    return counts, starts  # return counts, starts


def _balanced_row_counts(weights: np.ndarray, size: int, min_rows: int) -> tuple[np.ndarray, np.ndarray]:
    """Return (counts, starts) targeting equal cumulative weight per rank."""
    nrows = int(weights.size)
    # Do not allocate more slabs than rows.
    effective_size = min(size, nrows)
    # Guard against zero weights: treat them as 1 to keep partitions non-empty.
    safe_weights = np.where(weights <= 0.0, 1.0, weights).astype(np.float64)
    total = float(safe_weights.sum())
    target = total / max(1, effective_size)
    prefix = np.cumsum(safe_weights)
    boundaries: List[int] = []
    prev_start = 0

    for i in range(effective_size - 1):
        desired = target * float(i + 1)
        cut = int(np.searchsorted(prefix, desired, side="right"))
        # Enforce minimum rows on each side of the cut.
        cut = max(cut, prev_start + min_rows)
        remaining_rows = nrows - cut
        remaining_slots = effective_size - (i + 1)
        min_required = remaining_slots * min_rows
        if remaining_rows < min_required:
            cut = nrows - min_required
        boundaries.append(cut)
        prev_start = cut

    # Build starts/counts for the effective ranks.
    starts_eff = np.zeros(effective_size, dtype=np.int32)
    for i in range(1, effective_size):
        starts_eff[i] = int(boundaries[i - 1])
    counts_eff = np.zeros(effective_size, dtype=np.int32)
    for i in range(effective_size):
        end = boundaries[i] if i < len(boundaries) else nrows
        start = int(starts_eff[i])
        counts_eff[i] = max(0, end - start)

    # Pad to full size (zero rows for unused ranks) to keep communicator alignment.
    if effective_size < size:
        pad = size - effective_size
        starts_pad = np.full(pad, nrows, dtype=np.int32)
        counts_pad = np.zeros(pad, dtype=np.int32)
        starts_eff = np.concatenate([starts_eff, starts_pad])
        counts_eff = np.concatenate([counts_eff, counts_pad])

    return counts_eff, starts_eff


def build_partition_from_active(
    active_mask: np.ndarray, size: int, strategy: str = "auto", min_rows_per_rank: int = 1
) -> PartitionPlan:
    """Build a slab partition, optionally weighted by active cells per row."""
    nrows = int(active_mask.shape[0])
    # Short-circuit serial runs.
    if size <= 1:
        return PartitionPlan(counts=np.array([nrows], dtype=np.int32), starts=np.array([0], dtype=np.int32))

    strat = (strategy or "auto").lower().strip()
    if strat in ("auto", "balanced", "weighted"):
        weights = np.sum(active_mask.astype(np.int32), axis=1)
        counts, starts = _balanced_row_counts(weights=weights, size=size, min_rows=min_rows_per_rank)
    else:
        # Even slabs remain available for reproducibility comparisons.
        counts, starts = slab_counts_starts(nrows, size)
    return PartitionPlan(counts=counts, starts=starts)


def build_partition_from_weights(
    row_weights: np.ndarray, size: int, strategy: str = "auto", min_rows_per_rank: int = 1
) -> PartitionPlan:
    """Build a slab partition from arbitrary positive row weights."""
    nrows = int(row_weights.size)
    if size <= 1:
        return PartitionPlan(counts=np.array([nrows], dtype=np.int32), starts=np.array([0], dtype=np.int32))

    strat = (strategy or "auto").lower().strip()
    weights = np.asarray(row_weights, dtype=np.float64)
    weights = np.where(weights <= 0.0, 1.0, weights)
    if strat in ("auto", "balanced", "weighted"):
        counts, starts = _balanced_row_counts(weights=weights, size=size, min_rows=min_rows_per_rank)
    else:
        counts, starts = slab_counts_starts(nrows, size)
    return PartitionPlan(counts=counts, starts=starts)


def rank_of_row_partition(rows: np.ndarray, plan: PartitionPlan) -> np.ndarray:
    """Map each row index in `rows` to its owning rank for a partition plan."""
    return plan.owner_of_rows(rows)


def alltoallv_float64(comm, sendbuf_by_rank: List[np.ndarray]) -> np.ndarray:  # define function alltoallv_float64
    """Alltoallv exchange for variable-sized float64 buffers."""  # execute statement
    # Number of ranks.
    size = comm.Get_size()  # set size
    # Build sendcounts.
    sendcounts = np.array([b.size for b in sendbuf_by_rank], dtype=np.int64)  # set sendcounts
    # Build send displacements.
    senddispls = np.zeros(size, dtype=np.int64)  # set senddispls
    senddispls[1:] = np.cumsum(sendcounts[:-1])  # execute statement
    # Flatten outgoing buffers.
    sendflat = np.concatenate(sendbuf_by_rank).astype(np.float64) if sendcounts.sum() else np.zeros(0, dtype=np.float64)  # set sendflat
    # Receive counts from other ranks.
    recvcounts = np.zeros(size, dtype=np.int64)  # set recvcounts
    comm.Alltoall(sendcounts, recvcounts)  # execute statement
    # Build receive displacements.
    recvdispls = np.zeros(size, dtype=np.int64)  # set recvdispls
    recvdispls[1:] = np.cumsum(recvcounts[:-1])  # execute statement
    # Allocate receive buffer.
    recvflat = np.empty(int(recvcounts.sum()), dtype=np.float64)  # set recvflat
    # Perform exchange.
    comm.Alltoallv(  # execute statement
        [sendflat, sendcounts, senddispls, MPI.DOUBLE],  # execute statement
        [recvflat, recvcounts, recvdispls, MPI.DOUBLE],  # execute statement
    )  # execute statement
    # Return received flat buffer.
    return recvflat  # return recvflat


def migrate_particles_partition(comm, particles: Particles, plan: PartitionPlan) -> Particles:
    """Migrate particles between ranks based on a partition plan."""
    # Rank and size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Destination rank for each particle.
    dest = rank_of_row_partition(particles.r, plan)  # set dest
    # Keep those that remain local.
    keep_mask = dest == rank  # set keep_mask
    local = Particles(
        r=particles.r[keep_mask],
        c=particles.c[keep_mask],
        vol=particles.vol[keep_mask],
        tau=particles.tau[keep_mask],
    )
    # Identify actual destinations to avoid all-to-all when only neighbors receive payloads.
    send_targets = [int(x) for x in np.unique(dest) if int(x) != rank and int(x) < size and plan.counts[int(x)] > 0]  # compute unique targets

    # Compute neighbor ranks in the slab topology (previous/next non-empty slabs).
    prev_rank = next((r for r in range(rank - 1, -1, -1) if plan.counts[r] > 0), None)  # find previous active rank
    next_rank = next((r for r in range(rank + 1, size) if plan.counts[r] > 0), None)  # find next active rank
    neighbor_set = {r for r in (prev_rank, next_rank) if r is not None}  # build neighbor set

    # Prefer lightweight neighbor exchanges when all sends target adjacent slabs.
    use_neighbor_exchange = bool(send_targets) and all(dst in neighbor_set for dst in send_targets)  # decide exchange path

    if not use_neighbor_exchange:
        # Build per-destination send buffers for the fallback all-to-all path.
        sendbuf_by_rank: List[np.ndarray] = []
        for dst in range(size):
            if dst == rank:
                sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))
                continue
            m = dest == dst
            if np.any(m):
                buf = pack_particles_to_float64(
                    Particles(r=particles.r[m], c=particles.c[m], vol=particles.vol[m], tau=particles.tau[m])
                ).ravel()
                sendbuf_by_rank.append(buf)
            else:
                sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))
        recvflat = alltoallv_float64(comm, sendbuf_by_rank)  # set recvflat
    else:
        # Build neighbor-only buffers to reduce P^2 communication.
        send_buffers = {}  # map neighbor -> flat payload
        for dst in neighbor_set:
            mask_dst = dest == dst  # mask particles bound for this neighbor
            if np.any(mask_dst):  # if any particles go to this neighbor
                send_buffers[dst] = pack_particles_to_float64(
                    Particles(r=particles.r[mask_dst], c=particles.c[mask_dst], vol=particles.vol[mask_dst], tau=particles.tau[mask_dst])
                ).ravel()  # pack payload
            else:
                send_buffers[dst] = np.zeros(0, dtype=np.float64)  # empty payload

        recv_chunks: List[np.ndarray] = []  # collect incoming payloads
        for nbr in neighbor_set:
            send_buf = send_buffers.get(nbr, np.zeros(0, dtype=np.float64))  # payload to neighbor
            send_count = np.array([send_buf.size], dtype=np.int64)  # exchange payload sizes first
            recv_count = np.zeros(1, dtype=np.int64)  # buffer for incoming size
            comm.Sendrecv(send_count, dest=nbr, sendtag=0, recvbuf=recv_count, source=nbr, recvtag=0)  # exchange sizes

            # Allocate receive buffer using the announced size.
            if int(recv_count[0]) > 0:
                recv_buf = np.empty(int(recv_count[0]), dtype=np.float64)  # allocate receive payload
            else:
                recv_buf = np.zeros(0, dtype=np.float64)  # empty receive payload

            # Exchange particle payloads.
            comm.Sendrecv(send_buf, dest=nbr, sendtag=1, recvbuf=recv_buf, source=nbr, recvtag=1)  # exchange payloads
            if recv_buf.size:  # if we received any data
                recv_chunks.append(recv_buf)  # stash chunk

        recvflat = np.concatenate(recv_chunks) if recv_chunks else np.zeros(0, dtype=np.float64)

    # If nothing received, return local only.
    if recvflat.size == 0:
        return local
    # Sanity check: payload must be multiple of 4 floats.
    if recvflat.size % 4 != 0:
        raise RuntimeError("Received particle payload is not divisible by 4")
    # Reshape and unpack.
    received = unpack_particles_from_float64(recvflat.reshape((-1, 4)))  # set received
    # Concatenate local and received.
    from .particles import concat_particles  # import .particles import concat_particles

    return concat_particles(local, received)


def scatter_field_partition(
    comm, plan: PartitionPlan, full: Optional[np.ndarray], nrows: int, ncols: int, dtype
) -> np.ndarray:
    """Scatter full (nrows,ncols) 2D array from rank0 to slabs on all ranks."""
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    counts = plan.counts.astype(np.int64)
    starts = plan.starts.astype(np.int64)
    # Local bounds.
    r0, r1 = plan.bounds(rank)  # set r0, r1
    slab_h = r1 - r0  # set slab_h
    # Element counts/displacements for the flattened arrays.
    sendcounts = (counts * ncols).astype(np.int64)  # set sendcounts
    displs = (starts * ncols).astype(np.int64)  # set displs
    # Allocate local slab.
    local = np.empty((slab_h, ncols), dtype=dtype)  # set local
    # Prepare rank0 send buffer.
    if rank == 0:
        if full is None:
            # Initialize zeros when no full-field state is provided (e.g., fresh run with no restart).
            sendbuf = np.zeros((nrows, ncols), dtype=dtype).ravel()
        else:
            arr = np.asarray(full)
            if arr.shape != (nrows, ncols):
                raise ValueError(f"Scatter source has shape {arr.shape}, expected {(nrows, ncols)}")
            sendbuf = arr.astype(dtype, copy=False).ravel()
    else:
        sendbuf = None
    # Map dtype to MPI datatype.
    mpitype = MPI._typedict[np.dtype(dtype).char]  # set mpitype
    # Scatter.
    comm.Scatterv([sendbuf, sendcounts, displs, mpitype], local.ravel(), root=0)  # execute statement
    # Return local slab.
    return local  # return local


def gather_field_partition_to_rank0(comm, plan: PartitionPlan, slab: np.ndarray, nrows: int, ncols: int) -> Optional[np.ndarray]:
    """Gather slabs from all ranks to a full array on rank0."""
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    # Slab counts/starts.
    counts = plan.counts.astype(np.int64)
    starts = plan.starts.astype(np.int64)
    # Receive counts/displacements.
    recvcounts = (counts * ncols).astype(np.int64)  # set recvcounts
    displs = (starts * ncols).astype(np.int64)  # set displs
    # Allocate on rank0 only.
    full_flat = np.empty((nrows * ncols,), dtype=slab.dtype) if rank == 0 else None  # set full_flat
    # MPI datatype.
    mpitype = MPI._typedict[np.dtype(slab.dtype).char]  # set mpitype
    # Gather.
    comm.Gatherv(slab.ravel(), [full_flat, recvcounts, displs, mpitype], root=0)  # execute statement
    # Return full array on rank0.
    if rank != 0:
        return None
    return full_flat.reshape((nrows, ncols))


def gather_particles_to_rank0(comm, p_local: Particles) -> Particles:  # define function gather_particles_to_rank0
    """Gather particles from all ranks to rank0 (using Alltoallv pattern)."""  # execute statement
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Pack local particles.
    buf = pack_particles_to_float64(p_local).ravel()  # set buf
    # Build send list: only destination rank 0 gets the payload.
    sendbuf_by_rank = [buf if dst == 0 else np.zeros(0, dtype=np.float64) for dst in range(size)]  # set sendbuf_by_rank
    # Exchange.
    recvflat = alltoallv_float64(comm, sendbuf_by_rank)  # set recvflat
    # Non-root returns empty.
    if rank != 0 or recvflat.size == 0:  # check condition rank != 0 or recvflat.size == 0:
        return empty_particles()  # return empty_particles()
    # Unpack on root.
    return unpack_particles_from_float64(recvflat.reshape((-1, 4)))  # return unpack_particles_from_float64(recvflat.reshape((-1, 4)))


def scatter_particles_from_rank0(comm, plan: PartitionPlan, p_all: Optional[Particles]) -> Particles:
    """Scatter particles from rank0 to owning ranks based on row ownership."""
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Prepare send buffers.
    if rank == 0 and p_all is not None and p_all.r.size > 0:
        dest = rank_of_row_partition(p_all.r, plan)
        sendbuf_by_rank: List[np.ndarray] = []
        for dst in range(size):
            m = dest == dst
            if np.any(m):
                sendbuf_by_rank.append(
                    pack_particles_to_float64(Particles(r=p_all.r[m], c=p_all.c[m], vol=p_all.vol[m], tau=p_all.tau[m])).ravel()
                )
            else:
                sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))
    else:
        sendbuf_by_rank = [np.zeros(0, dtype=np.float64) for _ in range(size)]
    recvflat = alltoallv_float64(comm, sendbuf_by_rank)
    if recvflat.size == 0:
        return empty_particles()
    if recvflat.size % 4 != 0:
        raise RuntimeError("Received particle payload not divisible by 4")
    return unpack_particles_from_float64(recvflat.reshape((-1, 4)))



def rebalance_particles_even(comm, p_local: Particles) -> Particles:
    """Evenly redistribute particles across ranks (count-balanced, location-agnostic)."""
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size

    local_count = int(p_local.r.size)  # set local_count
    total_count = int(comm.allreduce(local_count, op=MPI.SUM))  # set total_count
    if total_count == 0:  # check condition total_count == 0:
        return empty_particles()  # return empty_particles()

    particles_all = gather_particles_to_rank0(comm, p_local)  # collect full set on rank0

    counts = np.full(size, total_count // max(1, size), dtype=np.int64)  # set counts
    counts[: (total_count % max(1, size))] += 1  # execute statement

    if rank == 0:  # check condition rank == 0:
        flat = pack_particles_to_float64(particles_all).ravel()  # set flat
        sendbuf_by_rank: list[np.ndarray] = []  # set sendbuf_by_rank
        offset = 0  # set offset
        for c in counts:  # loop over c in counts:
            nvals = int(c) * 4  # set nvals
            if nvals > 0:  # check condition nvals > 0:
                sendbuf_by_rank.append(flat[offset : offset + nvals])  # execute statement
            else:  # fallback branch
                sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))  # execute statement
            offset += nvals  # execute statement
    else:  # fallback branch
        sendbuf_by_rank = [np.zeros(0, dtype=np.float64) for _ in range(size)]  # set sendbuf_by_rank

    recvflat = alltoallv_float64(comm, sendbuf_by_rank)  # exchange evenly sized chunks
    if recvflat.size == 0:  # check condition recvflat.size == 0:
        return empty_particles()  # return empty_particles()
    if recvflat.size % 4 != 0:  # check condition recvflat.size % 4 != 0:
        raise RuntimeError("Even rebalance payload not divisible by 4")  # raise RuntimeError

    return unpack_particles_from_float64(recvflat.reshape((-1, 4)))  # unpack received particles
