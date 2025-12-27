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


def get_comm() -> tuple[Any, int, int]:  # define function get_comm
    """Return (comm, rank, size) for MPI or serial."""  # execute statement
    # If mpi4py is available, use COMM_WORLD.
    if HAVE_MPI:  # check condition HAVE_MPI:
        comm = MPI.COMM_WORLD  # set comm
        return comm, comm.Get_rank(), comm.Get_size()  # return comm, comm.Get_rank(), comm.Get_size()
    # Otherwise emulate a single-rank communicator as (None,0,1).
    return None, 0, 1  # return None, 0, 1


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


def slab_bounds(nrows: int, size: int, rank: int) -> tuple[int, int]:  # define function slab_bounds
    """Return (r0,r1) bounds of the slab owned by `rank`."""  # execute statement
    # Get counts and starts.
    counts, starts = slab_counts_starts(nrows, size)  # set counts, starts
    # Start row.
    r0 = int(starts[rank])  # set r0
    # End row (exclusive).
    r1 = int(starts[rank] + counts[rank])  # set r1
    return r0, r1  # return r0, r1


def rank_of_row(r: np.ndarray, nrows: int, size: int) -> np.ndarray:  # define function rank_of_row
    """Map each row index in r to its owning rank."""  # execute statement
    # Build ends array (start+count for each rank).
    counts, starts = slab_counts_starts(nrows, size)  # set counts, starts
    ends = starts + counts  # set ends
    # searchsorted finds the first end > r.
    return np.searchsorted(ends, r, side="right").astype(np.int32)  # return np.searchsorted(ends, r, side="right").astype(np.int32)


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


def migrate_particles_slab(comm, particles: Particles, nrows: int) -> Particles:  # define function migrate_particles_slab
    """Migrate particles between ranks based on slab ownership."""  # execute statement
    # Rank and size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Destination rank for each particle.
    dest = rank_of_row(particles.r, nrows, size)  # set dest
    # Keep those that remain local.
    keep_mask = (dest == rank)  # set keep_mask
    local = Particles(  # set local
        r=particles.r[keep_mask],  # set r
        c=particles.c[keep_mask],  # set c
        vol=particles.vol[keep_mask],  # set vol
        tau=particles.tau[keep_mask],  # set tau
    )  # execute statement
    # Build per-destination send buffers.
    sendbuf_by_rank: List[np.ndarray] = []  # execute statement
    for dst in range(size):  # loop over dst in range(size):
        # Skip self-destination.
        if dst == rank:  # check condition dst == rank:
            sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))  # execute statement
            continue  # continue loop
        # Mask for particles going to dst.
        m = (dest == dst)  # set m
        if np.any(m):  # check condition np.any(m):
            # Pack to float64 and flatten.
            buf = pack_particles_to_float64(Particles(  # set buf
                r=particles.r[m], c=particles.c[m], vol=particles.vol[m], tau=particles.tau[m]  # set r
            )).ravel()  # execute statement
            sendbuf_by_rank.append(buf)  # execute statement
        else:  # fallback branch
            sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))  # execute statement
    # Exchange buffers.
    recvflat = alltoallv_float64(comm, sendbuf_by_rank)  # set recvflat
    # If nothing received, return local only.
    if recvflat.size == 0:  # check condition recvflat.size == 0:
        return local  # return local
    # Sanity check: payload must be multiple of 4 floats.
    if recvflat.size % 4 != 0:  # check condition recvflat.size % 4 != 0:
        raise RuntimeError("Received particle payload is not divisible by 4")  # raise RuntimeError("Received particle payload is not divisible by 4")
    # Reshape and unpack.
    received = unpack_particles_from_float64(recvflat.reshape((-1, 4)))  # set received
    # Concatenate local and received.
    from .particles import concat_particles  # import .particles import concat_particles
    return concat_particles(local, received)  # return concat_particles(local, received)


def scatter_field_slab(comm, full: Optional[np.ndarray], nrows: int, ncols: int, dtype) -> np.ndarray:  # define function scatter_field_slab
    """Scatter full (nrows,ncols) 2D array from rank0 to slabs on all ranks."""  # execute statement
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Slab counts/starts.
    counts, starts = slab_counts_starts(nrows, size)  # set counts, starts
    # Local bounds.
    r0, r1 = slab_bounds(nrows, size, rank)  # set r0, r1
    slab_h = r1 - r0  # set slab_h
    # Element counts/displacements for the flattened arrays.
    sendcounts = (counts.astype(np.int64) * ncols).astype(np.int64)  # set sendcounts
    displs = (starts.astype(np.int64) * ncols).astype(np.int64)  # set displs
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


def gather_field_slab_to_rank0(comm, slab: np.ndarray, nrows: int, ncols: int) -> Optional[np.ndarray]:  # define function gather_field_slab_to_rank0
    """Gather slabs from all ranks to a full array on rank0."""  # execute statement
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Slab counts/starts.
    counts, starts = slab_counts_starts(nrows, size)  # set counts, starts
    # Receive counts/displacements.
    recvcounts = (counts.astype(np.int64) * ncols).astype(np.int64)  # set recvcounts
    displs = (starts.astype(np.int64) * ncols).astype(np.int64)  # set displs
    # Allocate on rank0 only.
    full_flat = np.empty((nrows * ncols,), dtype=slab.dtype) if rank == 0 else None  # set full_flat
    # MPI datatype.
    mpitype = MPI._typedict[np.dtype(slab.dtype).char]  # set mpitype
    # Gather.
    comm.Gatherv(slab.ravel(), [full_flat, recvcounts, displs, mpitype], root=0)  # execute statement
    # Return full array on rank0.
    if rank != 0:  # check condition rank != 0:
        return None  # return None
    return full_flat.reshape((nrows, ncols))  # return full_flat.reshape((nrows, ncols))


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


def scatter_particles_from_rank0(comm, p_all: Optional[Particles], nrows: int) -> Particles:  # define function scatter_particles_from_rank0
    """Scatter particles from rank0 to owning ranks based on row ownership."""  # execute statement
    # Rank/size.
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size
    # Prepare send buffers.
    if rank == 0 and p_all is not None and p_all.r.size > 0:  # check condition rank == 0 and p_all is not None and p_all.r.size > 0:
        # Compute destination for each particle.
        dest = rank_of_row(p_all.r, nrows, size)  # set dest
        sendbuf_by_rank: List[np.ndarray] = []  # execute statement
        # Build payload per destination.
        for dst in range(size):  # loop over dst in range(size):
            m = (dest == dst)  # set m
            if np.any(m):  # check condition np.any(m):
                sendbuf_by_rank.append(pack_particles_to_float64(Particles(  # execute statement
                    r=p_all.r[m], c=p_all.c[m], vol=p_all.vol[m], tau=p_all.tau[m]  # set r
                )).ravel())  # execute statement
            else:  # fallback branch
                sendbuf_by_rank.append(np.zeros(0, dtype=np.float64))  # execute statement
    else:  # fallback branch
        # Non-root ranks send empty buffers.
        sendbuf_by_rank = [np.zeros(0, dtype=np.float64) for _ in range(size)]  # set sendbuf_by_rank
    # Exchange.
    recvflat = alltoallv_float64(comm, sendbuf_by_rank)  # set recvflat
    # Unpack.
    if recvflat.size == 0:  # check condition recvflat.size == 0:
        return empty_particles()  # return empty_particles()
    if recvflat.size % 4 != 0:  # check condition recvflat.size % 4 != 0:
        raise RuntimeError("Received particle payload not divisible by 4")  # raise RuntimeError("Received particle payload not divisible by 4")
    return unpack_particles_from_float64(recvflat.reshape((-1, 4)))  # return unpack_particles_from_float64(recvflat.reshape((-1, 4)))
