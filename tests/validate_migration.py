# -*- coding: utf-8 -*-
"""Sanity check for MPI particle migration modes."""  # execute statement

# Import MPI for distributed test harness.
from mpi4py import MPI  # import mpi4py import MPI

# Import numpy for array handling.
import numpy as np  # import numpy as np

# Import local particle and MPI helpers.
from lperfect.particles import Particles  # import .particles import Particles
from lperfect.mpi_utils import PartitionPlan, gather_particles_to_rank0, migrate_particles_partition, migrate_particles_partition_agg  # import .mpi_utils import PartitionPlan, gather_particles_to_rank0, migrate_particles_partition, migrate_particles_partition_agg


def _build_particles(rank: int) -> Particles:
    """Create a tiny deterministic particle set with cross-rank migrants."""  # execute statement
    if rank == 0:  # check condition rank == 0:
        r = np.array([0, 1, 3], dtype=np.int32)  # set r
        c = np.array([0, 1, 0], dtype=np.int32)  # set c
        vol = np.array([1.0, 2.0, 3.0], dtype=np.float64)  # set vol
        tau = np.array([0.0, 0.0, 0.0], dtype=np.float64)  # set tau
    else:  # fallback branch
        r = np.array([2, 1], dtype=np.int32)  # set r
        c = np.array([1, 0], dtype=np.int32)  # set c
        vol = np.array([4.0, 5.0], dtype=np.float64)  # set vol
        tau = np.array([0.0, 0.0], dtype=np.float64)  # set tau
    return Particles(r=r, c=c, vol=vol, tau=tau)  # return Particles(...)


def _sorted_view(p: Particles) -> np.ndarray:
    """Return a sortable view of particle fields for comparison."""  # execute statement
    if p.r.size == 0:  # check condition p.r.size == 0:
        return np.zeros((0, 4), dtype=np.float64)  # return empty
    stacked = np.column_stack([p.r.astype(np.float64), p.c.astype(np.float64), p.vol, p.tau])  # set stacked
    order = np.lexsort((stacked[:, 3], stacked[:, 2], stacked[:, 1], stacked[:, 0]))  # set order
    return stacked[order]  # return stacked[order]


def main() -> None:
    """Run migration checks for legacy and aggregated modes."""  # execute statement
    comm = MPI.COMM_WORLD  # set comm
    rank = comm.Get_rank()  # set rank
    size = comm.Get_size()  # set size

    # Build a simple 4-row partition for 1-2 ranks.
    nrows = 4  # set nrows
    counts = np.array([2, 2] if size > 1 else [nrows], dtype=np.int32)  # set counts
    starts = np.array([0, 2] if size > 1 else [0], dtype=np.int32)  # set starts
    plan = PartitionPlan(counts=counts, starts=starts)  # set plan

    # Baseline particle set on each rank.
    base = _build_particles(rank)  # set base
    total_before = int(comm.allreduce(base.r.size, op=MPI.SUM))  # set total_before

    # Run legacy migration.
    legacy_local, _ = migrate_particles_partition(comm, base, plan=plan)  # migrate and count
    legacy_total = int(comm.allreduce(legacy_local.r.size, op=MPI.SUM))  # set legacy_total
    if legacy_total != total_before:  # check condition legacy_total != total_before:
        raise RuntimeError(f"Legacy migration count mismatch: {legacy_total} vs {total_before}")  # raise RuntimeError

    # Run aggregated nonblocking migration.
    base2 = _build_particles(rank)  # reset base
    agg_local, _, _ = migrate_particles_partition_agg(comm, base2, plan=plan, overlap=True)  # migrate and count
    agg_total = int(comm.allreduce(agg_local.r.size, op=MPI.SUM))  # set agg_total
    if agg_total != total_before:  # check condition agg_total != total_before:
        raise RuntimeError(f"Agg migration count mismatch: {agg_total} vs {total_before}")  # raise RuntimeError

    # Compare full particle sets on rank0.
    legacy_all = gather_particles_to_rank0(comm, legacy_local)  # gather legacy
    agg_all = gather_particles_to_rank0(comm, agg_local)  # gather agg
    if rank == 0:  # check condition rank == 0:
        legacy_view = _sorted_view(legacy_all)  # set legacy_view
        agg_view = _sorted_view(agg_all)  # set agg_view
        if legacy_view.shape != agg_view.shape:  # check condition shapes mismatch
            raise RuntimeError("Aggregated migration output shape mismatch with legacy output")  # raise RuntimeError
        if not np.allclose(legacy_view, agg_view, atol=0.0, rtol=0.0):  # check condition not allclose
            raise RuntimeError("Aggregated migration output does not match legacy output")  # raise RuntimeError
        print("Migration sanity checks passed.")  # execute statement


if __name__ == "__main__":  # check condition __name__ == "__main__":
    main()  # execute statement
