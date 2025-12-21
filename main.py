#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""LPERFECT entry point.

This file is intentionally small:
- parse CLI
- load+merge configuration
- initialize MPI (optional)
- load+broadcast domain
- run the simulation

All real logic lives in the `lperfect/` package.
"""

# Import logging (for module-level logger).
import logging

# Import config helpers.
from lperfect.config import default_config, load_json, deep_update

# Import CLI parser.
from lperfect.cli import parse_args

# Import logging configuration.
from lperfect.logging_utils import setup_logging

# Import MPI utilities.
from lperfect.mpi_utils import get_comm

# Import domain I/O.
from lperfect.domain import read_domain_netcdf_rank0, bcast_domain

# Import simulation driver.
from lperfect.simulation import run_simulation

# Import rain cache close.
from lperfect.rain import xr_close_cache


def main() -> None:
    """Program entry point."""
    # Initialize MPI or serial communicator.
    comm, rank, size = get_comm()

    # Parse command line.
    args = parse_args()

    # Configure logging with rank info.
    setup_logging(args.log_level, rank)

    # Create a logger.
    logger = logging.getLogger("lperfect")

    # Load defaults.
    cfg = default_config()

    # Merge config file on top of defaults.
    cfg = deep_update(cfg, load_json(args.config))

    # Apply CLI overrides (if provided).
    if args.restart_in is not None:
        cfg["restart"]["in"] = args.restart_in
    if args.restart_out is not None:
        cfg["restart"]["out"] = args.restart_out
    if args.out_nc is not None:
        cfg["output"]["out_netcdf"] = args.out_nc

    # Enforce NetCDF-only.
    if cfg.get("domain", {}).get("mode", "netcdf") != "netcdf":
        raise ValueError("LPERFECT is NetCDF-only: domain.mode must be 'netcdf'.")

    # Load domain (rank0) and broadcast if MPI.
    if size > 1:
        if rank == 0:
            dom0 = read_domain_netcdf_rank0(cfg)
            logger.info("Domain loaded (rank0): %s", dom0.dem.shape)
        else:
            dom0 = None
        dom = bcast_domain(comm, dom0)  # type: ignore[arg-type]
    else:
        dom = read_domain_netcdf_rank0(cfg)
        logger.info("Domain loaded (serial): %s", dom.dem.shape)

    # Run simulation.
    run_simulation(comm, rank, size, cfg, dom)

    # Close cached NetCDF rain datasets.
    xr_close_cache()

    # Final log on root.
    if rank == 0:
        logger.info("LPERFECT finished.")


if __name__ == "__main__":
    main()
