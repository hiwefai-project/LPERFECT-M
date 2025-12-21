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
    # Initialize MPI or fall back to a serial communicator.
    comm, rank, size = get_comm()

    # Parse command-line arguments into a structured namespace.
    args = parse_args()

    # Configure logging (include rank so MPI logs are distinguishable).
    setup_logging(args.log_level, rank)

    # Create a named logger for this application.
    logger = logging.getLogger("lperfect")

    # Load the built-in default configuration dictionary.
    cfg = default_config()

    # Load the user-provided config file and merge it onto the defaults.
    cfg = deep_update(cfg, load_json(args.config))

    # Apply CLI overrides (only if options were explicitly supplied).
    if args.restart_in is not None:
        # Override the restart input path.
        cfg["restart"]["in"] = args.restart_in
    if args.restart_out is not None:
        # Override the restart output path.
        cfg["restart"]["out"] = args.restart_out
    if args.out_nc is not None:
        # Override the NetCDF output path.
        cfg["output"]["out_netcdf"] = args.out_nc
    if args.device is not None:
        # Override the compute device (e.g., "cpu" or "cuda").
        cfg.setdefault("compute", {})["device"] = args.device

    # Enforce NetCDF-only domain input.
    if cfg.get("domain", {}).get("mode", "netcdf") != "netcdf":
        # Raise a clear error early if an unsupported domain mode is used.
        raise ValueError("LPERFECT is NetCDF-only: domain.mode must be 'netcdf'.")

    # Load domain (rank0) and broadcast if running under MPI.
    if size > 1:
        # In MPI, only rank 0 does I/O to avoid contention.
        if rank == 0:
            # Read the NetCDF domain file on the root rank.
            dom0 = read_domain_netcdf_rank0(cfg)
            # Log the domain shape for visibility.
            logger.info("Domain loaded (rank0): %s", dom0.dem.shape)
        else:
            # Non-root ranks hold a placeholder until broadcast.
            dom0 = None
        # Broadcast the domain object from rank 0 to all ranks.
        dom = bcast_domain(comm, dom0)  # type: ignore[arg-type]
    else:
        # In serial, read the domain directly.
        dom = read_domain_netcdf_rank0(cfg)
        # Log the domain shape for visibility.
        logger.info("Domain loaded (serial): %s", dom.dem.shape)

    # Run the main simulation driver with the loaded configuration and domain.
    run_simulation(comm, rank, size, cfg, dom)

    # Close cached NetCDF rain datasets to free resources.
    xr_close_cache()

    # Emit a final log line on rank 0.
    if rank == 0:
        # Signal successful completion.
        logger.info("LPERFECT finished.")


if __name__ == "__main__":
    main()
