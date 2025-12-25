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

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import logging (for module-level logger).
import logging

# Import stdlib helpers.
import importlib.util
import sys
from copy import deepcopy
from pathlib import Path
import re
from typing import Any, Dict, List

# Import lightweight config helpers early for shared utilities.
from lperfect.config import deep_update, default_config, load_json


def _require_numpy() -> None:
    """Validate that NumPy is available before importing LPERFECT modules."""
    if importlib.util.find_spec("numpy") is None:
        raise ModuleNotFoundError(
            "NumPy is required to run LPERFECT. Activate your virtual environment "
            f"or install it with '{sys.executable} -m pip install numpy'."
        )


def _slugify(label: str) -> str:
    """Return a filesystem-friendly slug for a domain label."""
    slug = re.sub(r"[^A-Za-z0-9]+", "-", label).strip("-").lower()
    return slug or "domain"


def _tag_path(path: str | None, label: str) -> str | None:
    """Append a domain label suffix to a path while preserving extension."""
    if not path:
        return None
    p = Path(str(path))
    suffix = p.suffix or ""
    return str(p.with_name(f"{p.stem}_{label}{suffix}"))


def _normalize_domain_entries(cfg: Dict[str, Any], default_domain: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Return a list of domain configuration dictionaries."""
    domains_cfg = cfg.get("domains", None)
    if domains_cfg is None:
        domains_cfg = default_domain

    if isinstance(domains_cfg, dict):
        domains_list = [domains_cfg]
    elif isinstance(domains_cfg, list):
        domains_list = domains_cfg
    else:
        raise ValueError("domain/domains config must be an object or list of objects")

    resolved: List[Dict[str, Any]] = []
    for idx, dom_cfg in enumerate(domains_list):
        merged = deepcopy(default_domain)
        if not isinstance(dom_cfg, dict):
            raise ValueError("Each domain entry must be an object")
        merged = deep_update(merged, dom_cfg)
        name = merged.get("name") or merged.get("id") or f"domain_{idx + 1}"
        merged.setdefault("name", name)
        resolved.append(merged)
    return resolved


def _prepare_domain_run_config(
    base_cfg: Dict[str, Any],
    dom_cfg: Dict[str, Any],
    label_slug: str,
    ndomains: int,
) -> Dict[str, Any]:
    """Clone the base config and apply domain-specific overrides and tagging."""

    run_cfg = deepcopy(base_cfg)
    run_cfg["domain"] = dom_cfg
    run_cfg.pop("domains", None)

    # Allow per-domain overrides for output/restart while keeping defaults.
    run_cfg["output"] = deep_update(run_cfg.get("output", {}), dom_cfg.get("output", {}))
    run_cfg["restart"] = deep_update(run_cfg.get("restart", {}), dom_cfg.get("restart", {}))

    # Enforce distinct artifacts when running multiple domains.
    if ndomains > 1:
        if "output" not in dom_cfg or "out_netcdf" not in dom_cfg.get("output", {}):
            run_cfg["output"]["out_netcdf"] = _tag_path(run_cfg["output"].get("out_netcdf"), label_slug)
        if "restart" not in dom_cfg or "out" not in dom_cfg.get("restart", {}):
            run_cfg["restart"]["out"] = _tag_path(run_cfg["restart"].get("out"), label_slug)
        if "restart" not in dom_cfg or "in" not in dom_cfg.get("restart", {}):
            run_cfg["restart"]["in"] = _tag_path(run_cfg["restart"].get("in"), label_slug)
        if "output" not in dom_cfg or "outflow_geojson" not in dom_cfg.get("output", {}):
            run_cfg["output"]["outflow_geojson"] = _tag_path(run_cfg["output"].get("outflow_geojson"), label_slug)

    return run_cfg


def main() -> None:
    """Program entry point."""
    _require_numpy()

    # Import CLI parser.
    from lperfect.cli import parse_args

    # Import logging configuration.
    from lperfect.logging_utils import setup_logging

    # Import MPI utilities.
    from lperfect.mpi_utils import get_comm

    # Import domain I/O.
    from lperfect.domain import read_domain_netcdf_rank0, bcast_domain

    # Import simulation driver.
    from lperfect.simulation import run_simulation, run_nested_simulations

    # Import rain cache close.
    from lperfect.rain import xr_close_cache
    # Import shared-memory config helper.
    from lperfect.shared_memory import SharedMemoryConfig
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
    if args.outflow_geojson is not None:
        # Override outflow GeoJSON export path.
        cfg["output"]["outflow_geojson"] = args.outflow_geojson
    if args.device is not None:
        # Override the compute device (e.g., "cpu" or "cuda").
        cfg.setdefault("compute", {})["device"] = args.device
    if args.travel_time_mode is not None:
        cfg["model"]["travel_time_mode"] = args.travel_time_mode
    if args.travel_time_hill_vel is not None:
        cfg["model"].setdefault("travel_time_auto", {})["hillslope_velocity_ms"] = args.travel_time_hill_vel
    if args.travel_time_channel_vel is not None:
        cfg["model"].setdefault("travel_time_auto", {})["channel_velocity_ms"] = args.travel_time_channel_vel
    if args.travel_time_min is not None:
        cfg["model"].setdefault("travel_time_auto", {})["min_s"] = args.travel_time_min
    if args.travel_time_max is not None:
        cfg["model"].setdefault("travel_time_auto", {})["max_s"] = args.travel_time_max

    # Normalize domains (support single domain object or list under cfg['domains']).
    default_domain_cfg = cfg.get("domain", {})
    domains = _normalize_domain_entries(cfg, default_domain_cfg)
    ndomains = len(domains)
    if ndomains == 0:
        raise ValueError("At least one domain must be configured.")

    for idx, dom_cfg in enumerate(domains):
        domain_label = str(dom_cfg.get("name", f"domain_{idx + 1}"))
        label_slug = _slugify(domain_label)
        run_cfg = _prepare_domain_run_config(cfg, dom_cfg, label_slug, ndomains)

        # Enforce NetCDF-only domain input.
        if run_cfg.get("domain", {}).get("mode", "netcdf") != "netcdf":
            raise ValueError("LPERFECT is NetCDF-only: domain.mode must be 'netcdf'.")

        # Load domain (rank0) and broadcast if running under MPI.
        if size > 1:
            # In MPI, only rank 0 does I/O to avoid contention.
            if rank == 0:
                dom0 = read_domain_netcdf_rank0(run_cfg)
                logger.info("Domain '%s' loaded (rank0): %s", domain_label, dom0.dem.shape)
            else:
                dom0 = None
            dom = bcast_domain(comm, dom0)  # type: ignore[arg-type]
        else:
            dom = read_domain_netcdf_rank0(run_cfg)
            logger.info("Domain '%s' loaded (serial): %s", domain_label, dom.dem.shape)

        if rank == 0:
            logger.info(
                "Starting simulation for domain '%s' (%d/%d)",
                domain_label,
                idx + 1,
                ndomains,
            )

        # Run the main simulation driver with the loaded configuration and domain.
        run_simulation(comm, rank, size, run_cfg, dom, domain_label=domain_label)

        # Close cached NetCDF rain datasets to free resources between domains.
        xr_close_cache()

    # Emit a final log line on rank 0.
    if rank == 0:
        # Signal successful completion.
        logger.info("LPERFECT finished across %d domain(s).", ndomains)


if __name__ == "__main__":
    main()
