# -*- coding: utf-8 -*-
"""Command line interface for LPERFECT."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import argparse for CLI parsing.
import argparse  # import argparse


def parse_args() -> argparse.Namespace:  # define function parse_args
    """Parse command-line arguments."""  # execute statement
    # Create argument parser with a program name.
    ap = argparse.ArgumentParser(prog="lperfect")  # set ap
    # Configuration file path.
    ap.add_argument("--config", default="config.json", help="Path to configuration JSON file.")  # execute statement
    # Logging level.
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])  # execute statement
    # Common operational overrides for pipelines.
    ap.add_argument("--restart-in", default=None, help="Restart NetCDF to resume from.")  # execute statement
    ap.add_argument("--restart-out", default=None, help="Restart NetCDF path to write.")  # execute statement
    ap.add_argument("--out-nc", default=None, help="Output NetCDF path (rank 0 only).")  # execute statement
    ap.add_argument("--outflow-geojson", default=None, help="GeoJSON path for outflow impact points.")  # execute statement
    ap.add_argument("--device", default=None, choices=["cpu", "gpu"], help="Compute device override.")  # execute statement
    ap.add_argument(
        "--mpi-mode",
        default=None,
        choices=["auto", "enabled", "disabled"],
        help="Force MPI on/off or auto-detect based on launcher.",
    )  # execute statement
    ap.add_argument(
        "--mpi-decomposition",
        default=None,
        choices=["auto", "balanced", "even"],
        help="MPI domain decomposition strategy (balanced=active-cell weighted slabs).",
    )  # execute statement
    ap.add_argument(
        "--mpi-min-rows",
        default=None,
        type=int,
        help="Minimum number of rows per MPI rank (prevents tiny slabs).",
    )  # execute statement
    ap.add_argument(
        "--travel-time-mode",
        default=None,
        choices=["fixed", "auto"],
        help="Travel time calculation mode override.",
    )  # execute statement
    ap.add_argument(
        "--travel-time-hill-vel",
        default=None,
        type=float,
        help="Hillslope velocity (m/s) when using --travel-time-mode auto.",
    )  # execute statement
    ap.add_argument(
        "--travel-time-channel-vel",
        default=None,
        type=float,
        help="Channel velocity (m/s) when using --travel-time-mode auto.",
    )  # execute statement
    ap.add_argument(
        "--travel-time-min",
        default=None,
        type=float,
        help="Minimum hop travel time (s) when using --travel-time-mode auto.",
    )  # execute statement
    ap.add_argument(
        "--travel-time-max",
        default=None,
        type=float,
        help="Maximum hop travel time (s) when using --travel-time-mode auto.",
    )  # execute statement
    ap.add_argument(
        "--parallel-metrics",
        action="store_true",
        help="Enable parallelization evaluation metrics (overrides metrics.parallelization.enabled).",
    )  # execute statement
    ap.add_argument(
        "--parallel-metrics-output",
        default=None,
        help="Optional JSON file to write GPT-friendly parallelization metrics.",
    )  # execute statement
    ap.add_argument(
        "--parallel-metrics-format",
        default=None,
        choices=["detailed", "compact"],
        help="JSON format for parallelization metrics (detailed=pretty-printed, compact=minified).",
    )  # execute statement
    ap.add_argument(
        "--ai-metrics",
        action="store_true",
        help="Enable GPT-friendly hydrology + compute metrics (overrides metrics.assistant.enabled).",
    )  # execute statement
    ap.add_argument(
        "--ai-metrics-output",
        default=None,
        help="Optional JSON file for AI-assistant hydrology + compute metrics.",
    )  # execute statement
    ap.add_argument(
        "--ai-metrics-format",
        default=None,
        choices=["detailed", "compact"],
        help="JSON format for AI-assistant metrics (detailed=pretty-printed, compact=minified).",
    )  # execute statement
    # Return parsed args.
    return ap.parse_args()  # return ap.parse_args()
