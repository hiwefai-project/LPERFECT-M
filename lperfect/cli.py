# -*- coding: utf-8 -*-
"""Command line interface for LPERFECT."""

# Import argparse for CLI parsing.
import argparse


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    # Create argument parser with a program name.
    ap = argparse.ArgumentParser(prog="lperfect")
    # Configuration file path.
    ap.add_argument("--config", default="config.json", help="Path to configuration JSON file.")
    # Logging level.
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    # Common operational overrides for pipelines.
    ap.add_argument("--restart-in", default=None, help="Restart NetCDF to resume from.")
    ap.add_argument("--restart-out", default=None, help="Restart NetCDF path to write.")
    ap.add_argument("--out-nc", default=None, help="Output NetCDF path (rank 0 only).")
    ap.add_argument("--device", default=None, choices=["cpu", "gpu"], help="Compute device override.")
    # Return parsed args.
    return ap.parse_args()
