# -*- coding: utf-8 -*-
"""Shared-memory parallel helpers (thread-based)."""

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

from __future__ import annotations

# Import dataclass for structured config.
from dataclasses import dataclass

# Import typing primitives.
from typing import List, Tuple

# Import stdlib helpers.
import math
import os


@dataclass(frozen=True)
class SharedMemoryConfig:
    """Configuration for optional shared-memory parallelism."""

    enabled: bool
    workers: int
    min_particles_per_worker: int
    chunk_size: int

    @classmethod
    def from_dict(cls, cfg: dict) -> "SharedMemoryConfig":
        """Build config from raw dictionary."""
        enabled = bool(cfg.get("enabled", False))
        raw_workers = cfg.get("workers", None)
        # Default: use hardware concurrency if available.
        workers = int(raw_workers) if raw_workers not in (None, "") else max(1, os.cpu_count() or 1)
        workers = max(1, workers)
        min_particles = int(cfg.get("min_particles_per_worker", 5000))
        chunk_size = int(cfg.get("chunk_size", 65536))
        return cls(
            enabled=enabled,
            workers=workers,
            min_particles_per_worker=max(1, min_particles),
            chunk_size=max(1, chunk_size),
        )


def should_parallelize(n_items: int, cfg: SharedMemoryConfig | None) -> bool:
    """Return True when shared-memory parallelism should be used."""
    return (
        cfg is not None
        and cfg.enabled
        and cfg.workers > 1
        and n_items >= cfg.min_particles_per_worker
    )


def chunk_bounds(n_items: int, chunk_size: int) -> List[Tuple[int, int]]:
    """Return list of (start, end) chunk bounds covering [0, n_items)."""
    if n_items <= 0:
        return []
    if chunk_size <= 0:
        return [(0, n_items)]
    nchunks = int(math.ceil(n_items / float(chunk_size)))
    out: List[Tuple[int, int]] = []
    for i in range(nchunks):
        start = i * chunk_size
        end = min(n_items, start + chunk_size)
        out.append((start, end))
    return out
