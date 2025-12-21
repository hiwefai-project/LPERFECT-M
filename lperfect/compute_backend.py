# -*- coding: utf-8 -*-
"""Optional array backend helpers (CPU/GPU)."""

# Import importlib for optional dependency checks.
import importlib
import importlib.util

# Import typing primitives.
from typing import Any, Optional

# Import numpy.
import numpy as np

_CUPY_SPEC = importlib.util.find_spec("cupy")
if _CUPY_SPEC is not None:
    cupy = importlib.import_module("cupy")
else:
    cupy = None


def gpu_available() -> bool:
    """Return True if CuPy is available for GPU execution."""
    return cupy is not None


def normalize_device(device: Optional[str]) -> str:
    """Normalize the device string to 'cpu' or 'gpu'."""
    if device is None:
        return "cpu"
    dev = str(device).lower().strip()
    if dev not in ("cpu", "gpu"):
        raise ValueError(f"Unknown device '{device}'. Use 'cpu' or 'gpu'.")
    return dev


def get_array_module(device: Optional[str]) -> Any:
    """Return the array module (numpy or cupy) for the requested device."""
    dev = normalize_device(device)
    if dev == "gpu" and cupy is not None:
        return cupy
    return np


def to_device(arr: Any, xp: Any) -> Any:
    """Move array-like data to the selected backend."""
    if xp is np:
        return np.asarray(arr)
    return xp.asarray(arr)


def to_numpy(arr: Any) -> Any:
    """Move backend arrays to NumPy."""
    if cupy is not None and isinstance(arr, cupy.ndarray):
        return cupy.asnumpy(arr)
    return arr
