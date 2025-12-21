# -*- coding: utf-8 -*-
"""Runoff generation models (currently SCS Curve Number)."""

# Import numpy.
import numpy as np

# Import optional backend helpers.
from .compute_backend import get_array_module, to_device, to_numpy


def scs_cn_cumulative_runoff_mm(
    P_cum_mm: np.ndarray,
    CN: np.ndarray,
    ia_ratio: float,
    device: str | None = None,
) -> np.ndarray:
    """Compute cumulative runoff Q (mm) from cumulative precipitation P (mm) using SCS-CN.

    Formulae
    --------
    S  = 25400/CN - 254  [mm]
    Ia = ia_ratio * S     [mm]

    Q = 0                              if P <= Ia
    Q = (P - Ia)^2 / (P - Ia + S)       if P > Ia

    Notes
    -----
    CN must be in (0,100]. Invalid CN yields Q=0.
    """
    # Choose backend.
    xp = get_array_module(device)

    # Ensure float arrays.
    P = to_device(P_cum_mm, xp).astype(xp.float64)
    CNv = to_device(CN, xp).astype(xp.float64)

    # Initialize runoff to zero.
    Q = xp.zeros_like(P, dtype=xp.float64)

    # Mask valid cells.
    ok = (CNv > 0.0) & (CNv <= 100.0) & xp.isfinite(CNv) & xp.isfinite(P)
    if not bool(xp.any(ok)):
        return to_numpy(Q)

    # Potential retention S.
    S = xp.zeros_like(P, dtype=xp.float64)
    S[ok] = (25400.0 / CNv[ok]) - 254.0

    # Initial abstraction.
    Ia = ia_ratio * S

    # Condition for runoff.
    cond = ok & (P > Ia) & (S > 0.0)
    if bool(xp.any(cond)):
        num = (P[cond] - Ia[cond]) ** 2
        den = (P[cond] - Ia[cond] + S[cond])
        Q[cond] = num / den

    return to_numpy(Q)
