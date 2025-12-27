# -*- coding: utf-8 -*-
"""Runoff generation models (currently SCS Curve Number)."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import numpy.
import numpy as np  # import numpy as np

# Import optional backend helpers.
from .compute_backend import get_array_module, to_device, to_numpy  # import .compute_backend import get_array_module, to_device, to_numpy


def scs_cn_cumulative_runoff_mm(  # define function scs_cn_cumulative_runoff_mm
    P_cum_mm: np.ndarray,  # execute statement
    CN: np.ndarray,  # execute statement
    ia_ratio: float,  # execute statement
    device: str | None = None,  # execute statement
) -> np.ndarray:  # execute statement
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
    # Choose backend and dtype (float32 to reduce memory when inputs are float32).
    xp = get_array_module(device)  # set xp
    inputs_float32 = getattr(P_cum_mm, "dtype", None) == np.float32 and getattr(CN, "dtype", None) == np.float32
    target_dtype = xp.float32 if inputs_float32 else xp.float64

    # Ensure float arrays.
    P = to_device(P_cum_mm, xp).astype(target_dtype, copy=False)  # set P
    CNv = to_device(CN, xp).astype(target_dtype, copy=False)  # set CNv

    # Initialize runoff to zero.
    Q = xp.zeros_like(P, dtype=target_dtype)  # set Q

    # Mask valid cells.
    ok = (CNv > 0.0) & (CNv <= 100.0) & xp.isfinite(CNv) & xp.isfinite(P)  # set ok
    if not bool(xp.any(ok)):  # check condition not bool(xp.any(ok)):
        return to_numpy(Q)  # return to_numpy(Q)

    # Potential retention S.
    S = xp.zeros_like(P, dtype=target_dtype)  # set S
    S[ok] = (25400.0 / CNv[ok]) - 254.0  # execute statement

    # Initial abstraction.
    Ia = ia_ratio * S  # set Ia

    # Condition for runoff.
    cond = ok & (P > Ia) & (S > 0.0)  # set cond
    if bool(xp.any(cond)):  # check condition bool(xp.any(cond)):
        num = (P[cond] - Ia[cond]) ** 2  # set num
        den = (P[cond] - Ia[cond] + S[cond])  # set den
        Q[cond] = num / den  # execute statement

    return to_numpy(Q)  # return to_numpy(Q)
