# -*- coding: utf-8 -*-
"""Runoff generation models (currently SCS Curve Number)."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import numpy.
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np  # import numpy as np

# Import optional backend helpers.
from .compute_backend import (
    get_array_module,
    normalize_device,
    to_device,
    to_numpy,
)  # import .compute_backend import get_array_module, normalize_device, to_device, to_numpy


@dataclass
class CurveNumberParams:
    """Precomputed, reusable parameters for the SCS-CN runoff model."""

    valid_mask: Any
    S: Any
    Ia: Any
    device: str
    dtype: Any


def precompute_scs_cn_params(
    CN: np.ndarray,
    ia_ratio: float,
    device: Optional[str] = None,
) -> CurveNumberParams:
    """Return precomputed CN-derived fields (S, Ia, mask) for reuse across steps."""

    xp = get_array_module(device)
    dev = normalize_device(device)
    inputs_float32 = getattr(CN, "dtype", None) == np.float32
    target_dtype = xp.float32 if inputs_float32 else xp.float64

    CNv = to_device(CN, xp).astype(target_dtype, copy=False)
    valid_mask = (CNv > 0.0) & (CNv <= 100.0) & xp.isfinite(CNv)

    S = xp.zeros_like(CNv, dtype=target_dtype)
    if bool(xp.any(valid_mask)):
        S[valid_mask] = (25400.0 / CNv[valid_mask]) - 254.0

    Ia = ia_ratio * S

    return CurveNumberParams(valid_mask=valid_mask, S=S, Ia=Ia, device=dev, dtype=target_dtype)


def scs_cn_cumulative_runoff_mm(  # define function scs_cn_cumulative_runoff_mm
    P_cum_mm: np.ndarray,  # execute statement
    CN: np.ndarray,  # execute statement
    ia_ratio: float,  # execute statement
    device: str | None = None,  # execute statement
    params: CurveNumberParams | None = None,
    workspace: np.ndarray | None = None,
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
    `workspace` allows reuse of a preallocated array for temporary storage.
    """
    dev = normalize_device(params.device if params is not None else device)
    xp = get_array_module(dev)

    if params is not None:
        if normalize_device(params.device) != dev:
            raise ValueError("precomputed CN parameters must be built for the requested device")
        S = params.S
        Ia = params.Ia
        valid_mask = params.valid_mask
        target_dtype = params.dtype
        if S.shape != P_cum_mm.shape:
            raise ValueError("precomputed CN parameters shape must match precipitation slab shape")
    else:
        inputs_float32 = getattr(P_cum_mm, "dtype", None) == np.float32 and getattr(CN, "dtype", None) == np.float32
        target_dtype = xp.float32 if inputs_float32 else xp.float64
        CNv = to_device(CN, xp).astype(target_dtype, copy=False)
        valid_mask = (CNv > 0.0) & (CNv <= 100.0) & xp.isfinite(CNv)
        S = xp.zeros_like(CNv, dtype=target_dtype)
        if bool(xp.any(valid_mask)):
            S[valid_mask] = (25400.0 / CNv[valid_mask]) - 254.0
        Ia = ia_ratio * S

    # Ensure float arrays.
    P = to_device(P_cum_mm, xp).astype(target_dtype, copy=False)  # set P

    # Initialize runoff to zero.
    Q = xp.zeros_like(P, dtype=target_dtype)  # set Q

    ok = valid_mask & xp.isfinite(P)
    if not bool(xp.any(ok)):
        return to_numpy(Q)

    cond = ok & (P > Ia) & (S > 0.0)
    if bool(xp.any(cond)):
        num = (P[cond] - Ia[cond]) ** 2
        den = (P[cond] - Ia[cond] + S[cond])
        Q[cond] = num / den

    cond = ok & (Q > 0.0) & (S > 0.0)
    if bool(xp.any(cond)):
        denom = workspace
        if denom is not None:
            if denom.shape != Q.shape:
                raise ValueError("workspace array must match precipitation shape")
            if denom.dtype != target_dtype:
                raise ValueError("workspace dtype must match target dtype")
        else:
            denom = xp.empty_like(Q, dtype=target_dtype)
        denom.fill(1.0)
        xp.add(Q, S, out=denom, where=cond, casting="unsafe")  # denom = P - Ia + S
        xp.multiply(Q, Q, out=Q, where=cond)  # Q = (P - Ia)^2
        xp.divide(Q, denom, out=Q, where=cond)
    return to_numpy(Q)  # return to_numpy(Q)
