# -*- coding: utf-8 -*-
"""Hydrogeological risk index computation."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import numpy.
import numpy as np  # import numpy as np

# Import local D8 utility to build downstream indices.
from .d8 import build_downstream_index  # import .d8 import build_downstream_index


def compute_flow_accum_area_m2(d8: np.ndarray, encoding: str, cell_area_m2: float | np.ndarray, active_mask: np.ndarray) -> np.ndarray:  # define function compute_flow_accum_area_m2
    """Compute upstream contributing area (m^2) via topological traversal."""  # execute statement
    valid, ds_r, ds_c = build_downstream_index(d8, encoding)  # set valid, ds_r, ds_c
    H, W = d8.shape  # set H, W

    acc_dtype = np.float32 if (not np.isscalar(cell_area_m2) and getattr(cell_area_m2, "dtype", None) == np.float32) else np.float64
    if np.isscalar(cell_area_m2):  # check condition np.isscalar(cell_area_m2):
        acc = np.where(active_mask, float(cell_area_m2), 0.0).astype(acc_dtype)  # set acc
    else:  # fallback branch
        acc = np.where(active_mask, cell_area_m2, 0.0).astype(acc_dtype)  # set acc

    indeg = np.zeros((H, W), dtype=np.int32)  # set indeg

    rr, cc = np.nonzero(active_mask & valid)  # set rr, cc
    rds = ds_r[rr, cc]  # set rds
    cds = ds_c[rr, cc]  # set cds
    np.add.at(indeg, (rds, cds), 1)  # execute statement

    q_r, q_c = np.nonzero(active_mask & (indeg == 0))  # set q_r, q_c
    stack_r = q_r.astype(np.int32).tolist()  # set stack_r
    stack_c = q_c.astype(np.int32).tolist()  # set stack_c

    while stack_r:  # loop while stack_r:
        r = stack_r.pop()  # set r
        c = stack_c.pop()  # set c
        if not (active_mask[r, c] and valid[r, c]):  # check condition not (active_mask[r, c] and valid[r, c]):
            continue  # continue loop
        rd = int(ds_r[r, c])  # set rd
        cd = int(ds_c[r, c])  # set cd
        if rd < 0 or cd < 0:  # check condition rd < 0 or cd < 0:
            continue  # continue loop
        acc[rd, cd] += acc[r, c]  # execute statement
        indeg[rd, cd] -= 1  # execute statement
        if indeg[rd, cd] == 0 and active_mask[rd, cd]:  # check condition indeg[rd, cd] == 0 and active_mask[rd, cd]:
            stack_r.append(rd)  # execute statement
            stack_c.append(cd)  # execute statement

    return acc  # return acc


def robust_normalize(a: np.ndarray, mask: np.ndarray, p_low: float, p_high: float) -> np.ndarray:  # define function robust_normalize
    """Robust normalization to [0,1] using percentiles."""  # execute statement
    target_dtype = np.float32 if getattr(a, "dtype", None) == np.float32 else np.float64
    x = np.asarray(a, dtype=target_dtype)  # set x
    valid = mask & np.isfinite(x)  # set valid
    if not np.any(valid):  # check condition not np.any(valid):
        return np.where(mask, target_dtype(0.0), np.nan).astype(target_dtype)  # return np.where(mask, 0.0, np.nan).astype(target_dtype)

    x_masked = np.where(valid, x, np.nan)  # set x_masked
    lo = target_dtype(np.nanpercentile(x_masked, p_low))  # set lo
    hi = target_dtype(np.nanpercentile(x_masked, p_high))  # set hi

    if (not np.isfinite(lo)) or (not np.isfinite(hi)) or hi <= lo:  # check condition (not np.isfinite(lo)) or (not np.isfinite(hi)) or hi <= lo:
        positives = x_masked[x_masked > target_dtype(0.0)]  # set positives
        if positives.size == 0:  # check condition positives.size == 0:
            return np.where(mask, target_dtype(0.0), np.nan).astype(target_dtype)  # return np.where(mask, 0.0, np.nan).astype(target_dtype)
        lo = target_dtype(np.nanpercentile(positives, p_low))  # set lo
        hi = target_dtype(np.nanpercentile(positives, p_high))  # set hi

    if (not np.isfinite(lo)) or (not np.isfinite(hi)) or hi <= lo:  # check condition (not np.isfinite(lo)) or (not np.isfinite(hi)) or hi <= lo:
        return np.where(mask, target_dtype(0.0), np.nan).astype(target_dtype)  # return np.where(mask, 0.0, np.nan).astype(target_dtype)

    y = (x - lo) / (hi - lo)  # set y
    return np.where(valid, np.clip(y, target_dtype(0.0), target_dtype(1.0)), np.nan).astype(target_dtype)  # return np.where(valid, np.clip(y, 0.0, 1.0), np.nan).astype(target_dtype)


def compute_risk_index(runoff_cum_mm: np.ndarray, flow_accum_m2: np.ndarray, active_mask: np.ndarray,  # define function compute_risk_index
                       balance: float, p_low: float, p_high: float) -> np.ndarray:  # execute statement
    """Combine normalized runoff and flow accumulation into a unitless risk index."""  # execute statement
    alpha = float(np.clip(balance, 0.0, 1.0))  # set alpha
    r1 = robust_normalize(runoff_cum_mm, active_mask, p_low, p_high)  # set r1
    r2 = robust_normalize(flow_accum_m2, active_mask, p_low, p_high)  # set r2
    blended = alpha * r1 + (1.0 - alpha) * (r1 * r2)  # set blended
    return np.where(active_mask, blended, np.nan).astype(np.float32)  # return np.where(active_mask, blended, np.nan).astype(np.float32)
