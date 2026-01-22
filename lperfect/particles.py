# -*- coding: utf-8 -*-
"""Particle container and helpers."""  # execute statement

# NOTE: Rain NetCDF inputs follow cdl/rain_time_dependent.cdl (CF-1.10).

# Import dataclass for a simple structured object.
from dataclasses import dataclass  # import dataclasses import dataclass

# Import numpy for arrays.
import numpy as np  # import numpy as np


@dataclass  # apply decorator
class Particles:  # define class Particles
    """Structure-of-arrays particle container.

    Attributes
    ----------
    r, c : np.ndarray (int32)
        Row and column indices of each particle.
    vol : np.ndarray (float64)
        Particle volume in cubic meters.
    tau : np.ndarray (float64)
        Time-to-next-hop in seconds (<=0 means eligible to move).
    """
    r: np.ndarray  # execute statement
    c: np.ndarray  # execute statement
    vol: np.ndarray  # execute statement
    tau: np.ndarray  # execute statement


@dataclass  # apply decorator
class ParticleBuffer:  # define class ParticleBuffer
    """Growable structure-of-arrays buffer for particle fields."""  # execute statement

    r: np.ndarray  # execute statement
    c: np.ndarray  # execute statement
    vol: np.ndarray  # execute statement
    tau: np.ndarray  # execute statement
    size: int  # execute statement

    @classmethod
    def from_particles(cls, p: Particles) -> "ParticleBuffer":
        """Create a buffer initialized with existing particles."""  # execute statement
        n = int(p.r.size)  # set n
        return cls(  # return cls(
            r=p.r.copy(),  # set r
            c=p.c.copy(),  # set c
            vol=p.vol.copy(),  # set vol
            tau=p.tau.copy(),  # set tau
            size=n,  # set size
        )  # execute statement

    @property
    def capacity(self) -> int:
        """Return current allocated capacity."""  # execute statement
        return int(self.r.size)  # return int(self.r.size)

    def _ensure_capacity(self, extra: int) -> None:
        """Ensure arrays can fit `extra` additional particles."""  # execute statement
        needed = self.size + int(extra)  # set needed
        if needed <= self.capacity:  # check condition needed <= self.capacity:
            return  # return
        new_cap = max(needed, max(1, self.capacity) * 2)  # set new_cap
        self.r = _grow_array(self.r, new_cap, dtype=self.r.dtype)  # execute statement
        self.c = _grow_array(self.c, new_cap, dtype=self.c.dtype)  # execute statement
        self.vol = _grow_array(self.vol, new_cap, dtype=self.vol.dtype)  # execute statement
        self.tau = _grow_array(self.tau, new_cap, dtype=self.tau.dtype)  # execute statement

    def append_arrays(self, r: np.ndarray, c: np.ndarray, vol: np.ndarray, tau: np.ndarray) -> None:
        """Append particle arrays to the buffer."""  # execute statement
        n_new = int(r.size)  # set n_new
        if n_new == 0:  # check condition n_new == 0:
            return  # return
        self._ensure_capacity(n_new)  # execute statement
        end = self.size + n_new  # set end
        self.r[self.size : end] = r  # execute statement
        self.c[self.size : end] = c  # execute statement
        self.vol[self.size : end] = vol  # execute statement
        self.tau[self.size : end] = tau  # execute statement
        self.size = end  # set size

    def to_particles(self) -> Particles:
        """Return a trimmed Particles view of the buffer."""  # execute statement
        return Particles(  # return Particles(
            r=self.r[: self.size].copy(),  # set r
            c=self.c[: self.size].copy(),  # set c
            vol=self.vol[: self.size].copy(),  # set vol
            tau=self.tau[: self.size].copy(),  # set tau
        )  # execute statement


def _grow_array(arr: np.ndarray, new_size: int, dtype: np.dtype) -> np.ndarray:
    """Grow a 1D numpy array to a new size."""  # execute statement
    expanded = np.empty(int(new_size), dtype=dtype)  # set expanded
    expanded[: arr.size] = arr  # execute statement
    return expanded  # return expanded


def empty_particles() -> Particles:  # define function empty_particles
    """Create an empty particle container."""  # execute statement
    return Particles(  # return Particles(
        r=np.zeros(0, dtype=np.int32),  # set r
        c=np.zeros(0, dtype=np.int32),  # set c
        vol=np.zeros(0, dtype=np.float64),  # set vol
        tau=np.zeros(0, dtype=np.float64),  # set tau
    )  # execute statement


def concat_particles(a: Particles, b: Particles) -> Particles:  # define function concat_particles
    """Concatenate two particle containers."""  # execute statement
    # If one is empty, return the other.
    if a.r.size == 0:  # check condition a.r.size == 0:
        return b  # return b
    if b.r.size == 0:  # check condition b.r.size == 0:
        return a  # return a
    # Concatenate fields.
    return Particles(  # return Particles(
        r=np.concatenate([a.r, b.r]),  # set r
        c=np.concatenate([a.c, b.c]),  # set c
        vol=np.concatenate([a.vol, b.vol]),  # set vol
        tau=np.concatenate([a.tau, b.tau]),  # set tau
    )  # execute statement


def pack_particles_to_float64(p: Particles) -> np.ndarray:  # define function pack_particles_to_float64
    """Pack particles into float64 matrix (N,4) for MPI transfer."""  # execute statement
    # Allocate (N,4).
    buf = np.empty((p.r.size, 4), dtype=np.float64)  # set buf
    # Copy fields to columns (cast ints to float for a single MPI datatype).
    buf[:, 0] = p.r.astype(np.float64)  # execute statement
    buf[:, 1] = p.c.astype(np.float64)  # execute statement
    buf[:, 2] = p.vol.astype(np.float64)  # execute statement
    buf[:, 3] = p.tau.astype(np.float64)  # execute statement
    # Return packed matrix.
    return buf  # return buf


def unpack_particles_from_float64(buf: np.ndarray) -> Particles:  # define function unpack_particles_from_float64
    """Unpack float64 matrix (N,4) back to Particles."""  # execute statement
    # Handle empty case.
    if buf.size == 0:  # check condition buf.size == 0:
        return empty_particles()  # return empty_particles()
    # Build Particles with proper dtypes.
    return Particles(  # return Particles(
        r=buf[:, 0].astype(np.int32),  # set r
        c=buf[:, 1].astype(np.int32),  # set c
        vol=buf[:, 2].astype(np.float64),  # set vol
        tau=buf[:, 3].astype(np.float64),  # set tau
    )  # execute statement
