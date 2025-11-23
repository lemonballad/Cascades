"""
Numba-optimized versions of core computation functions.

Install with: pip install cascades[fast]

These functions provide significant speedup for larger basis sets.
Falls back to pure numpy if numba is not available.
"""

import numpy as np
from numpy.typing import NDArray

try:
    from numba import jit, prange
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False
    # Create no-op decorator
    def jit(*args, **kwargs):
        def decorator(func):
            return func
        return decorator
    prange = range


@jit(nopython=True, cache=True)
def _fcfac2_tc_fast(d: float, nquanta: int) -> NDArray[np.float64]:
    """
    Numba-optimized Franck-Condon overlap integral calculation.

    10-50x faster than pure Python for nquanta > 5.
    """
    s = d**2 / 2
    o_0_0 = np.exp(-s / 2)

    size = nquanta + 1
    ovlp = np.zeros((size, size), dtype=np.float64)

    # Precompute factorials
    fact = np.zeros(size, dtype=np.float64)
    fact[0] = 1.0
    for i in range(1, size):
        fact[i] = fact[i-1] * i

    for mm in range(size):
        for nn in range(size):
            m = mm
            n = nn

            num = 0.0
            k_max = min(mm, nn) + 1
            for k in range(k_max):
                mk = m - k
                nk = n - k

                o_mk_g = o_0_0 * ((-1)**mk) * (s**(mk / 2)) / np.sqrt(fact[mk])
                o_g_nk = o_0_0 * (s**(nk / 2)) / np.sqrt(fact[nk])

                num += (1.0 / fact[k]
                       * np.sqrt(fact[m] * fact[n] / fact[mk] / fact[nk])
                       * o_mk_g * o_g_nk)

            ovlp[mm, nn] = num / o_0_0

    return ovlp


@jit(nopython=True, cache=True, parallel=True)
def _fcinfo_tc_fast(
    base: NDArray[np.int64],
    disp: NDArray[np.float64],
    nmode: int,
    nquanta: int
) -> NDArray[np.float64]:
    """
    Numba-optimized product of overlap integrals.

    Uses parallel loops for additional speedup on multi-core systems.
    """
    nstates = base.shape[0]

    # Compute overlap integrals for each mode
    ovlpall = np.zeros((nmode, nquanta + 1, nquanta + 1), dtype=np.float64)
    for im in range(nmode):
        ovlpall[im] = _fcfac2_tc_fast(disp[im], nquanta)

    # Compute products of overlaps for all state pairs
    fcall = np.zeros((nstates, nstates), dtype=np.float64)

    for iq_o in prange(nstates):
        for iq_i in range(nstates):
            num = 1.0
            for iv in range(nmode):
                j1 = base[iq_o, iv]
                j2 = base[iq_i, iv]
                num *= ovlpall[iv, j1, j2]
            fcall[iq_o, iq_i] = num

    return fcall


def fcfac2_tc_fast(d: float, nquanta: int) -> NDArray[np.float64]:
    """
    Fast Franck-Condon overlap integral calculation.

    Uses Numba JIT compilation if available, otherwise falls back
    to the pure numpy implementation.

    Parameters
    ----------
    d : float
        Unitless mode displacement.
    nquanta : int
        Maximum number of vibrational quanta.

    Returns
    -------
    NDArray[np.float64]
        Overlap integral matrix, shape (nquanta+1, nquanta+1).
    """
    if HAS_NUMBA:
        return _fcfac2_tc_fast(d, nquanta)
    else:
        from cascades.core.franck_condon import fcfac2_tc
        return fcfac2_tc(d, nquanta)


def fcinfo_tc_fast(
    base: NDArray[np.int64],
    disp: NDArray[np.float64],
    nmode: int,
    nquanta: int
) -> NDArray[np.float64]:
    """
    Fast product of overlap integrals calculation.

    Uses Numba JIT compilation with parallel loops if available.

    Parameters
    ----------
    base : NDArray[np.int64]
        Matrix of quantum numbers for each basis state.
    disp : NDArray[np.float64]
        Unitless displacements for each vibrational mode.
    nmode : int
        Number of vibrational modes.
    nquanta : int
        Maximum number of quanta.

    Returns
    -------
    NDArray[np.float64]
        Product of overlaps for each state pair.
    """
    if HAS_NUMBA:
        return _fcinfo_tc_fast(base, disp, nmode, nquanta)
    else:
        from cascades.core.franck_condon import fcinfo_tc
        return fcinfo_tc(base, disp, nmode, nquanta)


def check_numba_available() -> bool:
    """Check if Numba is available for acceleration."""
    return HAS_NUMBA
