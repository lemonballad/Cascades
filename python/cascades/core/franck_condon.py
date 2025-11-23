"""
Franck-Condon overlap integral calculations.

Compute vibrational overlap integrals using Heller's method for
spectroscopy simulations.
"""

import numpy as np
from numpy.typing import NDArray
from math import factorial, sqrt


def fcfac2_tc(
    d: float,
    nquanta: int,
    wn: float = 1.0,
    wm: float = 1.0
) -> NDArray[np.float64]:
    """
    Compute Franck-Condon overlap integrals by Heller's method.

    Calculates <m|n> overlap integrals between vibrational states using
    the method from Myers et al. JCP 77, 3857 (1982), Equation 7.

    Parameters
    ----------
    d : float
        Unitless mode displacement.
    nquanta : int
        Maximum number of vibrational quanta.
    wn : float, optional
        Ground state vibrational frequency (default 1.0).
    wm : float, optional
        Excited state vibrational frequency (default 1.0).

    Returns
    -------
    NDArray[np.float64]
        Overlap integral matrix, shape (nquanta+1, nquanta+1).
        First index is excited state, second is ground state.

    Notes
    -----
    The Huang-Rhys factor is S = d^2/2.
    """
    # Huang-Rhys factor
    s = d**2 / 2
    # <0|0> overlap
    o_0_0 = np.exp(-s / 2)

    size = nquanta + 1
    ovlp = np.zeros((size, size), dtype=np.float64)

    for mm in range(size):
        for nn in range(size):
            m = mm  # excited state quantum number
            n = nn  # ground state quantum number

            num = 0.0
            for k in range(min(mm, nn) + 1):
                o_mk_g = o_0_0 * (-1)**(m - k) * s**((m - k) / 2) / sqrt(factorial(m - k))
                o_g_nk = o_0_0 * s**((n - k) / 2) / sqrt(factorial(n - k))

                num += (
                    1.0 / factorial(k)
                    * sqrt(factorial(m) * factorial(n) / factorial(m - k) / factorial(n - k))
                    * o_mk_g * o_g_nk
                )

            # First index (mm) is excited state, second (nn) is ground state
            ovlp[mm, nn] = num / o_0_0

    return ovlp


def fcinfo_tc(
    base: NDArray[np.int64],
    disp: NDArray[np.float64],
    nmode: int,
    nquanta: int
) -> NDArray[np.float64]:
    """
    Calculate products of overlap integrals for all basis state pairs.

    For each pair of basis states, computes the product of single-mode
    overlap integrals across all vibrational modes.

    Parameters
    ----------
    base : NDArray[np.int64]
        Matrix of quantum numbers for each basis state, shape (nstates, nmode).
    disp : NDArray[np.float64]
        Unitless displacements for each vibrational mode.
    nmode : int
        Number of vibrational modes.
    nquanta : int
        Maximum number of quanta (determines overlap matrix size).

    Returns
    -------
    NDArray[np.float64]
        Product of overlaps for each state pair, shape (nstates, nstates).

    Notes
    -----
    Currently assumes ground and excited state frequencies are equal (wn = wm = 1).
    """
    nstates = base.shape[0]

    # Compute overlap integrals for each mode
    ovlpall_tc = np.zeros((nmode, nquanta + 1, nquanta + 1), dtype=np.float64)
    for im in range(nmode):
        ovlpall_tc[im] = fcfac2_tc(disp[im], nquanta)

    # Compute products of overlaps for all state pairs
    fcall = np.zeros((nstates, nstates), dtype=np.float64)

    for iq_o in range(nstates):
        for iq_i in range(nstates):
            num = 1.0
            for iv in range(nmode):
                j1 = base[iq_o, iv]
                j2 = base[iq_i, iv]
                num *= ovlpall_tc[iv, j1, j2]
            fcall[iq_o, iq_i] = num

    return fcall
