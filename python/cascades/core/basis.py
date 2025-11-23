"""
Vibrational basis state generation for spectroscopy simulations.

Functions for creating quantum basis sets and computing vibrational energies.
"""

import numpy as np
from numpy.typing import NDArray
from math import comb


def rcbasis(nquanta: int, nmodes: int, target: int) -> NDArray[np.int64]:
    """
    Recursively generate integer partitions with sum constraint.

    Creates all combinations of non-negative integers that sum to target,
    distributed across nmodes columns, with each value <= nquanta.

    Parameters
    ----------
    nquanta : int
        Maximum value for any single element.
    nmodes : int
        Number of integers in each partition.
    target : int
        Required sum of integers in each row.

    Returns
    -------
    NDArray[np.int64]
        Matrix where each row is a valid partition, shape (n_partitions, nmodes).

    Examples
    --------
    >>> rcbasis(2, 2, 2)
    array([[0, 2],
           [1, 1],
           [2, 0]])
    """
    if nmodes == 1:
        if target <= nquanta:
            return np.array([[target]], dtype=np.int64)
        return np.empty((0, 1), dtype=np.int64)

    rows = []
    for ii in range(min(nquanta, target) + 1):
        sub = rcbasis(nquanta, nmodes - 1, target - ii)
        if sub.size > 0:
            prefix = np.full((sub.shape[0], 1), ii, dtype=np.int64)
            rows.append(np.hstack([prefix, sub]))

    if not rows:
        return np.empty((0, nmodes), dtype=np.int64)
    return np.vstack(rows)


def basis_tc(
    nmode: int,
    nquanta: int,
    wvib: NDArray[np.float64]
) -> tuple[NDArray[np.int64], NDArray[np.float64]]:
    """
    Create vibrational basis set and compute state energies.

    Generates all vibrational basis states up to a given number of total quanta
    and calculates the energy of each state.

    Parameters
    ----------
    nmode : int
        Number of vibrational modes to include in basis set.
    nquanta : int
        Maximum total quanta across all modes.
    wvib : NDArray[np.float64]
        Vibrational frequencies in cm^-1, length >= nmode.

    Returns
    -------
    base : NDArray[np.int64]
        Matrix of quantum numbers, shape (nstates, nmode).
        Each row represents a basis state.
    wviball : NDArray[np.float64]
        Vibrational energy for each basis state in cm^-1.

    Raises
    ------
    ValueError
        If nmode > len(wvib) or wvib is empty.

    Examples
    --------
    >>> base, energies = basis_tc(2, 2, np.array([100.0, 200.0]))
    >>> base
    array([[0, 0],
           [0, 1],
           [1, 0],
           [0, 2],
           [1, 1],
           [2, 0]])
    """
    assert nmode <= len(wvib), (
        f"nmode ({nmode}) exceeds available modes in wvib ({len(wvib)})"
    )
    assert len(wvib) > 0, "wvib cannot be empty"

    # Generate basis states
    rows = [rcbasis(nq, nmode, nq) for nq in range(nquanta + 1)]
    base = np.vstack(rows)

    nstates = comb(nquanta + nmode, nmode)
    assert base.shape[0] == nstates

    # Compute energies: E = sum((n + 0.5) * w) for each mode
    wviball = np.sum((base + 0.5) * wvib[:nmode], axis=1)

    return base, wviball
