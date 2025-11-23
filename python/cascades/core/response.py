"""
Response function calculations for cascade spectroscopy.

Compute third-order cascade and fifth-order direct response functions
for 2D resonance Raman spectroscopy simulations.

Reference: T. P. Cheshire and A. M. Moran, J. Chem. Phys. 151, 104203 (2019)
https://doi.org/10.1063/1.5115401

Pathway Notation
----------------
Each response function term corresponds to a specific Liouville pathway through
the density matrix. The pathways are characterized by:

- Lp: Absorption-side electronic transition (ket interaction, +k)
- Lm: Emission-side electronic transition (bra interaction, -k)
- Lc: Vibrational coherence propagation

Third-order cascade pathways (f[0-3], fc[0-3]):
    f[0], fc[0]: |g,m><g,m| → |e,n><g,m| → |g,k><g,m| → |e,l><g,m|
    f[1], fc[1]: |g,m><g,m| → |g,m><e,n| → |g,m><g,k| → |e,l><g,k|
    f[2], fc[2]: |g,m><g,m| → |e,n><g,m| → |e,n><g,k| → |e,n><g,l|
    f[3], fc[3]: |g,m><g,m| → |g,m><e,n| → |g,k><e,n| → |e,k><g,l|

Fifth-order direct pathways (r[0-15]):
    r[0-3]:  Type I   - Both vibrational coherences on ket side
    r[4-7]:  Type II  - Both vibrational coherences on bra side
    r[8-11]: Type III - First coherence on ket, second on bra
    r[12-15]: Type IV - First coherence on bra, second on ket
"""

import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
from scipy.signal import convolve


@dataclass
class LaserParameters:
    """Laser parameters for spectroscopy simulation."""
    dt: float  # Time step in fs
    nt: int    # Number of time points
    w_L: float # Laser frequency in cm^-1


@dataclass
class MaterialParameters:
    """Material parameters for spectroscopy simulation."""
    gamma_eg: float           # Electronic dephasing rate in cm^-1
    gamma_vib: float          # Vibrational dephasing rate in cm^-1
    weg: float                # Electronic energy gap in cm^-1
    wvib: NDArray[np.float64] # Vibrational frequencies in cm^-1


def cascade_2drr_res(
    e_vib: NDArray[np.float64],
    nquanta: int,
    ovlp: NDArray[np.float64],
    laser_params: LaserParameters,
    material_params: MaterialParameters,
    return_spectrum: bool = False
) -> tuple[float, float, float] | tuple[float, float, float, NDArray[np.float64], NDArray[np.float64]]:
    """
    Compute third-order cascade and fifth-order direct response functions.

    Calculates the cascade-to-direct signal ratio for 2D resonance Raman
    spectroscopy. Includes both sequential and parallel cascade contributions.

    Parameters
    ----------
    e_vib : NDArray[np.float64]
        Vibrational energies for each basis state in cm^-1.
    nquanta : int
        Number of vibrational quanta (determines basis size).
    ovlp : NDArray[np.float64]
        Overlap integrals between basis states, shape (nstates, nstates).
    laser_params : LaserParameters
        Laser parameters (dt, nt, w_L).
    material_params : MaterialParameters
        Material parameters (gamma_eg, gamma_vib, weg, wvib).
    return_spectrum : bool, optional
        If True, also return frequency grid and full parallel cascade spectrum.
        Default is False.

    Returns
    -------
    ratio : float
        Cascade-to-direct signal ratio (at first vibrational frequency).
    cas : float
        Absolute cascade signal (at first vibrational frequency).
    dir_ : float
        Absolute direct signal.
    ww : NDArray[np.float64], optional
        Frequency grid (cm^-1). Only returned if return_spectrum=True.
    par_spectrum : NDArray[np.float64], optional
        Full parallel cascade spectrum vs frequency. Only returned if
        return_spectrum=True. Use with ww to extract signals at other
        vibrational frequencies.

    Notes
    -----
    kT is fixed at 200 cm^-1 for Boltzmann population calculations.
    Speed of light: c = 2.998e-5 cm/fs.

    The primary output uses the first vibrational mode (wvib[0]). To extract
    cascade signals at other vibrational frequencies, use return_spectrum=True
    and find the appropriate index in ww.
    """
    # Boltzmann populations
    kT = 200.0  # cm^-1
    boltz_factor = np.exp(-e_vib / kT)
    partition_func = np.sum(boltz_factor)
    boltz_pop = boltz_factor / partition_func

    # Unpack parameters
    dt = laser_params.dt
    nt = laser_params.nt
    w_L = laser_params.w_L

    gamma_eg = material_params.gamma_eg
    gamma_vib = material_params.gamma_vib
    weg = material_params.weg
    wvib = material_params.wvib

    iq = nquanta
    c = 2.998e-5  # Speed of light in cm/fs

    # Vibrational energy gaps
    wi, wf = np.meshgrid(e_vib, e_vib)
    w = wf - wi

    # Initialize response function terms
    r = [np.zeros(iq, dtype=np.complex128) for _ in range(16)]
    f = [np.zeros(iq, dtype=np.complex128) for _ in range(4)]
    fc = [np.zeros(iq, dtype=np.complex128) for _ in range(4)]

    # Frequency grid
    dw = 1.0 / nt
    ww = np.arange(-0.5, 0.5, dw) / dt / c
    nw = len(ww)

    fpa = [np.zeros((nw, iq), dtype=np.complex128) for _ in range(4)]
    fp = [np.zeros((nw, iq), dtype=np.complex128) for _ in range(4)]

    # 3D grids for parallel cascades
    W, Wf, Wi = np.meshgrid(ww, e_vib, e_vib, indexing='ij')
    W = -W + Wf - Wi

    # Lineshape functions
    Lp = 1.0 / (w_L - weg - w + 1j * gamma_eg)
    Lm = 1.0 / (-w_L + weg - w + 1j * gamma_eg)
    Lc0 = 1.0 / (-wvib[0] - w + 1j * gamma_vib)
    Lc0C = 1.0 / (wvib[0] + w + 1j * gamma_vib)
    Lc1 = Lc0.copy()
    Lc2 = Lc0.copy()
    LcP = 1.0 / (-W + 1j * gamma_vib)

    # Zero out diagonal terms
    diag_mask = (w == 0)
    Lc0[diag_mask] = 0
    Lc0C[diag_mask] = 0
    Lc1[diag_mask] = 0
    Lc2[diag_mask] = 0
    # For 3D array LcP, broadcast the 2D mask
    LcP[:, diag_mask] = 0

    # Main computation loop - assuming m=0 (ground state)
    m = 0
    bp_m = boltz_pop[m]

    for n in range(iq):
        for k in range(iq):
            for l in range(iq):
                # Third-order auxiliary functions
                f[0][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * Lc0[k, m] * Lp[l, m]
                f[1][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * Lc0[m, k] * Lp[l, k]
                fc[0][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * Lc0C[k, m] * Lp[l, m]
                fc[1][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * Lc0C[m, k] * Lp[l, k]

                f[2][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * Lc0[n, k] * Lp[n, l]
                f[3][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * Lc0[k, n] * Lp[k, l]
                fc[2][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * Lc0C[n, k] * Lp[n, l]
                fc[3][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * Lc0C[k, n] * Lp[k, l]

                # Parallel cascade terms
                fpa[0][:, n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * LcP[:, k, m] * Lp[l, m]
                fpa[1][:, n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * LcP[:, m, k] * Lp[l, k]
                fpa[2][:, n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * LcP[:, n, k] * Lp[n, l]
                fpa[3][:, n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * LcP[:, k, n] * Lp[k, l]

                fp[0][:, n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * LcP[:, k, m] * Lc0[k, m] * Lp[l, m]
                fp[1][:, n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * LcP[:, m, k] * Lc0[m, k] * Lp[l, k]
                fp[2][:, n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * LcP[:, n, k] * Lc0[n, k] * Lp[n, l]
                fp[3][:, n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * LcP[:, k, n] * Lc0[k, n] * Lp[k, l]

                # Fifth-order direct terms (16 Liouville pathways)
                # See module docstring for pathway classification
                for u in range(iq):
                    for v in range(iq):
                        # Type I: Both vibrational coherences evolve on ket side
                        # r[0]: Lp-Lc1-Lp-Lc2-Lp (all ket-side)
                        r[0][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, u] * ovlp[v, u] * ovlp[v, m] * Lp[n, m] * Lc1[k, m] * Lp[l, m] * Lc2[u, m] * Lp[v, m]
                        # r[1]: Lp-Lc1-Lm-Lc2-Lp (ket-bra-ket)
                        r[1][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, u] * ovlp[v, k] * ovlp[v, u] * Lp[n, m] * Lc1[k, m] * Lm[k, l] * Lc2[k, u] * Lp[v, u]
                        # r[2]: Lm-Lc1-Lm-Lc2-Lp (all bra-side except final)
                        r[2][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, u] * ovlp[v, m] * ovlp[v, u] * Lm[m, n] * Lc1[m, k] * Lm[m, l] * Lc2[m, u] * Lp[v, u]
                        # r[3]: Lm-Lc1-Lp-Lc2-Lp (bra-ket-ket)
                        r[3][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, u] * ovlp[v, u] * ovlp[v, k] * Lm[m, n] * Lc1[m, k] * Lp[l, k] * Lc2[u, k] * Lp[v, k]

                        # Type II: Both vibrational coherences evolve on bra side
                        # r[4]: Lp-Lc1-Lm-Lc2-Lp
                        r[4][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[u, l] * ovlp[k, v] * ovlp[u, v] * Lp[n, m] * Lc1[n, k] * Lm[l, k] * Lc2[u, k] * Lp[u, v]
                        # r[5]: Lp-Lc1-Lp-Lc2-Lp
                        r[5][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[u, l] * ovlp[u, v] * ovlp[n, v] * Lp[n, m] * Lc1[n, k] * Lp[n, l] * Lc2[n, u] * Lp[n, v]
                        # r[6]: Lm-Lc1-Lm-Lc2-Lp
                        r[6][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[u, l] * ovlp[n, v] * ovlp[u, v] * Lm[m, n] * Lc1[k, n] * Lm[l, n] * Lc2[u, n] * Lp[u, v]
                        # r[7]: Lm-Lc1-Lp-Lc2-Lp
                        r[7][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[u, l] * ovlp[u, v] * ovlp[k, v] * Lm[m, n] * Lc1[k, n] * Lp[k, l] * Lc2[k, u] * Lp[k, v]

                        # Type III: First coherence ket-side, second bra-side
                        # r[8]: Lp-Lc1-Lp-Lc2-Lp
                        r[8][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[u, m] * ovlp[u, v] * ovlp[l, v] * Lp[n, m] * Lc1[k, m] * Lp[l, m] * Lc2[l, u] * Lp[l, v]
                        # r[9]: Lp-Lc1-Lm-Lc2-Lp
                        r[9][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[u, k] * ovlp[l, v] * ovlp[u, v] * Lp[n, m] * Lc1[k, m] * Lm[k, l] * Lc2[u, l] * Lp[u, v]
                        # r[10]: Lm-Lc1-Lp-Lc2-Lp
                        r[10][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[u, k] * ovlp[u, v] * ovlp[l, v] * Lm[m, n] * Lc1[m, k] * Lp[l, k] * Lc2[l, u] * Lp[l, v]
                        # r[11]: Lm-Lc1-Lm-Lc2-Lp
                        r[11][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[u, m] * ovlp[l, v] * ovlp[u, v] * Lm[m, n] * Lc1[m, k] * Lm[m, l] * Lc2[u, l] * Lp[u, v]

                        # Type IV: First coherence bra-side, second ket-side
                        # r[12]: Lp-Lc1-Lm-Lc2-Lp
                        r[12][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, u] * ovlp[v, l] * ovlp[v, u] * Lp[n, m] * Lc1[n, k] * Lm[l, k] * Lc2[l, u] * Lp[v, u]
                        # r[13]: Lp-Lc1-Lp-Lc2-Lp
                        r[13][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, u] * ovlp[v, u] * ovlp[v, l] * Lp[n, m] * Lc1[n, k] * Lp[n, l] * Lc2[u, l] * Lp[v, l]
                        # r[14]: Lm-Lc1-Lm-Lc2-Lp
                        r[14][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, u] * ovlp[v, l] * ovlp[v, u] * Lm[m, n] * Lc1[k, n] * Lm[l, n] * Lc2[l, u] * Lp[v, u]
                        # r[15]: Lm-Lc1-Lp-Lc2-Lp
                        r[15][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, u] * ovlp[v, u] * ovlp[v, l] * Lm[m, n] * Lc1[k, n] * Lp[k, l] * Lc2[u, l] * Lp[v, l]

    # Sum over states
    f_sum = [np.sum(fi) for fi in f]
    fc_sum = [np.sum(fci) for fci in fc]
    fpa_sum = [np.sum(fpai, axis=1) for fpai in fpa]
    fp_sum = [np.sum(fpi, axis=1) for fpi in fp]
    r_sum = [np.sum(ri) for ri in r]

    # Direct fifth-order signal
    direct = sum(r_sum)

    # Sequential cascades
    seq1 = sum(fi**2 for fi in f_sum)
    seq1 += 2 * (f_sum[0] * f_sum[1] + f_sum[0] * f_sum[2] + f_sum[0] * f_sum[3] +
                 f_sum[1] * f_sum[2] + f_sum[1] * f_sum[3] + f_sum[2] * f_sum[3])

    seq2 = sum(fci * fi for fci, fi in zip(fc_sum, f_sum))
    seq2 += (fc_sum[0] * f_sum[1] + fc_sum[0] * f_sum[2] + fc_sum[0] * f_sum[3] +
             fc_sum[1] * f_sum[0] + fc_sum[1] * f_sum[2] + fc_sum[1] * f_sum[3] +
             fc_sum[2] * f_sum[0] + fc_sum[2] * f_sum[1] + fc_sum[2] * f_sum[3] +
             fc_sum[3] * f_sum[0] + fc_sum[3] * f_sum[1] + fc_sum[3] * f_sum[2])

    seq = seq1 + seq2

    # Parallel cascades
    wr = 1.0 / nw
    fpa_total = sum(fpa_sum)
    fp_total = sum(fp_sum)
    par1 = convolve(fpa_total, fp_total, mode='same') * wr
    par = 2 * par1

    # Find index closest to vibrational frequency
    iomega = np.argmin(np.abs(ww - wvib[0]))
    cascade = seq + par[iomega]

    cas = np.abs(cascade)
    dir_ = np.abs(direct)
    ratio = cas / dir_

    if return_spectrum:
        # Return full parallel cascade spectrum for multi-frequency extraction
        par_spectrum = np.abs(seq + par)
        return float(ratio), float(cas), float(dir_), ww, par_spectrum

    return float(ratio), float(cas), float(dir_)
