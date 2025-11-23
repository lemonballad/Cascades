"""
Off-resonance response functions for 2DRR spectroscopy.

Compute cascade and direct signals for off-resonance conditions
including solvent contributions.

Reference: T. P. Cheshire and A. M. Moran, J. Chem. Phys. 151, 104203 (2019)
https://doi.org/10.1063/1.5115401

This module implements Section IV.C of the paper, which discusses cascades
involving both solute and solvent molecules. The cascade signal includes:
1. Solute-solute sequential cascades (dominant)
2. Solute-solvent cross-cascades (first molecule solute, second solvent)
3. Solvent-solute cross-cascades (first molecule solvent, second solute)
4. Parallel cascades for all combinations
"""

import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
from scipy.signal import convolve

from cascades.core.response import LaserParameters, MaterialParameters


@dataclass
class OffResMaterialParameters(MaterialParameters):
    """Material parameters for off-resonance simulation with solvent."""
    wsolv: float  # Solvent vibrational frequency (cm^-1)


def cascade_2drr_offres(
    e_vib: NDArray[np.float64],
    nquanta: int,
    ovlp: NDArray[np.float64],
    laser_params: LaserParameters,
    material_params: OffResMaterialParameters
) -> tuple[float, float, float]:
    """
    Compute cascade and direct response for off-resonance 2DRR.

    Includes solvent vibrational mode contributions to the cascade signal.

    Parameters
    ----------
    e_vib : NDArray[np.float64]
        Vibrational energies for each basis state (cm^-1).
    nquanta : int
        Number of vibrational quanta.
    ovlp : NDArray[np.float64]
        Overlap integrals between basis states.
    laser_params : LaserParameters
        Laser parameters.
    material_params : OffResMaterialParameters
        Material parameters including solvent frequency.

    Returns
    -------
    ratio : float
        Cascade-to-direct signal ratio.
    cascade : float
        Absolute cascade signal.
    direct : float
        Absolute direct signal.
    """
    # Unpack parameters
    dt = laser_params.dt
    nt = laser_params.nt
    w_L = laser_params.w_L

    gamma_eg = material_params.gamma_eg
    gamma_vib = material_params.gamma_vib
    weg = material_params.weg
    wvib = material_params.wvib[0] if hasattr(material_params.wvib, '__len__') else material_params.wvib
    wSolv = material_params.wsolv

    # Solvent vibrational energy levels
    wsolv = (np.arange(nquanta + 1) + 0.5) * wSolv

    # Boltzmann populations
    kT = 200.0  # cm^-1
    boltz_factor = np.exp(-e_vib / kT)
    boltz_factor_solv = np.exp(-wsolv / kT)
    partition_func = np.sum(boltz_factor)
    partition_func_solv = np.sum(boltz_factor_solv)
    boltz_pop = boltz_factor / partition_func
    boltz_pop_solv = boltz_factor_solv / partition_func_solv

    iq = nquanta
    c = 2.998e-5  # Speed of light (cm/fs)

    # Vibrational energy gaps - Solute
    wi, wf = np.meshgrid(e_vib, e_vib)
    w = wf - wi

    # Vibrational energy gaps - Solvent
    wSi, wSf = np.meshgrid(wsolv, wsolv)
    wS = wSf - wSi

    # Initialize response functions
    r = [np.zeros(iq, dtype=np.complex128) for _ in range(16)]
    f_sol = [np.zeros(iq, dtype=np.complex128) for _ in range(4)]
    f_solC = [np.zeros(iq, dtype=np.complex128) for _ in range(4)]
    f_solv = [np.zeros(iq, dtype=np.complex128) for _ in range(2)]
    f_solvC = [np.zeros(iq, dtype=np.complex128) for _ in range(2)]

    # Frequency grid
    dw = 1.0 / nt
    ww = np.arange(-0.5, 0.5, dw) / dt / c
    nw = len(ww)
    iomega = np.argmin(np.abs(ww - wvib))

    fPA = [np.zeros((nw, iq), dtype=np.complex128) for _ in range(4)]
    fPB = [np.zeros((nw, iq), dtype=np.complex128) for _ in range(4)]
    fPSA = [np.zeros((nw, iq), dtype=np.complex128) for _ in range(2)]
    fPSB = [np.zeros((nw, iq), dtype=np.complex128) for _ in range(2)]

    # 3D grids
    W, Wf_3d, Wi_3d = np.meshgrid(ww, e_vib, e_vib, indexing='ij')
    W = -W + Wf_3d - Wi_3d
    WS, WSf, WSi = np.meshgrid(ww, wsolv, wsolv, indexing='ij')
    WS = -WS + WSf - WSi

    # Lineshape functions
    Lp = 1.0 / (w_L - weg - w + 1j * gamma_eg)
    Lm = 1.0 / (-w_L + weg - w + 1j * gamma_eg)
    Lc0 = 1.0 / (-wvib - w + 1j * gamma_vib)
    Lc0C = 1.0 / (-wvib + w + 1j * gamma_vib)
    Lc1 = Lc0.copy()
    Lc2 = Lc0.copy()
    LcS = 1.0 / (-wvib - wS + 1j * gamma_vib)
    LcSC = 1.0 / (-wvib + wS + 1j * gamma_vib)
    LcP = 1.0 / (W + 1j * gamma_vib)
    LcPS = 1.0 / (WS + 1j * gamma_vib)

    # Main computation loop
    for m in range(iq):
        bp_m = boltz_pop[m]

        for n in range(iq):
            # Solvent contributions (only for adjacent states)
            if abs(n - m) == 1:
                f_solv[0][m] += boltz_pop_solv[m] * LcS[n, m]
                f_solv[1][m] += boltz_pop_solv[m] * LcS[m, n]
                f_solvC[0][m] += boltz_pop_solv[m] * LcSC[n, m]
                f_solvC[1][m] += boltz_pop_solv[m] * LcSC[m, n]
                fPSA[0][:, m] += boltz_pop_solv[m] * LcPS[:, n, m]
                fPSA[1][:, m] += boltz_pop_solv[m] * LcPS[:, m, n]
                fPSB[0][:, m] += boltz_pop_solv[m] * LcPS[:, n, m] * LcS[n, m]
                fPSB[1][:, m] += boltz_pop_solv[m] * LcPS[:, m, n] * LcS[m, n]

            for k in range(iq):
                for l in range(iq):
                    # Solute third-order terms
                    if k != m:
                        f_sol[0][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * Lc0[k, m] * Lp[l, m]
                        f_sol[1][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * Lc0[m, k] * Lp[l, k]
                        f_solC[0][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * Lc0C[k, m] * Lp[l, m]
                        f_solC[1][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * Lc0C[m, k] * Lp[l, k]
                        fPA[0][:, m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * LcP[:, k, m] * Lp[l, m]
                        fPA[1][:, m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * LcP[:, m, k] * Lp[l, k]
                        fPB[0][:, m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * Lp[n, m] * LcP[:, k, m] * Lc0[m, k] * Lp[l, m]
                        fPB[1][:, m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * Lm[m, n] * LcP[:, m, k] * Lc0C[m, k] * Lp[l, k]

                    if k != n:
                        f_sol[2][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * Lc0[n, k] * Lp[n, l]
                        f_sol[3][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * Lc0[k, n] * Lp[k, l]
                        f_solC[2][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * Lc0C[n, k] * Lp[n, l]
                        f_solC[3][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * Lc0C[k, n] * Lp[k, l]
                        fPA[2][:, m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * LcP[:, n, k] * Lp[n, l]
                        fPA[3][:, m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * LcP[:, k, n] * Lp[k, l]
                        fPB[2][:, m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * Lp[n, m] * LcP[:, n, k] * Lc0C[n, k] * Lp[n, l]
                        fPB[3][:, m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * Lm[m, n] * LcP[:, k, n] * Lc0C[k, n] * Lp[k, l]

                    # Fifth-order direct terms
                    for u in range(iq):
                        for v in range(iq):
                            if m != u and m != k:
                                r[0][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, u] * ovlp[v, u] * ovlp[v, m] * Lp[n, m] * Lc1[k, m] * Lp[l, m] * Lc2[u, m] * Lp[v, m]
                            if k != u and m != k:
                                r[1][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, u] * ovlp[v, k] * ovlp[v, u] * Lp[n, m] * Lc1[k, m] * Lm[k, l] * Lc2[k, u] * Lp[v, u]
                            if m != u and m != k:
                                r[2][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, u] * ovlp[v, m] * ovlp[v, u] * Lm[m, n] * Lc1[m, k] * Lm[m, l] * Lc2[m, u] * Lp[v, u]
                            if k != u and m != k:
                                r[3][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, u] * ovlp[v, u] * ovlp[v, k] * Lm[m, n] * Lc1[m, k] * Lp[l, k] * Lc2[u, k] * Lp[v, k]

                            if k != u and n != k:
                                r[4][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[u, l] * ovlp[k, v] * ovlp[u, v] * Lp[n, m] * Lc1[n, k] * Lm[l, k] * Lc2[u, k] * Lp[u, v]
                            if n != u and n != k:
                                r[5][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[u, l] * ovlp[u, v] * ovlp[n, v] * Lp[n, m] * Lc1[n, k] * Lp[n, l] * Lc2[n, u] * Lp[n, v]
                                r[6][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[u, l] * ovlp[n, v] * ovlp[u, v] * Lm[m, n] * Lc1[k, n] * Lm[l, n] * Lc2[u, n] * Lp[u, v]
                            if k != u and n != k:
                                r[7][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[u, l] * ovlp[u, v] * ovlp[k, v] * Lm[m, n] * Lc1[k, n] * Lp[k, l] * Lc2[k, u] * Lp[k, v]

                            if l != u and m != k:
                                r[8][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[u, m] * ovlp[u, v] * ovlp[l, v] * Lp[n, m] * Lc1[k, m] * Lp[l, m] * Lc2[l, u] * Lp[l, v]
                                r[9][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[u, k] * ovlp[l, v] * ovlp[u, v] * Lp[n, m] * Lc1[k, m] * Lm[k, l] * Lc2[u, l] * Lp[u, v]
                                r[10][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[u, k] * ovlp[u, v] * ovlp[l, v] * Lm[m, n] * Lc1[m, k] * Lp[l, k] * Lc2[l, u] * Lp[l, v]
                                r[11][m] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[u, m] * ovlp[l, v] * ovlp[u, v] * Lm[m, n] * Lc1[m, k] * Lm[m, l] * Lc2[u, l] * Lp[u, v]
                            if l != u and n != k:
                                r[12][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, u] * ovlp[v, l] * ovlp[v, u] * Lp[n, m] * Lc1[n, k] * Lm[l, k] * Lc2[l, u] * Lp[v, u]
                                r[13][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, u] * ovlp[v, u] * ovlp[v, l] * Lp[n, m] * Lc1[n, k] * Lp[n, l] * Lc2[u, l] * Lp[v, l]
                                r[14][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, u] * ovlp[v, l] * ovlp[v, u] * Lm[m, n] * Lc1[k, n] * Lm[l, n] * Lc2[l, u] * Lp[v, u]
                                r[15][m] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, u] * ovlp[v, u] * ovlp[v, l] * Lm[m, n] * Lc1[k, n] * Lp[k, l] * Lc2[u, l] * Lp[v, l]

    # Sum over states
    f_sol_sum = [np.sum(f) for f in f_sol]
    f_solC_sum = [np.sum(f) for f in f_solC]
    fPA_sum = [np.sum(f, axis=1) for f in fPA]
    fPB_sum = [np.sum(f, axis=1) for f in fPB]
    f_solv_sum = [np.sum(f) for f in f_solv]
    f_solvC_sum = [np.sum(f) for f in f_solvC]
    fPSA_sum = [np.sum(f, axis=1) for f in fPSA]
    fPSB_sum = [np.sum(f, axis=1) for f in fPSB]
    r_sum = [np.sum(ri) for ri in r]

    # Direct fifth-order signal
    direct = sum(r_sum)

    # Sequential cascades
    # 1. Solute-solute sequential (dominant contribution)
    seq_solute = sum(fi**2 for fi in f_sol_sum)
    seq_solute += 2 * (f_sol_sum[0] * f_sol_sum[1] + f_sol_sum[0] * f_sol_sum[2] +
                       f_sol_sum[0] * f_sol_sum[3] + f_sol_sum[1] * f_sol_sum[2] +
                       f_sol_sum[1] * f_sol_sum[3] + f_sol_sum[2] * f_sol_sum[3])
    # Conjugate cross-terms for solute
    seq_solute += sum(fci * fi for fci, fi in zip(f_solC_sum, f_sol_sum))
    seq_solute += (f_solC_sum[0] * f_sol_sum[1] + f_solC_sum[0] * f_sol_sum[2] +
                   f_solC_sum[0] * f_sol_sum[3] + f_solC_sum[1] * f_sol_sum[0] +
                   f_solC_sum[1] * f_sol_sum[2] + f_solC_sum[1] * f_sol_sum[3] +
                   f_solC_sum[2] * f_sol_sum[0] + f_solC_sum[2] * f_sol_sum[1] +
                   f_solC_sum[2] * f_sol_sum[3] + f_solC_sum[3] * f_sol_sum[0] +
                   f_solC_sum[3] * f_sol_sum[1] + f_solC_sum[3] * f_sol_sum[2])

    # 2. Solute-solvent cross-cascades
    # First cascade through solute, second through solvent
    seq_cross1 = sum(f_sol_sum) * sum(f_solv_sum)
    seq_cross1 += sum(f_solC_sum) * sum(f_solvC_sum)

    # 3. Solvent-solute cross-cascades
    # First cascade through solvent, second through solute
    seq_cross2 = sum(f_solv_sum) * sum(f_sol_sum)
    seq_cross2 += sum(f_solvC_sum) * sum(f_solC_sum)

    seq = seq_solute + seq_cross1 + seq_cross2

    # Parallel cascades
    wr = 1.0 / nw
    # Solute-solute parallel
    fPA_total = sum(fPA_sum)
    fPB_total = sum(fPB_sum)
    par_solute = 2 * convolve(fPA_total, fPB_total, mode='same') * wr

    # Solute-solvent parallel cross-cascades
    fPSA_total = sum(fPSA_sum)
    fPSB_total = sum(fPSB_sum)
    par_cross1 = 2 * convolve(fPA_total, fPSB_total, mode='same') * wr
    par_cross2 = 2 * convolve(fPSA_total, fPB_total, mode='same') * wr

    par = par_solute + par_cross1 + par_cross2

    cascade = seq + par[iomega]
    ratio = np.abs(cascade) / np.abs(direct)

    return float(ratio), float(np.abs(cascade)), float(np.abs(direct))
