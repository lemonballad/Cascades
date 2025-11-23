"""
FSRS (Femtosecond Stimulated Raman Scattering) response functions.

Compute cascade and direct signals for resonant FSRS spectroscopy.
"""

import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass


@dataclass
class FSRSLaserParameters:
    """Laser parameters for FSRS simulation."""
    w_ap: float        # Actinic pulse frequency (cm^-1)
    w_rp: float        # Raman pulse frequency (cm^-1)
    LAMBDA_ap: float   # Actinic pulse spectral width (cm^-1)
    LAMBDA_rp: float   # Raman pulse spectral width (cm^-1)
    dt: float          # Time step (fs)
    nt: int            # Number of time points


@dataclass
class FSRSMaterialParameters:
    """Material parameters for FSRS simulation."""
    gamma_eg: float           # Electronic dephasing rate (cm^-1)
    gamma_vib: float          # Vibrational dephasing rate (cm^-1)
    weg: float                # Electronic energy gap (cm^-1)
    wvib: float               # Vibrational frequency (cm^-1)


def cascade_fsrs_res(
    e_vib: NDArray[np.float64],
    nquanta: int,
    ovlp: NDArray[np.float64],
    laser_params: FSRSLaserParameters,
    material_params: FSRSMaterialParameters
) -> tuple[float, float, float]:
    """
    Compute cascade and direct response functions for resonant FSRS.

    Parameters
    ----------
    e_vib : NDArray[np.float64]
        Vibrational energies for each basis state (cm^-1).
    nquanta : int
        Number of vibrational quanta.
    ovlp : NDArray[np.float64]
        Overlap integrals between basis states.
    laser_params : FSRSLaserParameters
        Laser parameters (actinic/Raman pulses).
    material_params : FSRSMaterialParameters
        Material parameters.

    Returns
    -------
    ratio : float
        Cascade-to-direct signal ratio.
    cascade : float
        Absolute cascade signal.
    direct : float
        Absolute direct signal.
    """
    # Boltzmann populations
    kT = 200.0  # cm^-1
    boltz_factor = np.exp(-e_vib / kT)
    partition_func = np.sum(boltz_factor)
    boltz_pop = boltz_factor / partition_func

    iq = nquanta
    c = 2.998e-5  # Speed of light (cm/fs)
    r2w = 2 * np.pi * c

    # Vibrational energy gaps
    wi, wf = np.meshgrid(e_vib, e_vib)
    w = wf - wi

    # Unpack parameters
    gamma_eg = material_params.gamma_eg
    gamma_vib = material_params.gamma_vib
    weg = material_params.weg
    wvib = material_params.wvib

    w_ap = laser_params.w_ap
    w_rp = laser_params.w_rp
    LAMBDA_ap = laser_params.LAMBDA_ap
    LAMBDA_rp = laser_params.LAMBDA_rp

    # Damping
    damp = gamma_vib * r2w
    DAMP = damp * np.ones_like(w)
    DAMP[w == 0] = 0

    # Time delays
    tau1 = 100 / c / wvib
    tau2 = 0

    # Signal frequency
    w_t = w_rp - wvib

    # Initialize response functions
    r = [np.zeros(iq, dtype=np.complex128) for _ in range(16)]
    fA = [np.zeros(iq, dtype=np.complex128) for _ in range(4)]
    fB = [np.zeros(iq, dtype=np.complex128) for _ in range(4)]

    # Lineshape functions
    L_ap_p = 1.0 / (w_ap - weg - w + 1j * gamma_eg)
    L_ap_m = 1.0 / (-w_ap + weg - w + 1j * gamma_eg)
    L_rp_m = 1.0 / (-w_rp + weg - w + 1j * gamma_eg)
    L_t = 1.0 / (w_t - weg - w + 1j * gamma_eg)

    # D and J functions
    D0 = np.exp((-1j * w * r2w - DAMP) * tau1) * 2 * LAMBDA_ap / (w**2 + LAMBDA_ap**2)
    D = np.exp((-1j * w * r2w - damp) * tau1) * 2 * LAMBDA_ap / (w**2 + LAMBDA_ap**2)
    J = 1.0 / (w_t - w_rp - w + 1j * (gamma_vib - LAMBDA_rp)) * np.exp((1j * (w_rp - w_t) - LAMBDA_rp) * tau2 * r2w)

    D[w == 0] = D0[w == 0]
    J[w == 0] = 0

    # Main computation loop
    m = 0
    bp_m = boltz_pop[m]

    for n in range(iq):
        for k in range(iq):
            for l in range(iq):
                # Third-order cascade auxiliary functions (A)
                fA[0][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * L_t[n, m] * J[n, k] * L_t[n, l]
                fA[1][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * L_rp_m[m, n] * J[k, n] * L_t[k, l]
                fA[2][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * L_t[n, m] * J[k, m] * L_t[l, m]
                fA[3][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * L_rp_m[m, n] * J[m, k] * L_t[l, k]

                # Third-order cascade auxiliary functions (B)
                fB[0][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, l] * L_ap_p[n, m] * L_t[n, l] * D[n, k]
                fB[1][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, l] * L_ap_m[m, n] * L_t[k, l] * D[k, n]
                fB[2][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, m] * L_ap_p[n, m] * L_t[l, m] * D[k, m]
                fB[3][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, k] * L_ap_m[m, n] * L_t[l, k] * D[m, k]

                # Fifth-order direct terms
                for u in range(iq):
                    for v in range(iq):
                        r[0][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, u] * ovlp[v, u] * ovlp[v, m] * L_ap_p[n, m] * D[k, m] * L_t[l, m] * J[u, m] * L_t[v, m]
                        r[1][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, u] * ovlp[v, k] * ovlp[v, u] * L_ap_p[n, m] * D[k, m] * L_rp_m[k, l] * J[k, u] * L_t[v, u]
                        r[2][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[l, u] * ovlp[v, m] * ovlp[v, u] * L_ap_m[m, n] * D[m, k] * L_rp_m[m, l] * J[m, u] * L_t[v, u]
                        r[3][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[l, u] * ovlp[v, u] * ovlp[v, k] * L_ap_m[m, n] * D[m, k] * L_t[l, k] * J[u, k] * L_t[v, k]

                        r[4][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[u, l] * ovlp[k, v] * ovlp[u, v] * L_ap_p[n, m] * D[n, k] * L_rp_m[l, k] * J[u, k] * L_t[u, v]
                        r[5][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[u, l] * ovlp[u, v] * ovlp[n, v] * L_ap_p[n, m] * D[n, k] * L_t[n, l] * J[n, u] * L_t[n, v]
                        r[6][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[u, l] * ovlp[n, v] * ovlp[u, v] * L_ap_m[m, n] * D[k, n] * L_rp_m[l, n] * J[u, n] * L_t[u, v]
                        r[7][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[u, l] * ovlp[u, v] * ovlp[k, v] * L_ap_m[m, n] * D[k, n] * L_t[k, l] * J[k, u] * L_t[k, v]

                        r[8][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[u, m] * ovlp[u, v] * ovlp[l, v] * L_ap_p[n, m] * D[k, m] * L_t[l, m] * J[l, u] * L_t[l, v]
                        r[9][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[u, k] * ovlp[l, v] * ovlp[u, v] * L_ap_p[n, m] * D[k, m] * L_rp_m[k, l] * J[u, l] * L_t[u, v]
                        r[10][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, m] * ovlp[u, k] * ovlp[u, v] * ovlp[l, v] * L_ap_m[m, n] * D[m, k] * L_t[l, k] * J[l, u] * L_t[l, v]
                        r[11][n] += bp_m * ovlp[n, m] * ovlp[n, k] * ovlp[l, k] * ovlp[u, m] * ovlp[l, v] * ovlp[u, v] * L_ap_m[m, n] * D[m, k] * L_rp_m[m, l] * J[u, l] * L_t[u, v]

                        r[12][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, u] * ovlp[v, l] * ovlp[v, u] * L_ap_p[n, m] * D[n, k] * L_rp_m[l, k] * J[l, u] * L_t[v, u]
                        r[13][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, u] * ovlp[v, u] * ovlp[v, l] * L_ap_p[n, m] * D[n, k] * L_t[n, l] * J[u, l] * L_t[v, l]
                        r[14][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[k, l] * ovlp[n, u] * ovlp[v, l] * ovlp[v, u] * L_ap_m[m, n] * D[k, n] * L_rp_m[l, n] * J[l, u] * L_t[v, u]
                        r[15][n] += bp_m * ovlp[n, m] * ovlp[k, m] * ovlp[n, l] * ovlp[k, u] * ovlp[v, u] * ovlp[v, l] * L_ap_m[m, n] * D[k, n] * L_t[k, l] * J[u, l] * L_t[v, l]

    # Sum over states
    fA_sum = [np.sum(f) for f in fA]
    fB_sum = [np.sum(f) for f in fB]
    r_sum = [np.sum(ri) for ri in r]

    # Cascade signal (sequential)
    cascade = (1j)**6 * sum(fA_sum) * sum(fB_sum)

    # Direct fifth-order signal
    direct = (1j)**5 * sum(r_sum)

    ratio = np.abs(cascade) / np.abs(direct)

    return float(ratio), float(np.abs(cascade)), float(np.abs(direct))
