"""
Myoglobin spectroscopic parameters.

Parameters from: B. P. Molesky, Z. Guo, T. P. Cheshire, A. M. Moran,
"Two-Dimensional Resonance Raman Spectroscopy of Oxygen- and Water-Ligated
Myoglobin" J. Chem. Phys., 145, 034203 (2016)
"""

import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass


@dataclass
class MyoglobinParameters:
    """Spectroscopic parameters for myoglobin."""

    gamma_vib: float             # Vibrational dephasing rate (cm^-1)
    gamma_eg: float              # Electronic homog FWHM (cm^-1)
    weg: float                   # Electronic origin (cm^-1)
    wvib: NDArray[np.float64]    # Vibrational frequencies (cm^-1)
    disp: NDArray[np.float64]    # Dimensionless displacements
    mu_eg: float                 # Transition dipole (Debye)
    n_w_t: float                 # Refractive index
    l: float                     # Path length (m)
    w_t: float                   # Signal frequency (cm^-1)


def myoglobin_parameters() -> MyoglobinParameters:
    """
    Get spectral simulation parameters for myoglobin.

    Returns
    -------
    MyoglobinParameters
        Dataclass containing all spectroscopic parameters.

    References
    ----------
    B. P. Molesky, Z. Guo, T. P. Cheshire, A. M. Moran,
    "Two-Dimensional Resonance Raman Spectroscopy of Oxygen- and
    Water-Ligated Myoglobin" J. Chem. Phys., 145, 034203 (2016)
    """
    return MyoglobinParameters(
        gamma_vib=10,
        gamma_eg=750,
        weg=23250,
        wvib=np.array([220, 370, 674, 1356], dtype=np.float64),
        disp=np.array([0.47, 0.20, 0.26, 0.34], dtype=np.float64),
        mu_eg=8.8,
        n_w_t=1.39,
        l=0.22e-3,
        w_t=23250,  # Same as weg
    )
