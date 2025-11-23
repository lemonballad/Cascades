"""
p-Nitroaniline (PNA) spectroscopic parameters.

Parameters from: A. M. Moran, A. M. Kelley. "Solvent effects on ground and
excited electronic state structures of p-nitroaniline."
"""

import numpy as np
from numpy.typing import NDArray
from dataclasses import dataclass
from typing import Literal


Solvent = Literal[
    "cyclohexane",
    "1,4-dioxane",
    "dichloromethane",
    "acetonitrile",
    "methanol"
]


@dataclass
class PNAParameters:
    """Spectroscopic parameters for p-nitroaniline in a given solvent."""

    inhomog: float               # Electronic inhomog std dev (cm^-1)
    gamma_0: float               # Electronic homog FWHM (cm^-1)
    kappa: float                 # LAMBDA/DELTA ratio
    lambda_: float               # Classic reorganization energy: DEL^2/(2kT) (cm^-1)
    weg: float                   # Electronic origin (cm^-1)
    translen: float              # Transition length (Angstroms)
    refrac_index: float          # Refractive index
    wvib: NDArray[np.float64]    # Vibrational frequencies (cm^-1)
    disp: NDArray[np.float64]    # Dimensionless displacements


def pna_parameters(solvent: Solvent) -> PNAParameters:
    """
    Get spectral simulation parameters for p-nitroaniline in given solvent.

    Parameters
    ----------
    solvent : str
        Solvent name. Options: cyclohexane, 1,4-dioxane, dichloromethane,
        acetonitrile, methanol.

    Returns
    -------
    PNAParameters
        Dataclass containing all spectroscopic parameters.

    Raises
    ------
    ValueError
        If solvent is not recognized.

    References
    ----------
    A. M. Moran, A. M. Kelley. "Solvent effects on ground and excited
    electronic state structures of p-nitroaniline."
    """
    params = {
        "cyclohexane": PNAParameters(
            inhomog=0,
            gamma_0=2045,
            kappa=0.1,
            lambda_=1890,
            weg=27542,
            translen=1.1,
            refrac_index=1.427,
            wvib=np.array([859, 1112, 1338, 1498, 1599], dtype=np.float64),
            disp=np.array([0.88, 0.288, 1.433, 0.382, 0.496], dtype=np.float64),
        ),
        "1,4-dioxane": PNAParameters(
            inhomog=110,
            gamma_0=2473,
            kappa=0.05,
            lambda_=2700,
            weg=24100,
            translen=1.189,
            refrac_index=1.422,
            wvib=np.array([861, 1113, 1328, 1508, 1602], dtype=np.float64),
            disp=np.array([1.23, 0.53, 1.17, 0.223, 0.341], dtype=np.float64),
        ),
        "dichloromethane": PNAParameters(
            inhomog=60,
            gamma_0=2915,
            kappa=0.07,
            lambda_=3780,
            weg=23490,
            translen=1.172,
            refrac_index=1.424,
            wvib=np.array([862, 1112, 1326, 1504, 1602, 1625], dtype=np.float64),
            disp=np.array([1.049, 0.231, 1.184, 0.1727, 0.3195, 0.27], dtype=np.float64),
        ),
        "acetonitrile": PNAParameters(
            inhomog=0,
            gamma_0=3000,
            kappa=0.2,
            lambda_=4310,
            weg=22720,
            translen=1.199,
            refrac_index=1.344,
            wvib=np.array([861, 1112, 1326, 1507, 1599], dtype=np.float64),
            disp=np.array([0.84, 0.242, 0.976, 0.152, 0.292], dtype=np.float64),
        ),
        "methanol": PNAParameters(
            inhomog=750,
            gamma_0=1410,
            kappa=0.1,
            lambda_=900,
            weg=23800,
            translen=1.243,
            refrac_index=1.329,
            wvib=np.array([861, 1113, 1328, 1508, 1600], dtype=np.float64),
            disp=np.array([1.46, 0.65, 1.451, 0.266, 0.287], dtype=np.float64),
        ),
    }

    assert solvent in params, f"Invalid solvent: {solvent}. Valid options: {list(params.keys())}"
    return params[solvent]
