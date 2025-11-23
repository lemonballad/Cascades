"""
Main simulation script for 2D resonance Raman cascade calculations.

Example usage demonstrating how to run cascade simulations with the
converted Python codebase.
"""

import numpy as np
from cascades.core.basis import basis_tc
from cascades.core.franck_condon import fcinfo_tc
from cascades.core.response import cascade_2drr_res, LaserParameters, MaterialParameters
from cascades.parameters.pna import pna_parameters


def run_2drr_simulation(
    solvent: str = "methanol",
    nmode: int = 3,
    nquanta: int = 4,
    w_L: float = 24000.0,
    dt: float = 10.0,
    nt: int = 256
) -> dict:
    """
    Run a 2D resonance Raman cascade simulation.

    Parameters
    ----------
    solvent : str
        Solvent name for PNA parameters.
    nmode : int
        Number of vibrational modes to include.
    nquanta : int
        Maximum vibrational quanta.
    w_L : float
        Laser frequency in cm^-1.
    dt : float
        Time step in fs.
    nt : int
        Number of time points.

    Returns
    -------
    dict
        Results containing ratio, cascade signal, direct signal,
        and simulation parameters.
    """
    # Load parameters
    params = pna_parameters(solvent)

    # Generate basis set
    base, e_vib = basis_tc(nmode, nquanta, params.wvib)
    print(f"Generated {len(e_vib)} basis states for {nmode} modes with {nquanta} max quanta")

    # Calculate overlap integrals
    ovlp = fcinfo_tc(base, params.disp, nmode, nquanta)
    print(f"Computed overlap integrals: shape {ovlp.shape}")

    # Set up laser parameters
    laser_params = LaserParameters(dt=dt, nt=nt, w_L=w_L)

    # Set up material parameters
    material_params = MaterialParameters(
        gamma_eg=params.gamma_0,
        gamma_vib=10.0,  # Default vibrational dephasing
        weg=params.weg,
        wvib=params.wvib[:nmode]
    )

    # Calculate response
    print("Computing cascade and direct response functions...")
    ratio, cascade, direct = cascade_2drr_res(
        e_vib, nquanta, ovlp, laser_params, material_params
    )

    print(f"Results:")
    print(f"  Cascade signal: {cascade:.6e}")
    print(f"  Direct signal:  {direct:.6e}")
    print(f"  Ratio (C/D):    {ratio:.4f}")

    return {
        "ratio": ratio,
        "cascade": cascade,
        "direct": direct,
        "base": base,
        "e_vib": e_vib,
        "ovlp": ovlp,
        "params": params,
        "laser_params": laser_params,
        "material_params": material_params
    }


def run_detuning_scan(
    solvent: str = "methanol",
    nmode: int = 3,
    nquanta: int = 4,
    detunings: np.ndarray = None,
    dt: float = 10.0,
    nt: int = 256
) -> dict:
    """
    Scan cascade/direct ratio as a function of laser detuning.

    Parameters
    ----------
    solvent : str
        Solvent name for PNA parameters.
    nmode : int
        Number of vibrational modes.
    nquanta : int
        Maximum vibrational quanta.
    detunings : np.ndarray, optional
        Array of detuning values (w_L - w_eg) in cm^-1.
        Default: -2000 to 2000 in 50 steps.
    dt : float
        Time step in fs.
    nt : int
        Number of time points.

    Returns
    -------
    dict
        Results containing detunings and corresponding ratios.
    """
    if detunings is None:
        detunings = np.linspace(-2000, 2000, 21)

    # Load parameters
    params = pna_parameters(solvent)

    # Generate basis and overlaps (only need to do once)
    base, e_vib = basis_tc(nmode, nquanta, params.wvib)
    ovlp = fcinfo_tc(base, params.disp, nmode, nquanta)

    # Material parameters (constant)
    material_params = MaterialParameters(
        gamma_eg=params.gamma_0,
        gamma_vib=10.0,
        weg=params.weg,
        wvib=params.wvib[:nmode]
    )

    # Scan over detunings
    ratios = np.zeros(len(detunings))
    cascades = np.zeros(len(detunings))
    directs = np.zeros(len(detunings))

    for i, det in enumerate(detunings):
        w_L = params.weg + det
        laser_params = LaserParameters(dt=dt, nt=nt, w_L=w_L)

        ratio, cascade, direct = cascade_2drr_res(
            e_vib, nquanta, ovlp, laser_params, material_params
        )

        ratios[i] = ratio
        cascades[i] = cascade
        directs[i] = direct

        print(f"Detuning {det:+7.0f} cm^-1: ratio = {ratio:.4f}")

    return {
        "detunings": detunings,
        "ratios": ratios,
        "cascades": cascades,
        "directs": directs
    }


if __name__ == "__main__":
    # Example: single simulation
    print("=" * 60)
    print("Running single 2DRR simulation for PNA in methanol")
    print("=" * 60)

    results = run_2drr_simulation(
        solvent="methanol",
        nmode=3,
        nquanta=3,
        w_L=24000.0
    )

    print("\n" + "=" * 60)
    print("Running detuning scan")
    print("=" * 60)

    scan_results = run_detuning_scan(
        solvent="methanol",
        nmode=2,
        nquanta=3,
        detunings=np.linspace(-1000, 1000, 5)
    )
