# Cascades - Python Implementation

Python implementation of cascade artifact simulations for 2D resonance Raman spectroscopy.

This package calculates the ratio of cascade artifacts to direct fifth-order signals in 2DRR spectroscopy. It implements the computational methods from [Cheshire & Moran, J. Chem. Phys. 151, 104203 (2019)](https://doi.org/10.1063/1.5115401).

## Features

- **Resonant 2DRR**: Sequential and parallel cascade calculations with full Liouville pathway decomposition
- **FSRS**: Femtosecond Stimulated Raman Scattering response functions
- **Off-resonance**: Solute-solvent cascade contributions
- **Franck-Condon**: Overlap integrals using the Heller recursive algorithm
- **Optimization**: Optional Numba JIT compilation for 10-50x speedup
- **Type hints**: Full type annotations with numpy.typing

## Installation

```bash
cd python
pip install -e .
```

Or install dependencies only:

```bash
pip install -r requirements.txt
```

## Quick Start

```python
from cascades.simulations import run_2drr_simulation

# Run a simulation for PNA in methanol
results = run_2drr_simulation(
    solvent="methanol",
    nmode=3,
    nquanta=4,
    w_L=24000.0
)

print(f"Cascade/Direct ratio: {results['ratio']:.4f}")
```

## Package Structure

```
cascades/
├── core/              # Core computation modules
│   ├── basis.py       # Vibrational basis state generation
│   ├── franck_condon.py  # Overlap integral calculations
│   ├── response.py    # Resonant 2DRR cascade response functions
│   ├── fsrs.py        # FSRS (Femtosecond Stimulated Raman) response
│   ├── offres.py      # Off-resonance 2DRR with solvent effects
│   └── optimized.py   # Numba-accelerated functions (10-50x speedup)
├── parameters/        # Molecular system parameters
│   ├── pna.py         # p-Nitroaniline parameters
│   └── myoglobin.py   # Myoglobin parameters
├── visualization/     # Plotting utilities
│   └── plots.py       # 2D spectrum plotting
└── simulations/       # Main simulation scripts
    └── run_2drr.py    # 2DRR simulation entry point
```

## Performance Optimization

For 10-50x speedup on larger basis sets, install with Numba:

```bash
pip install cascades[fast]
# or: pip install numba>=0.56
```

Use the optimized functions:

```python
from cascades.core import fcfac2_tc_fast, fcinfo_tc_fast, check_numba_available

if check_numba_available():
    ovlp = fcinfo_tc_fast(base, disp, nmode, nquanta)  # JIT-compiled
else:
    ovlp = fcinfo_tc(base, disp, nmode, nquanta)  # Pure numpy fallback
```

## Available Spectroscopy Types

- **Resonant 2DRR**: `cascade_2drr_res` - Two-dimensional resonance Raman
- **FSRS**: `cascade_fsrs_res` - Femtosecond Stimulated Raman Scattering
- **Off-resonance 2DRR**: `cascade_2drr_offres` - Off-resonance with solvent contributions

## Usage Examples

### Basic Simulation

```python
import numpy as np
from cascades.core.basis import basis_tc
from cascades.core.franck_condon import fcinfo_tc
from cascades.core.response import cascade_2drr_res, LaserParameters, MaterialParameters
from cascades.parameters.pna import pna_parameters

# Load parameters for PNA in acetonitrile
params = pna_parameters("acetonitrile")

# Generate vibrational basis
nmode = 3
nquanta = 4
base, e_vib = basis_tc(nmode, nquanta, params.wvib)

# Calculate Franck-Condon overlaps
ovlp = fcinfo_tc(base, params.disp, nmode, nquanta)

# Set up simulation parameters
laser = LaserParameters(dt=10.0, nt=256, w_L=23000.0)
material = MaterialParameters(
    gamma_eg=params.gamma_0,
    gamma_vib=10.0,
    weg=params.weg,
    wvib=params.wvib[:nmode]
)

# Calculate cascade and direct signals
ratio, cascade, direct = cascade_2drr_res(e_vib, nquanta, ovlp, laser, material)
```

### FSRS Simulation

```python
from cascades.core.fsrs import cascade_fsrs_res, FSRSLaserParameters, FSRSMaterialParameters
from cascades.core.basis import basis_tc
from cascades.core.franck_condon import fcinfo_tc

# Set up basis
nmode, nquanta = 2, 3
wvib = np.array([1000.0, 1200.0])
disp = np.array([0.5, 0.3])
base, e_vib = basis_tc(nmode, nquanta, wvib)
ovlp = fcinfo_tc(base, disp, nmode, nquanta)

# FSRS parameters
laser = FSRSLaserParameters(
    w_ap=24000.0, w_rp=12000.0,
    LAMBDA_ap=100.0, LAMBDA_rp=50.0,
    dt=10.0, nt=256
)
material = FSRSMaterialParameters(
    gamma_eg=1000.0, gamma_vib=10.0,
    weg=23000.0, wvib=1000.0
)

ratio, cascade, direct = cascade_fsrs_res(e_vib, nquanta, ovlp, laser, material)
```

### Detuning Scan

```python
from cascades.simulations import run_detuning_scan
import numpy as np

results = run_detuning_scan(
    solvent="methanol",
    nmode=2,
    nquanta=3,
    detunings=np.linspace(-2000, 2000, 41)
)

# Plot results
import matplotlib.pyplot as plt
plt.plot(results['detunings'], results['ratios'])
plt.xlabel('Detuning (cm$^{-1}$)')
plt.ylabel('Cascade/Direct Ratio')
plt.show()
```

### Available Solvents

p-Nitroaniline parameters are available for:
- `cyclohexane`
- `1,4-dioxane`
- `dichloromethane`
- `acetonitrile`
- `methanol`

## Running Tests

```bash
cd python
pytest tests/ -v
```

## Physical Units

- Frequencies: cm^-1
- Time: femtoseconds (fs)
- Energies: cm^-1 (kT = 200 cm^-1)
- Displacements: unitless

## Dependencies

- numpy >= 1.20
- scipy >= 1.7
- matplotlib >= 3.4
- pytest >= 7.0 (for testing)

## References

### Primary Citation

T. P. Cheshire and A. M. Moran, "Susceptibility of two-dimensional resonance Raman spectroscopies to cascades involving solute and solvent molecules," J. Chem. Phys. **151**, 104203 (2019). [https://doi.org/10.1063/1.5115401](https://doi.org/10.1063/1.5115401)

### Related Work

- A. M. Moran and A. M. Kelley, "Solvent effects on ground and excited electronic state structures of p-nitroaniline," J. Chem. Phys. **115**, 912 (2001). [https://doi.org/10.1063/1.1378319](https://doi.org/10.1063/1.1378319)
- B. P. Molesky, Z. Guo, T. P. Cheshire, and A. M. Moran, "Two-Dimensional Resonance Raman Spectroscopy of Oxygen- and Water-Ligated Myoglobin," J. Chem. Phys. **145**, 034203 (2016). [https://doi.org/10.1063/1.4958625](https://doi.org/10.1063/1.4958625)
- A. B. Myers, R. A. Mathies, D. J. Tannor, and E. J. Heller, "Excited state geometry changes from preresonance Raman intensities: Isoprene and hexatriene," J. Chem. Phys. **77**, 3857 (1982). [https://doi.org/10.1063/1.444339](https://doi.org/10.1063/1.444339)
