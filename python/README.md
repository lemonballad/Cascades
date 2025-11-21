# Cascades - Python Implementation

Python implementation of cascade artifact simulations for 2D resonance Raman spectroscopy.

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
│   └── response.py    # Cascade response functions
├── parameters/        # Molecular system parameters
│   ├── pna.py         # p-Nitroaniline parameters
│   └── myoglobin.py   # Myoglobin parameters
├── visualization/     # Plotting utilities
│   └── plots.py       # 2D spectrum plotting
└── simulations/       # Main simulation scripts
    └── run_2drr.py    # 2DRR simulation entry point
```

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

- A. M. Moran, A. M. Kelley. "Solvent effects on ground and excited electronic state structures of p-nitroaniline."
- B. P. Molesky, Z. Guo, T. P. Cheshire, A. M. Moran, "Two-Dimensional Resonance Raman Spectroscopy of Oxygen- and Water-Ligated Myoglobin" J. Chem. Phys., 145, 034203 (2016)
- Myers et al. JCP 77, 3857 (1982) - Franck-Condon overlap integrals
