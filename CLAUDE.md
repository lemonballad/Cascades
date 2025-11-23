# CLAUDE.md

## Project Overview

**Cascades** is a computational spectroscopy research project investigating cascade artifacts in two-dimensional resonance Raman (2DRR) spectroscopy. The codebase simulates direct and cascade signals to understand how these third-order nonlinear optical artifacts contaminate spectroscopic measurements.

**Status:** Archive project - completed and documented for historical reference.

## Technology Stack

- **Primary Language:** Python 3.9+ with numpy/scipy (recommended)
- **Legacy:** MATLAB (v2016 or later required)
- **IDE:** JetBrains PyCharm/IntelliJ (configured for Python 3.9)

## Primary Citation

T. P. Cheshire and A. M. Moran, "Susceptibility of two-dimensional resonance Raman spectroscopies to cascades involving solute and solvent molecules," J. Chem. Phys. 151, 104203 (2019). https://doi.org/10.1063/1.5115401

## Directory Structure

```
Cascades/
├── README.md                    # Project documentation
├── python/                      # Python implementation (recommended)
│   ├── cascades/                # Main package
│   │   ├── core/                # Basis, Franck-Condon, response functions
│   │   ├── parameters/          # PNA, myoglobin parameters
│   │   ├── visualization/       # Plotting utilities
│   │   └── simulations/         # Main entry points
│   ├── tests/                   # Test suite (44 tests)
│   └── docs/                    # Sphinx documentation
└── matlab/                      # Legacy MATLAB code (248 .m files)
    ├── sub_task_1_consolidate/  # Main 2DRR/FSRS simulation (production code)
    │   └── fig_pics/            # Generated figures
    ├── sub_task_2_i3_cascade/   # I3 cascade analysis
    │   ├── sub_task_2_1/        # Enhanced cascade analysis
    │   └── sub_task_2_2/        # Additional cascade variants
    ├── sub_task_3/              # TC model analysis
    │   └── sub_task_3_!/        # TC variants
    ├── sub_task_4/              # Basic simulation framework
    ├── sub_task_5/              # Extended TC simulations
    │   └── sub_task_5_1/        # TC variant analysis
    └── sub_task_6/              # Additional simulations
```

## Key Scientific Concepts

- **2DRR:** Two-dimensional resonance Raman spectroscopy
- **FSRS:** Femtosecond Stimulated Raman Scattering
- **Cascades:** Third-order nonlinear optical artifacts contaminating measurements
- **Franck-Condon factors:** Vibrational overlap integrals
- **Model systems:** p-Nitroaniline (PNA), Myoglobin in various solvents

## Core Function Modules

### Computation Functions
- `basis_TC.m` - Create vibrational quantum basis states
- `fcfac2_TC.m` - Calculate Franck-Condon overlap integrals
- `fcinfo_TC.m` - Compute overlap integral products
- `cascade_2dRR_Res.m` - Calculate cascade response functions
- `response2.m` - General response function computation

### Parameter Functions
- `PNA_parameters.m` - p-Nitroaniline parameters (6 solvents)
- `myoglobin_parameters.m` - Protein system parameters

### Visualization Functions
- `Plot_2dRR*.m` - 2D Raman spectrum plotting
- `Plot_FSRS*.m` - FSRS spectrum plotting
- `TEST_PLOTS*.m` - Testing and validation plots
- `DEFENSE_*.m` - Publication-ready figures

## Coding Conventions

### Naming
- **Functions:** Snake_case with descriptive names (`cascade_2dRR_Res`, `PNA_parameters`)
- **Variables:** Lowercase with physics acronyms (`kT`, `gamma_vib`, `w_t`, `weg`)
- **Scripts:** Descriptive prefixes (`main_*.m`, `Plot_*.m`, `TEST_*.m`)

### Code Style
- Clear workspace at start (`clear`, `clc`)
- Section headers using `%% Section Name`
- Function headers with purpose, inputs, outputs documentation
- Physical units documented inline (frequencies in cm^-1, time in fs)

### File Organization
- One primary function per file
- Related functions grouped by functionality
- Parameter abstractions in separate `*_parameters.m` files

## Running Simulations

### Python (Recommended)

```bash
cd python
pip install -e .
python -c "from cascades.simulations import run_2drr_simulation; print(run_2drr_simulation('methanol'))"
```

Run tests: `python -m pytest tests/ -v`

### MATLAB (Legacy)

1. Open MATLAB (v2016+)
2. Navigate to desired sub-task directory
3. Run main scripts: `main_*.m`
4. Results saved to `.mat` files
5. Generate plots using `Plot_*.m` or `TEST_PLOTS*.m` scripts

## Common Tasks

### Adding a New Solvent
1. Edit `PNA_parameters.m` in the relevant sub-task
2. Add solvent-specific parameters (reorganization energy, dephasing rates, etc.)
3. Update main simulation script to include new solvent case

### Running Different Model Systems
- p-Nitroaniline: Use `PNA_parameters.m`
- Myoglobin: Use `myoglobin_parameters.m`

### Creating Publication Figures
- Use `DEFENSE_*.m` scripts for publication-ready output
- Figures exported as `.tif` files to `fig_pics/` directories

## Testing & Validation

No automated test framework. Validation is performed through:
- Visual inspection of generated spectral plots
- Comparison with published experimental data
- `TEST_PLOTS*.m` scripts for sanity checks

## Data Files

- **`.mat` files:** Binary MATLAB data caching simulation results
- **`.tif` files:** Generated figure outputs
- Located within respective sub-task directories

## Dependencies

### Python (Recommended)
- numpy >= 1.20
- scipy >= 1.7
- matplotlib >= 3.4
- pytest >= 7.0 (testing)
- numba >= 0.56 (optional, for performance)

### MATLAB (Legacy)
- MATLAB 2016 or later
- Standard MATLAB toolboxes (no external dependencies)

## Notes

- All frequency values are in cm^-1
- All time values are in femtoseconds (fs)
- Energies typically expressed in kT units
- Complex number operations used extensively for response functions
