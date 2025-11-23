# Cascades

![Figure 1](https://github.com/user-attachments/assets/697d605d-a020-4467-a59c-32c27f4616ac)

Computational spectroscopy project investigating cascade artifacts in two-dimensional resonance Raman (2DRR) spectroscopy using simulated data.

## Overview

This repository simulates direct and cascade signals to understand how third-order nonlinear optical artifacts contaminate spectroscopic measurements. Available in both MATLAB and Python.

## Implementations

### Python (Recommended)

Modern, well-documented Python package with type hints and test coverage.

```bash
cd python
pip install -e .
```

```python
from cascades.simulations import run_2drr_simulation

results = run_2drr_simulation(
    solvent="methanol",
    nmode=3,
    nquanta=4
)
print(f"Cascade/Direct ratio: {results['ratio']:.4f}")
```

**Features:**
- Resonant 2DRR, FSRS, and off-resonance calculations
- 5 solvents for p-nitroaniline, myoglobin parameters
- Visualization utilities
- Comprehensive test suite (44 tests)
- Sphinx documentation ready

See `python/README.md` for full documentation.

### MATLAB

Original implementation requiring MATLAB 2016+.

Located in `matlab/` with sub-tasks for different analyses.

## Folder Structure

```
Cascades/
├── python/           # Python implementation
│   ├── cascades/     # Main package
│   ├── tests/        # Test suite
│   └── docs/         # Sphinx documentation
├── matlab/           # Original MATLAB code
│   └── sub_task_*/   # Analysis sub-tasks
└── CLAUDE.md         # Project context for AI assistants
```

## Dependencies

### Python
- numpy >= 1.20
- scipy >= 1.7
- matplotlib >= 3.4
- pytest >= 7.0 (testing)

### MATLAB
- MATLAB 2016 or later

## Scientific Background

- **2DRR**: Two-dimensional resonance Raman spectroscopy
- **FSRS**: Femtosecond Stimulated Raman Scattering
- **Cascades**: Third-order nonlinear optical artifacts

## References

### Primary Citation

T. P. Cheshire and A. M. Moran, "Susceptibility of two-dimensional resonance Raman spectroscopies to cascades involving solute and solvent molecules," J. Chem. Phys. **151**, 104203 (2019). [https://doi.org/10.1063/1.5115401](https://doi.org/10.1063/1.5115401)

### Related Work

- A. M. Moran and A. M. Kelley, "Solvent effects on ground and excited electronic state structures of p-nitroaniline," J. Chem. Phys. **115**, 912 (2001). [https://doi.org/10.1063/1.1378319](https://doi.org/10.1063/1.1378319)
- B. P. Molesky, Z. Guo, T. P. Cheshire, and A. M. Moran, "Two-Dimensional Resonance Raman Spectroscopy of Oxygen- and Water-Ligated Myoglobin," J. Chem. Phys. **145**, 034203 (2016). [https://doi.org/10.1063/1.4958625](https://doi.org/10.1063/1.4958625)
- A. B. Myers, R. A. Mathies, D. J. Tannor, and E. J. Heller, "Excited state geometry changes from preresonance Raman intensities: Isoprene and hexatriene," J. Chem. Phys. **77**, 3857 (1982). [https://doi.org/10.1063/1.444339](https://doi.org/10.1063/1.444339)

## Status

Archive project - completed and documented for historical reference.

## License

None specified.
