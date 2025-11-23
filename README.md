# Cascades

![Figure 1](https://github.com/user-attachments/assets/697d605d-a020-4467-a59c-32c27f4616ac)

*Cascade artifacts in 2D resonance Raman spectroscopy. Sequential cascades (left) involve signal re-absorption by a second molecule. Parallel cascades (right) involve simultaneous emission from two molecules. From Cheshire & Moran, J. Chem. Phys. 151, 104203 (2019).*

Computational spectroscopy project investigating cascade artifacts in two-dimensional resonance Raman (2DRR) spectroscopy.

## Overview

In nonlinear optical spectroscopy, **cascade artifacts** occur when the signal field radiated by one molecule induces a four-wave mixing process in a second molecule. These third-order artifacts can contaminate fifth-order 2DRR measurements, leading to misinterpretation of molecular dynamics.

This repository provides computational tools to:
- Calculate cascade-to-direct signal ratios for different experimental conditions
- Simulate sequential and parallel cascade pathways
- Model solute-solute and solute-solvent cascades
- Evaluate the susceptibility of 2DRR and FSRS techniques to these artifacts

Available in both Python (recommended) and MATLAB implementations.

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

### Key Concepts

- **2DRR**: Two-dimensional resonance Raman spectroscopy - a fifth-order nonlinear technique that probes vibrational dynamics and structural heterogeneity
- **FSRS**: Femtosecond Stimulated Raman Scattering - ultrafast vibrational spectroscopy with temporal resolution
- **Cascades**: Third-order artifacts where signal from one molecule drives a response in another

### Cascade Types

1. **Sequential cascades**: The radiated field from one molecule is absorbed by a second molecule, which then radiates a new signal
2. **Parallel cascades**: Two molecules simultaneously emit signals that combine at the detector

### Key Findings from the Paper

- Parallel cascades involving two solute molecules can exceed the desired 2DRR signal when mode displacements are ≤1.0
- Solute-solvent cascades can be significant despite 4-6 orders of magnitude smaller Raman cross-sections, due to concentration differences
- The cascade-to-direct ratio depends strongly on laser detuning from electronic resonance

### Model Systems

- **p-Nitroaniline (PNA)**: Intramolecular charge-transfer molecule in various solvents (cyclohexane, dioxane, dichloromethane, acetonitrile, methanol)
- **Myoglobin**: Heme protein system for biological applications

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
