"""
Cascades: Simulation of cascade artifacts in 2D resonance Raman spectroscopy.

This package provides computational tools for simulating direct and cascade
signals in two-dimensional resonance Raman (2DRR) and Femtosecond Stimulated
Raman Scattering (FSRS) spectroscopy.
"""

from cascades.core import basis, franck_condon, response
from cascades.parameters import pna

__version__ = "0.1.0"
__all__ = ["basis", "franck_condon", "response", "pna"]
