"""Tests for off-resonance response function calculations."""

import numpy as np
import pytest
from cascades.core.offres import cascade_2drr_offres, OffResMaterialParameters
from cascades.core.response import LaserParameters
from cascades.core.basis import basis_tc
from cascades.core.franck_condon import fcinfo_tc


class TestOffResMaterialParameters:
    """Tests for OffResMaterialParameters dataclass."""

    def test_creation(self):
        """Can create OffResMaterialParameters with solvent frequency."""
        params = OffResMaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=np.array([1000.0]),
            wsolv=800.0
        )
        assert params.gamma_eg == 1000.0
        assert params.wsolv == 800.0

    def test_inherits_from_material_parameters(self):
        """OffResMaterialParameters includes base material parameters."""
        params = OffResMaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=np.array([1000.0]),
            wsolv=800.0
        )
        assert hasattr(params, 'gamma_eg')
        assert hasattr(params, 'wsolv')


class TestCascade2DRROffRes:
    """Tests for off-resonance 2DRR cascade response function."""

    @pytest.fixture
    def simple_setup(self):
        """Create a simple test setup."""
        nmode = 1
        nquanta = 2
        wvib = np.array([1000.0])
        disp = np.array([0.5])

        base, e_vib = basis_tc(nmode, nquanta, wvib)
        ovlp = fcinfo_tc(base, disp, nmode, nquanta)

        laser_params = LaserParameters(dt=10.0, nt=64, w_L=20000.0)
        material_params = OffResMaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=wvib,
            wsolv=800.0
        )

        return e_vib, nquanta, ovlp, laser_params, material_params

    def test_returns_three_values(self, simple_setup):
        """Function returns ratio, cascade, and direct values."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        result = cascade_2drr_offres(e_vib, nquanta, ovlp, laser_params, material_params)
        assert len(result) == 3

    def test_positive_signals(self, simple_setup):
        """Cascade and direct signals are positive."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_2drr_offres(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        assert cascade >= 0
        assert direct >= 0

    def test_ratio_calculation(self, simple_setup):
        """Ratio equals cascade/direct when direct is non-zero."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_2drr_offres(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        if direct > 0:
            np.testing.assert_almost_equal(ratio, cascade / direct)
