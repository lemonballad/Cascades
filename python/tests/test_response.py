"""Tests for response function calculations."""

import numpy as np
import pytest
from cascades.core.response import cascade_2drr_res, LaserParameters, MaterialParameters
from cascades.core.basis import basis_tc
from cascades.core.franck_condon import fcinfo_tc


class TestLaserParameters:
    """Tests for LaserParameters dataclass."""

    def test_creation(self):
        """Can create LaserParameters with valid values."""
        params = LaserParameters(dt=10.0, nt=256, w_L=24000.0)
        assert params.dt == 10.0
        assert params.nt == 256
        assert params.w_L == 24000.0


class TestMaterialParameters:
    """Tests for MaterialParameters dataclass."""

    def test_creation(self):
        """Can create MaterialParameters with valid values."""
        params = MaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=np.array([800.0, 1200.0])
        )
        assert params.gamma_eg == 1000.0
        assert len(params.wvib) == 2


class TestCascade2drrRes:
    """Tests for cascade response function calculation."""

    @pytest.fixture
    def simple_setup(self):
        """Create a simple test setup with minimal parameters."""
        nmode = 1
        nquanta = 2
        wvib = np.array([1000.0])
        disp = np.array([0.5])

        base, e_vib = basis_tc(nmode, nquanta, wvib)
        ovlp = fcinfo_tc(base, disp, nmode, nquanta)

        laser_params = LaserParameters(dt=10.0, nt=64, w_L=24000.0)
        material_params = MaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=wvib
        )

        return e_vib, nquanta, ovlp, laser_params, material_params

    def test_returns_three_values(self, simple_setup):
        """Function returns ratio, cascade, and direct values."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        result = cascade_2drr_res(e_vib, nquanta, ovlp, laser_params, material_params)
        assert len(result) == 3

    def test_positive_signals(self, simple_setup):
        """Cascade and direct signals are positive."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_2drr_res(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        assert cascade > 0
        assert direct > 0

    def test_ratio_is_float(self, simple_setup):
        """Ratio is a single float value."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_2drr_res(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        assert isinstance(ratio, float)

    def test_ratio_calculation(self, simple_setup):
        """Ratio equals cascade/direct."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_2drr_res(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        np.testing.assert_almost_equal(ratio, cascade / direct)
