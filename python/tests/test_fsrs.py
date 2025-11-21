"""Tests for FSRS response function calculations."""

import numpy as np
import pytest
from cascades.core.fsrs import cascade_fsrs_res, FSRSLaserParameters, FSRSMaterialParameters
from cascades.core.basis import basis_tc
from cascades.core.franck_condon import fcinfo_tc


class TestFSRSLaserParameters:
    """Tests for FSRSLaserParameters dataclass."""

    def test_creation(self):
        """Can create FSRSLaserParameters with valid values."""
        params = FSRSLaserParameters(
            w_ap=24000.0, w_rp=12000.0,
            LAMBDA_ap=100.0, LAMBDA_rp=50.0,
            dt=10.0, nt=256
        )
        assert params.w_ap == 24000.0
        assert params.w_rp == 12000.0
        assert params.LAMBDA_ap == 100.0


class TestFSRSMaterialParameters:
    """Tests for FSRSMaterialParameters dataclass."""

    def test_creation(self):
        """Can create FSRSMaterialParameters with valid values."""
        params = FSRSMaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=1000.0
        )
        assert params.gamma_eg == 1000.0
        assert params.wvib == 1000.0


class TestCascadeFSRSRes:
    """Tests for FSRS cascade response function."""

    @pytest.fixture
    def simple_setup(self):
        """Create a simple test setup."""
        nmode = 1
        nquanta = 2
        wvib = np.array([1000.0])
        disp = np.array([0.5])

        base, e_vib = basis_tc(nmode, nquanta, wvib)
        ovlp = fcinfo_tc(base, disp, nmode, nquanta)

        laser_params = FSRSLaserParameters(
            w_ap=24000.0, w_rp=12000.0,
            LAMBDA_ap=100.0, LAMBDA_rp=50.0,
            dt=10.0, nt=64
        )
        material_params = FSRSMaterialParameters(
            gamma_eg=1000.0,
            gamma_vib=10.0,
            weg=23000.0,
            wvib=1000.0
        )

        return e_vib, nquanta, ovlp, laser_params, material_params

    def test_returns_three_values(self, simple_setup):
        """Function returns ratio, cascade, and direct values."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        result = cascade_fsrs_res(e_vib, nquanta, ovlp, laser_params, material_params)
        assert len(result) == 3

    def test_positive_signals(self, simple_setup):
        """Cascade and direct signals are positive."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_fsrs_res(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        assert cascade > 0
        assert direct > 0

    def test_ratio_is_float(self, simple_setup):
        """Ratio is a single float value."""
        e_vib, nquanta, ovlp, laser_params, material_params = simple_setup
        ratio, cascade, direct = cascade_fsrs_res(
            e_vib, nquanta, ovlp, laser_params, material_params
        )
        assert isinstance(ratio, float)
