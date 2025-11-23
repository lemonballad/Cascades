"""Tests for parameter modules."""

import numpy as np
import pytest
from cascades.parameters.pna import pna_parameters, PNAParameters
from cascades.parameters.myoglobin import myoglobin_parameters, MyoglobinParameters


class TestPNAParameters:
    """Tests for p-nitroaniline parameters."""

    def test_valid_solvents(self):
        """All valid solvents return PNAParameters."""
        solvents = ["cyclohexane", "1,4-dioxane", "dichloromethane",
                    "acetonitrile", "methanol"]
        for solvent in solvents:
            params = pna_parameters(solvent)
            assert isinstance(params, PNAParameters)

    def test_invalid_solvent(self):
        """Invalid solvent raises error."""
        with pytest.raises(AssertionError):
            pna_parameters("water")

    def test_methanol_values(self):
        """Methanol parameters match expected values."""
        params = pna_parameters("methanol")
        assert params.inhomog == 750
        assert params.gamma_0 == 1410
        assert params.kappa == 0.1
        assert params.lambda_ == 900
        assert params.weg == 23800
        np.testing.assert_almost_equal(params.translen, 1.243)

    def test_wvib_array_type(self):
        """Vibrational frequencies are numpy arrays."""
        params = pna_parameters("cyclohexane")
        assert isinstance(params.wvib, np.ndarray)
        assert isinstance(params.disp, np.ndarray)

    def test_array_lengths_match(self):
        """wvib and disp arrays have same length."""
        for solvent in ["cyclohexane", "methanol", "dichloromethane"]:
            params = pna_parameters(solvent)
            assert len(params.wvib) == len(params.disp)

    def test_positive_values(self):
        """Physical parameters are positive."""
        params = pna_parameters("acetonitrile")
        assert params.gamma_0 > 0
        assert params.weg > 0
        assert params.translen > 0
        assert params.refrac_index > 0


class TestMyoglobinParameters:
    """Tests for myoglobin parameters."""

    def test_returns_dataclass(self):
        """Function returns MyoglobinParameters dataclass."""
        params = myoglobin_parameters()
        assert isinstance(params, MyoglobinParameters)

    def test_expected_values(self):
        """Parameters match expected values from literature."""
        params = myoglobin_parameters()
        assert params.gamma_vib == 10
        assert params.gamma_eg == 750
        assert params.weg == 23250
        assert params.mu_eg == 8.8

    def test_array_lengths(self):
        """wvib and disp arrays have same length."""
        params = myoglobin_parameters()
        assert len(params.wvib) == len(params.disp)
        assert len(params.wvib) == 4

    def test_signal_frequency(self):
        """Signal frequency equals electronic gap."""
        params = myoglobin_parameters()
        assert params.w_t == params.weg
