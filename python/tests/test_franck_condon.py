"""Tests for Franck-Condon overlap integral functions."""

import numpy as np
import pytest
from cascades.core.franck_condon import fcfac2_tc, fcinfo_tc
from cascades.core.basis import basis_tc


class TestFcfac2TC:
    """Tests for Franck-Condon overlap integral calculation."""

    def test_output_shape(self):
        """Output has correct shape (nquanta+1, nquanta+1)."""
        ovlp = fcfac2_tc(d=1.0, nquanta=3)
        assert ovlp.shape == (4, 4)

    def test_zero_displacement(self):
        """Zero displacement gives identity-like matrix."""
        ovlp = fcfac2_tc(d=0.0, nquanta=3)
        # <m|n> = delta_{mn} when d=0
        expected = np.eye(4)
        np.testing.assert_array_almost_equal(ovlp, expected, decimal=10)

    def test_ground_state_overlap(self):
        """<0|0> follows expected formula for non-zero displacement."""
        d = 1.0
        ovlp = fcfac2_tc(d=d, nquanta=2)
        # <0|0> = exp(-d^2/2) = exp(-0.5)
        expected_00 = np.exp(-d**2 / 2)
        np.testing.assert_almost_equal(ovlp[0, 0], expected_00)

    def test_symmetry_property(self):
        """Matrix has expected symmetry properties."""
        ovlp = fcfac2_tc(d=0.5, nquanta=4)
        # For equal frequencies, |<m|n>| = |<n|m>|
        np.testing.assert_array_almost_equal(np.abs(ovlp), np.abs(ovlp.T))

    def test_normalization(self):
        """Sum of squared overlaps approximately equals 1."""
        ovlp = fcfac2_tc(d=0.5, nquanta=10)
        # For sufficiently high nquanta, sum_n |<0|n>|^2 ≈ 1
        sum_sq = np.sum(ovlp[0, :]**2)
        np.testing.assert_almost_equal(sum_sq, 1.0, decimal=3)


class TestFcinfoTC:
    """Tests for product of overlap integrals calculation."""

    def test_output_shape(self):
        """Output has correct shape (nstates, nstates)."""
        wvib = np.array([100.0, 200.0])
        disp = np.array([0.5, 0.3])
        base, _ = basis_tc(2, 2, wvib)

        fcall = fcinfo_tc(base, disp, 2, 2)
        nstates = base.shape[0]
        assert fcall.shape == (nstates, nstates)

    def test_diagonal_elements(self):
        """Diagonal elements are products of self-overlaps."""
        wvib = np.array([100.0])
        disp = np.array([0.0])  # Zero displacement = identity
        base, _ = basis_tc(1, 2, wvib)

        fcall = fcinfo_tc(base, disp, 1, 2)
        # With zero displacement, diagonal should be 1
        np.testing.assert_array_almost_equal(np.diag(fcall), np.ones(3))

    def test_single_mode(self):
        """Single mode matches direct fcfac2_tc calculation."""
        wvib = np.array([100.0])
        disp = np.array([0.5])
        base, _ = basis_tc(1, 3, wvib)

        fcall = fcinfo_tc(base, disp, 1, 3)
        ovlp_direct = fcfac2_tc(0.5, 3)

        # fcall should match ovlp_direct for single mode
        np.testing.assert_array_almost_equal(fcall, ovlp_direct)
