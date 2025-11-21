"""Tests for basis state generation functions."""

import numpy as np
import pytest
from cascades.core.basis import rcbasis, basis_tc


class TestRcbasis:
    """Tests for recursive basis generation."""

    def test_single_mode_single_quantum(self):
        """Single mode with one quantum returns [[1]]."""
        result = rcbasis(1, 1, 1)
        expected = np.array([[1]])
        np.testing.assert_array_equal(result, expected)

    def test_single_mode_zero_quanta(self):
        """Single mode with zero quanta returns [[0]]."""
        result = rcbasis(0, 1, 0)
        expected = np.array([[0]])
        np.testing.assert_array_equal(result, expected)

    def test_two_modes_two_quanta(self):
        """Two modes with two quanta returns correct partitions."""
        result = rcbasis(2, 2, 2)
        expected = np.array([[0, 2], [1, 1], [2, 0]])
        np.testing.assert_array_equal(result, expected)

    def test_target_exceeds_nquanta(self):
        """Returns empty when target exceeds nquanta for single mode."""
        result = rcbasis(1, 1, 2)
        assert result.size == 0

    def test_three_modes(self):
        """Three modes produces correct number of partitions."""
        result = rcbasis(2, 3, 2)
        # Partitions of 2 into 3 parts: (0,0,2), (0,1,1), (0,2,0), (1,0,1), (1,1,0), (2,0,0)
        assert result.shape[0] == 6
        assert result.shape[1] == 3
        assert np.all(np.sum(result, axis=1) == 2)


class TestBasisTC:
    """Tests for basis_tc function."""

    def test_basic_output_shape(self):
        """Output has correct shape for given parameters."""
        wvib = np.array([100.0, 200.0, 300.0])
        base, wviball = basis_tc(2, 2, wvib)

        # nstates = C(nquanta + nmode, nmode) = C(4, 2) = 6
        assert base.shape == (6, 2)
        assert wviball.shape == (6,)

    def test_ground_state_energy(self):
        """Ground state (all zeros) has correct zero-point energy."""
        wvib = np.array([100.0, 200.0])
        base, wviball = basis_tc(2, 1, wvib)

        # Ground state is [0, 0], energy = 0.5*100 + 0.5*200 = 150
        ground_idx = np.where(np.all(base == 0, axis=1))[0][0]
        expected_energy = 0.5 * 100 + 0.5 * 200
        np.testing.assert_almost_equal(wviball[ground_idx], expected_energy)

    def test_nmode_exceeds_wvib_length(self):
        """Raises error when nmode > len(wvib)."""
        wvib = np.array([100.0])
        with pytest.raises(AssertionError):
            basis_tc(2, 1, wvib)

    def test_empty_wvib(self):
        """Raises error for empty wvib."""
        wvib = np.array([])
        with pytest.raises(AssertionError):
            basis_tc(1, 1, wvib)

    def test_energy_ordering(self):
        """Higher quantum states have higher energies."""
        wvib = np.array([100.0])
        base, wviball = basis_tc(1, 3, wvib)

        # Should be monotonically increasing
        sorted_energies = np.sort(wviball)
        np.testing.assert_array_equal(wviball, sorted_energies)
