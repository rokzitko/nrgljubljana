"""Analytic checks for the QSTZ Wick oracle and the finite-space contract."""

import unittest

import numpy as np

from validate_qstz_spin import FIXTURE, MATS, check_retention, compare, dimensions, inspect_fixture, spin_matsubara


class QSTZSpinTests(unittest.TestCase):
    def test_atomic_sign_and_normalization(self):
        for temperature in (0.2, 0.7):
            f = 1 / (1 + np.exp(0.15 / temperature))
            expected = [-f * (1 - f) / (2 * temperature), 0, 0, 0]
            compare(spin_matsubara(np.array([[0.15]]), temperature), expected, "free atomic spin")

    def test_degenerate_levels(self):
        # At zero energy the atomic four-state ensemble has <Sz^2>=1/8.
        compare(spin_matsubara(np.zeros((7, 7)), 0.2), [-1 / (8 * 0.2), 0, 0, 0], "degenerate atom")

    def test_symmetric_dimer(self):
        v = 0.3
        for temperature in (0.2, 0.7):
            omega = 2 * np.pi * temperature * np.arange(MATS)
            expected = -v * np.tanh(v / (2 * temperature)) / (2 * (omega ** 2 + 4 * v ** 2))
            f = 1 / (1 + np.exp(v / temperature))
            expected[0] -= f * (1 - f) / (4 * temperature)
            compare(spin_matsubara(np.array([[0, v], [v, 0]]), temperature), expected, "dimer")

    def test_prepared_seed_and_spectators(self):
        _, h = inspect_fixture(FIXTURE)
        self.assertEqual(h.shape, (7, 7))
        for temperature in (0.2, 0.7):
            full = spin_matsubara(h, temperature)
            central = spin_matsubara(h[np.ix_([0, 2, 5], [0, 2, 5])], temperature)
            compare(full, central, "spectator cancellation")
            self.assertTrue(np.all(full.real < 0))
            self.assertGreater(abs(full[1]), 1e-4)

    def test_dimensions(self):
        for charges, multiplets, states in (((0, 1, 0, -1), 126, 256),
                                            ((0, 1, 0, -1, 1, 0, -1), 6435, 16384)):
            counts = dimensions(charges)
            self.assertEqual(sum(counts.values()), multiplets)
            self.assertEqual(sum(ss * count for (_, ss, _), count in counts.items()), states)

    def test_retention_and_numeric_failures(self):
        log = ("Number of states (multiplicity taken into account): 256\n"
               "[NRG] Iteration 1/1\n"
               "nrkept=6435 nrkeptmult=16384 nrall=6435 nrallmult=16384\n"
               "[DM] Iteration 1/1\n")
        check_retention(log)
        for bad in ("", log.replace("nrkept=6435", "nrkept=100"), log.replace("[DM]", "[NRG]")):
            with self.assertRaises(ValueError):
                check_retention(bad)
        for actual in ([], [np.nan], [np.inf], [1.0]):
            with self.assertRaises(ValueError):
                compare(actual, [0.0], "invalid or mismatched output")


if __name__ == "__main__":
    unittest.main()
