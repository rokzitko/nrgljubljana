"""Analytic checks for the QSTZ Wick oracle and the finite-space contract."""

import unittest

import numpy as np

from ed_siam import AndersonModel, solve
from validate_qstz_spin import FIXTURE, MATS, check_chit, check_retention, compare, dimensions, inspect_fixture, spin_matsubara


def lehmann_chit(energies, operator, temperature):
    """T*chi from a full many-body eigenbasis, including degenerate matrix elements."""
    energies = np.asarray(energies)
    p = np.exp(-(energies - energies.min()) / temperature)
    p /= p.sum()
    result = 0.0
    for i, ei in enumerate(energies):
        for j, ej in enumerate(energies):
            weight = p[i] if abs(ei - ej) < 1e-12 else temperature * (p[i] - p[j]) / (ej - ei)
            result += abs(operator[i, j]) ** 2 * weight
    return result


class QSTZSpinTests(unittest.TestCase):
    def test_free_spin_curie_normalization(self):
        for spin, expected in ((0.5, 0.25), (1, 2 / 3)):
            sz = np.diag(np.arange(-spin, spin + 1))
            energies = np.zeros(len(sz))
            # A rotated spin component puts the weight off diagonal within the degenerate multiplet.
            rotation = np.fft.fft(np.eye(len(sz))) / np.sqrt(len(sz))
            for operator in (sz, rotation.conj().T @ sz @ rotation):
                for temperature in (0.2, 0.7):
                    compare(lehmann_chit(energies, operator, temperature), expected, "free-spin T*chi")
                    compare(lehmann_chit(energies, 2 * operator, temperature), 4 * expected, "Pauli normalization")

    def test_exchange_pair_lehmann(self):
        j = 0.43
        sx = np.array([[0, 0.5], [0.5, 0]])
        sy = np.array([[0, -0.5j], [0.5j, 0]])
        sz = np.diag([0.5, -0.5])
        energies, vectors = np.linalg.eigh(j * sum(np.kron(s, s) for s in (sx, sy, sz)))
        operator = vectors.conj().T @ np.kron(sz, np.eye(2)) @ vectors
        for temperature in (0.17, 0.63):
            q = np.exp(-j / temperature)
            expected = (q / 2 + temperature * (1 - q) / (2 * j)) / (1 + 3 * q)
            compare(lehmann_chit(energies, operator, temperature), expected, "exchange-pair T*chi")
            self.assertGreater(abs(expected - 0.25), 1e-3)

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

    def test_nonconserved_spin_against_many_body_lehmann(self):
        _, h = inspect_fixture(FIXTURE)
        for shell, sites in ((h[:4, :4], [0, 2]), (h, [0, 2, 5])):
            central = shell[np.ix_(sites, sites)]
            solution = solve(AndersonModel(central[0, 0], 0, central[0, 1],
                                          tuple(np.diag(central)[1:]), tuple(np.diag(central, 1)[1:])))
            states = np.arange(len(solution.energies))
            sz = ((states & 1) - ((states >> 1) & 1)) / 2
            operator = solution.vectors.conj().T @ (sz[:, None] * solution.vectors)
            differences = solution.energies[:, None] - solution.energies[None, :]
            self.assertGreater(np.linalg.norm(differences * operator), 0.1)
            for temperature in (0.2, 0.7, 0.5 / np.log(2), np.sqrt(2) * 0.5 / np.log(2)):
                expected = lehmann_chit(solution.energies, operator, temperature)
                compare(-temperature * spin_matsubara(shell, temperature)[0], expected, "Wick vs Lehmann T*chi")
                p = np.exp(-(solution.energies - solution.energies.min()) / temperature)
                equal_time = np.dot(p / p.sum(), np.sum(np.abs(operator) ** 2, axis=1))
                self.assertGreater(abs(expected - equal_time), 1e-4)

    def test_chit_curve_contract(self):
        _, h = inspect_fixture(FIXTURE)
        temperatures = np.array([0.5 / np.log(2), np.sqrt(2) * 0.5 / np.log(2)])
        values = [-t * spin_matsubara(shell, t)[0].real for shell, t in zip((h, h[:4, :4]), temperatures)]
        rows = np.column_stack((temperatures, values))
        self.assertEqual(check_chit(rows, h), 0)
        wrong_values = [np.array(values) * factor for factor in (-1, 3, 4)]
        wrong_values.extend((np.array(values) / temperatures, values[::-1],
                             [-0.2 * spin_matsubara(shell, 0.2)[0].real for shell in (h, h[:4, :4])]))
        for values in wrong_values:
            with self.assertRaises(ValueError):
                check_chit(np.column_stack((temperatures, values)), h)
        for bad in ([], rows[:1], np.vstack((rows, rows)), rows[:, :1], rows[::-1],
                    [[0, 0], [1, 0]], [[1, 0], [1, 0]], [[np.nan, 0], [1, 0]], [[1, np.inf], [2, 0]]):
            with self.assertRaises(ValueError):
                check_chit(bad, h)

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
