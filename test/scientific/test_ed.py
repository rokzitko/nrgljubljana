"""NumPy/stdlib self-checks for the independent Anderson-chain ED reference."""

from collections import Counter
from dataclasses import replace
from math import comb, log, pi, sqrt
import unittest

import numpy as np

from ed_siam import (
    AndersonModel,
    MAX_SPATIAL_ORBITALS,
    creation_operator,
    flat_band_model,
    solve,
    wilson_hopping,
)


class AndersonEDTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.model = AndersonModel(-0.41, 0.93, -0.27, (0.14, -0.22), (0.31,), B=0.19)
        cls.solution = solve(cls.model)

    @staticmethod
    def single_particle(model, spin):
        diagonal = [model.epsilon_d + (1 if spin == 0 else -1) * model.B / 2, *model.zeta]
        hopping = np.asarray([model.V, *model.t])
        return np.diag(diagonal) + np.diag(hopping, 1) + np.diag(hopping, -1)

    def test_creation_operators_and_fermion_car(self):
        operators = [creation_operator(3, orbital) for orbital in range(6)]
        identity = np.eye(64)
        for i, left in enumerate(operators):
            for j, right in enumerate(operators):
                with self.subTest(i=i, j=j):
                    np.testing.assert_array_equal(left @ right + right @ left, np.zeros((64, 64)))
                    np.testing.assert_array_equal(
                        left.T @ right + right @ left.T,
                        identity if i == j else np.zeros((64, 64)),
                    )
        self.assertEqual(operators[0][1, 0], 1)
        self.assertEqual(operators[1][3, 1], -1)
        self.assertEqual(operators[2][5, 1], -1)
        self.assertEqual(operators[2][7, 3], 1)
        self.assertEqual(operators[0][0, 1], 0)

    def test_hamiltonian_eigenpairs_and_sector_counts(self):
        model, solution = self.model, self.solution
        spatial_orbitals = len(model.zeta) + 1
        dimension = 4 ** spatial_orbitals
        identity = np.eye(dimension)
        creators = [creation_operator(spatial_orbitals, i) for i in range(2 * spatial_orbitals)]
        numbers = [operator @ operator.T for operator in creators]
        # An independent dense-operator construction checks the bit-built Hamiltonian.
        hamiltonian = (
            model.epsilon_d * (numbers[0] + numbers[1])
            + model.U * (numbers[0] @ numbers[1])
            + model.B / 2 * (numbers[0] - numbers[1])
        )
        for site, onsite in enumerate(model.zeta, start=1):
            hamiltonian += onsite * (numbers[2 * site] + numbers[2 * site + 1] - identity)
        for site, hopping in enumerate((model.V, *model.t)):
            for spin in (0, 1):
                left, right = creators[2 * site + spin], creators[2 * (site + 1) + spin]
                hamiltonian += hopping * (left @ right.T + right @ left.T)
        np.testing.assert_array_equal(hamiltonian, hamiltonian.T)
        np.testing.assert_allclose(solution.vectors.T @ solution.vectors, identity, atol=2e-14)
        np.testing.assert_allclose(
            hamiltonian @ solution.vectors,
            solution.vectors * solution.energies,
            atol=2e-14,
        )
        np.testing.assert_allclose(
            (solution.vectors * solution.energies) @ solution.vectors.T,
            hamiltonian,
            atol=2e-14,
        )
        expected_counts = {
            (up, down): comb(spatial_orbitals, up) * comb(spatial_orbitals, down)
            for up in range(spatial_orbitals + 1)
            for down in range(spatial_orbitals + 1)
        }
        self.assertEqual(Counter(solution.sectors), expected_counts)
        self.assertEqual(set(solution.sector_energies), set(expected_counts))
        self.assertEqual(solution.energies.shape, (dimension,))
        self.assertEqual(solution.vectors.shape, (dimension, dimension))
        self.assertIs(solution.model, model)
        bit_sectors = [
            tuple(sum((state >> (2 * site + spin)) & 1 for site in range(spatial_orbitals))
                  for spin in (0, 1))
            for state in range(dimension)
        ]
        for sector, count in expected_counts.items():
            rows = [state for state, bit_sector in enumerate(bit_sectors) if bit_sector == sector]
            columns = [i for i, eigen_sector in enumerate(solution.sectors) if eigen_sector == sector]
            other_rows = [state for state, bit_sector in enumerate(bit_sectors) if bit_sector != sector]
            self.assertEqual(len(solution.sector_energies[sector]), count)
            np.testing.assert_array_equal(solution.vectors[np.ix_(other_rows, columns)], 0)
            np.testing.assert_array_equal(solution.sector_energies[sector], solution.energies[columns])
            np.testing.assert_allclose(
                solution.sector_energies[sector],
                np.linalg.eigvalsh(hamiltonian[np.ix_(rows, rows)]),
                atol=2e-14,
            )

    def test_decoupled_atomic_limit(self):
        model = AndersonModel(-0.37, 0.91, 0.0, (0.0,), (), B=0.18)
        solution = solve(model)
        atomic_energies = np.array([0, model.epsilon_d + model.B / 2,
                                    model.epsilon_d - model.B / 2, 2 * model.epsilon_d + model.U])
        np.testing.assert_allclose(np.sort(solution.energies), np.sort(np.repeat(atomic_energies, 4)))
        self.assertLess(solution.energies.min(), solution.energies[0])
        frequencies = np.array([0.17j, 0.4 + 0.3j, -0.8 + 0.21j])
        for T in (0.07, 0.6, 4.0):
            with self.subTest(T=T):
                weights = np.exp(-(atomic_energies - atomic_energies.min()) / T)
                p = weights / weights.sum()
                expected_energy = float(p @ atomic_energies)
                expected_free = atomic_energies.min() - T * log(4 * weights.sum())
                expected_entropy = -float(p @ np.log(p)) + log(4)
                expected_heat = float(p @ ((atomic_energies - expected_energy) ** 2)) / T ** 2
                thermo = solution.thermodynamics(T)
                self.assertEqual(set(thermo), {"E", "F", "S", "C"})
                for name, expected in zip(("E", "F", "S", "C"),
                                          (expected_energy, expected_free, expected_entropy, expected_heat)):
                    self.assertAlmostEqual(thermo[name], expected, places=13)
                values = solution.expectations(T)
                self.assertEqual(set(values), {"I", "n_d", "n_d_ud"})
                self.assertAlmostEqual(values["I"], 1)
                self.assertAlmostEqual(values["n_d"], p[1] + p[2] + 2 * p[3])
                self.assertAlmostEqual(values["n_d_ud"], p[3])
                greens = solution.greens(T, frequencies)
                self.assertEqual(greens.shape, (2, len(frequencies)))
                for spin in (0, 1):
                    other_occupancy = p[2 - spin] + p[3]
                    onsite = model.epsilon_d + (1 if spin == 0 else -1) * model.B / 2
                    expected = ((1 - other_occupancy) / (frequencies - onsite)
                                + other_occupancy / (frequencies - onsite - model.U))
                    np.testing.assert_allclose(greens[spin], expected, rtol=2e-13, atol=2e-13)

    def test_noninteracting_greens_match_single_particle_resolvent(self):
        frequencies = np.array([0.17j, 1.2j, -1.1 + 0.13j, 0.7 + 0.5j, 10.0])
        for B in (0.0, 0.31, -0.23):
            model = replace(self.model, U=0.0, B=B)
            solution = solve(model)
            for T in (0.035, 0.4, 3.0):
                actual = solution.greens(T, frequencies)
                for spin in (0, 1):
                    with self.subTest(B=B, T=T, spin=spin):
                        one_body = self.single_particle(model, spin)
                        expected = np.array([
                            np.linalg.inv(frequency * np.eye(len(one_body)) - one_body)[0, 0]
                            for frequency in frequencies
                        ])
                        np.testing.assert_allclose(actual[spin], expected, rtol=3e-13, atol=3e-13)
            np.testing.assert_allclose(
                solution.greens(0.4, frequencies.conjugate()),
                solution.greens(0.4, frequencies).conjugate(),
                atol=2e-14,
            )

    def test_finite_temperature_occupation_and_spectral_moments(self):
        for U in (0.0, self.model.U):
            solution = self.solution if U == self.model.U else solve(replace(self.model, U=U))
            model = solution.model
            energies, vectors = solution.energies, solution.vectors
            dimension = len(energies)
            states = np.arange(dimension)
            differences = energies[:, None] - energies[None, :]
            transformed = [vectors.T @ creation_operator(3, spin) @ vectors for spin in (0, 1)]
            for T in (0.03, 0.5, 8.0):
                p = np.exp(-(energies - energies.min()) / T)
                p /= p.sum()
                bit_probabilities = (vectors ** 2) @ p
                direct_occupancy = [float(bit_probabilities @ ((states >> spin) & 1)) for spin in (0, 1)]
                values = solution.expectations(T)
                self.assertAlmostEqual(values["n_d"], sum(direct_occupancy), places=13)
                self.assertAlmostEqual(
                    values["n_d_ud"], float(bit_probabilities @ ((states & 3) == 3)), places=13
                )
                exponential = np.exp(-np.abs(differences) / T)
                fermi = np.where(differences >= 0, exponential / (1 + exponential), 1 / (1 + exponential))
                omega = 1e5
                asymptotic = solution.greens(T, np.array([1j * omega]))[:, 0]
                for spin in (0, 1):
                    with self.subTest(U=U, T=T, spin=spin):
                        weights = (p[:, None] + p[None, :]) * transformed[spin] ** 2
                        self.assertAlmostEqual(float(np.sum(weights * fermi)), direct_occupancy[spin], places=13)
                        self.assertAlmostEqual(float(weights.sum()), 1, places=13)
                        first_moment = (model.epsilon_d + (1 if spin == 0 else -1) * model.B / 2
                                        + U * direct_occupancy[1 - spin])
                        self.assertAlmostEqual(float(np.sum(weights * differences)), first_moment, places=13)
                        self.assertAlmostEqual(-omega * asymptotic[spin].imag, 1, delta=2e-9)
                        self.assertAlmostEqual(-omega ** 2 * asymptotic[spin].real, first_moment, delta=2e-9)
                        if U == 0:
                            one_energies, one_vectors = np.linalg.eigh(self.single_particle(model, spin))
                            one_exp = np.exp(-np.abs(one_energies) / T)
                            one_fermi = np.where(one_energies >= 0, one_exp / (1 + one_exp), 1 / (1 + one_exp))
                            expected = float((one_vectors[0] ** 2) @ one_fermi)
                            self.assertAlmostEqual(direct_occupancy[spin], expected, places=13)
                if U == 0:
                    self.assertAlmostEqual(values["n_d_ud"], direct_occupancy[0] * direct_occupancy[1], places=13)

    def test_thermodynamic_limits_and_stability(self):
        solution = self.solution
        T = 1e9
        values = solution.thermodynamics(T)
        self.assertAlmostEqual(values["S"], log(len(solution.energies)), places=13)
        self.assertAlmostEqual(values["E"], float(solution.energies.mean()), delta=2e-9)
        self.assertGreaterEqual(values["C"], 0)
        self.assertLess(values["C"], 1e-16)
        expected = solution.expectations(T)
        self.assertAlmostEqual(expected["n_d"], 1, delta=2e-9)
        self.assertAlmostEqual(expected["n_d_ud"], 0.25, delta=2e-9)
        atomic = solve(AndersonModel(-0.37, 0.91, 0.0, (0.0,), (), B=0.18))
        with np.errstate(over="raise", invalid="raise", divide="raise"):
            for low_T in (1e-12, 1e-300, np.nextafter(0.0, 1.0)):
                thermo = atomic.thermodynamics(low_T)
                self.assertEqual(thermo["E"], float(atomic.energies.min()))
                self.assertAlmostEqual(thermo["F"], float(atomic.energies.min()), delta=2e-12)
                self.assertAlmostEqual(thermo["S"], log(4), places=14)
                self.assertEqual(thermo["C"], 0)

    def test_numpy_scalar_temperatures_use_double_precision(self):
        solution = solve(AndersonModel(0, 0, 0, (0,), ()))
        for T in (np.float16(1), np.float32(1), np.float64(1)):
            with self.subTest(T=type(T).__name__):
                self.assertEqual(solution.thermodynamics(T), solution.thermodynamics(1.0))
                self.assertAlmostEqual(solution.thermodynamics(T)["F"], -log(16), places=14)

    def test_flat_band_model_and_wilson_hopping(self):
        parameters = dict(epsilon_d=-0.3, U=0.8, Gamma=0.07, D=1.4,
                          Lambda=2.3, z=1, bath_sites=3, B=-0.11)
        model = flat_band_model(parameters)
        self.assertEqual(model.epsilon_d, parameters["epsilon_d"])
        self.assertEqual(model.U, parameters["U"])
        self.assertEqual(model.B, parameters["B"])
        self.assertAlmostEqual(model.V, sqrt(2 * parameters["D"] * parameters["Gamma"] / pi), places=15)
        self.assertEqual(model.zeta, (0.0, 0.0, 0.0))
        self.assertEqual(model.t, tuple(wilson_hopping(n, 2.3, 1.4) for n in range(2)))
        single = flat_band_model({**parameters, "bath_sites": 1, "Gamma": 0})
        self.assertEqual(single.t, ())
        self.assertEqual(single.V, 0)
        without_field = {name: value for name, value in parameters.items() if name != "B"}
        self.assertEqual(flat_band_model(without_field).B, 0)
        for Lambda in (1.3, 2.0, 4.0, 10.0):
            q = 1 / Lambda
            for n in range(9):
                expected = (0.8 * (1 - q) / log(Lambda) * q ** (n / 2) * (1 - q ** (n + 1))
                            / sqrt((1 - q ** (2 * n + 1)) * (1 - q ** (2 * n + 3))))
                self.assertAlmostEqual(wilson_hopping(n, Lambda, 0.8), expected, places=14)
                self.assertAlmostEqual(wilson_hopping(n, Lambda), expected / 0.8, places=14)
        for n in range(5):
            continuum = (n + 1) / sqrt((2 * n + 1) * (2 * n + 3))
            self.assertAlmostEqual(wilson_hopping(n, np.nextafter(1.0, 2.0)), continuum, places=14)

    def test_invalid_models_and_operator_indices(self):
        for name in ("epsilon_d", "U", "V", "B"):
            for value in (np.nan, np.inf, -np.inf, 0j, "0.1", True, None):
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    replace(self.model, **{name: value})
        for changes in (
            {"zeta": (), "t": ()},
            {"zeta": (0,) * MAX_SPATIAL_ORBITALS, "t": (0,) * (MAX_SPATIAL_ORBITALS - 1)},
            {"t": ()}, {"t": (0, 0)}, {"zeta": None}, {"t": None},
            {"zeta": (0, np.nan)}, {"zeta": (0, 0j)}, {"zeta": (0, np.inf)},
            {"t": (np.inf,)}, {"t": (0j,)}, {"t": ("0",)},
        ):
            with self.subTest(changes=changes), self.assertRaises(ValueError):
                replace(self.model, **changes)
        largest = AndersonModel(0, 0, 0, (0,) * (MAX_SPATIAL_ORBITALS - 1),
                                (0,) * (MAX_SPATIAL_ORBITALS - 2))
        self.assertEqual(len(largest.zeta) + 1, MAX_SPATIAL_ORBITALS)
        normalized = AndersonModel(np.float64(-0.4), 1, 0.2, [0, 0.1], np.array([0.3]))
        self.assertEqual(normalized.zeta, (0.0, 0.1))
        self.assertEqual(normalized.t, (0.3,))
        for spatial, orbital in ((0, 0), (7, 0), (1.5, 0), (True, 0),
                                 (2, -1), (2, 4), (2, 0.5), (2, True)):
            with self.subTest(spatial=spatial, orbital=orbital), self.assertRaises(ValueError):
                creation_operator(spatial, orbital)

    def test_invalid_temperatures_and_frequencies(self):
        for T in (0, -0.1, np.nan, np.inf, 0j, "1", True):
            for method in (self.solution.thermodynamics, self.solution.expectations):
                with self.subTest(T=T, method=method.__name__), self.assertRaises(ValueError):
                    method(T)
            with self.subTest(T=T, method="greens"), self.assertRaises(ValueError):
                self.solution.greens(T, np.array([1j]))
        for frequencies in (np.array(1j), np.array([[1j]]), np.array([np.nan]),
                            np.array([complex(0, np.inf)]), np.array(["1j"]), np.array([True])):
            with self.subTest(frequencies=frequencies), self.assertRaises(ValueError):
                self.solution.greens(1, frequencies)
        self.assertEqual(self.solution.greens(1, np.array([], dtype=complex)).shape, (2, 0))

    def test_invalid_flat_band_parameters(self):
        parameters = dict(epsilon_d=-0.3, U=0.8, Gamma=0.07, D=1.4, Lambda=2.3, z=1, bath_sites=1)
        for name in parameters:
            with self.subTest(missing=name), self.assertRaises(ValueError):
                flat_band_model({key: value for key, value in parameters.items() if key != name})
        for name, values in {
            "Gamma": (-0.1, np.inf, np.nan, 0j),
            "D": (0, -1, np.inf, np.nan, 1j),
            "Lambda": (0, 1, -2, np.inf, np.nan, 2j),
            "z": (0, 0.5, 2, np.nan, np.inf, 1j, True),
            "bath_sites": (0, -1, MAX_SPATIAL_ORBITALS, 1.5, True),
            "epsilon_d": (np.inf,), "U": (np.nan,), "B": (0j,),
        }.items():
            for value in values:
                with self.subTest(name=name, value=value), self.assertRaises(ValueError):
                    flat_band_model({**parameters, name: value})
        for n in (-1, 0.5, True):
            with self.subTest(n=n), self.assertRaises(ValueError):
                wilson_hopping(n, 2)
        for Lambda in (1, 0, -1, np.nan, np.inf, 2j):
            with self.subTest(Lambda=Lambda), self.assertRaises(ValueError):
                wilson_hopping(0, Lambda)
        for D in (0, -1, np.nan, np.inf, 1j):
            with self.subTest(D=D), self.assertRaises(ValueError):
                wilson_hopping(0, 2, D)


if __name__ == "__main__":
    unittest.main()
