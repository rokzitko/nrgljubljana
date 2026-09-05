"""Independent dense ED for a real, spin-conserving Anderson chain.

Spin orbitals are ordered d_up, d_down, f0_up, f0_down, f1_up, ... .
The Hamiltonian is

    epsilon_d * n_d + U * n_d_up * n_d_down
    + B / 2 * (n_d_up - n_d_down) + sum_j zeta_j * (n_j - 1)
    + V * sum_spin (d^dagger f0 + h.c.)
    + sum_j,spin t_j * (f_j^dagger f_{j+1} + h.c.).

All particle-number sectors enter the grand-canonical ensemble, with k_B = 1
and chemical potential zero. Temperatures must be finite and strictly positive.
This module uses only the Python standard library and NumPy.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import exp, expm1, isfinite, log, pi, sqrt
from numbers import Integral, Real

import numpy as np


MAX_SPATIAL_ORBITALS = 6


def _finite_real(value: Real, name: str) -> float:
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, Real):
        raise ValueError(f"{name} must be a finite real number")
    try:
        result = float(value)
    except (OverflowError, ValueError) as error:
        raise ValueError(f"{name} must be a finite real number") from error
    if not isfinite(result):
        raise ValueError(f"{name} must be a finite real number")
    return result


@dataclass(frozen=True)
class AndersonModel:
    epsilon_d: float
    U: float
    V: float
    zeta: tuple[float, ...]
    t: tuple[float, ...]
    B: float = 0.0

    def __post_init__(self) -> None:
        for name in ("epsilon_d", "U", "V", "B"):
            object.__setattr__(self, name, _finite_real(getattr(self, name), name))
        for name in ("zeta", "t"):
            try:
                values = tuple(getattr(self, name))
            except TypeError as error:
                raise ValueError(f"{name} must be a sequence of finite real numbers") from error
            object.__setattr__(
                self, name, tuple(_finite_real(value, name) for value in values)
            )
        if not self.zeta:
            raise ValueError("at least one bath site is required")
        if len(self.zeta) + 1 > MAX_SPATIAL_ORBITALS:
            raise ValueError(
                f"dense ED is limited to {MAX_SPATIAL_ORBITALS} spatial orbitals "
                "including the impurity"
            )
        if len(self.t) != len(self.zeta) - 1:
            raise ValueError("len(t) must equal len(zeta) - 1")


def _create(state: int, orbital: int) -> tuple[int, int] | None:
    bit = 1 << orbital
    if state & bit:
        return None
    sign = 1 - 2 * ((state & (bit - 1)).bit_count() % 2)
    return state | bit, sign


def _hop(state: int, target: int, source: int) -> tuple[int, int] | None:
    bit = 1 << source
    if not state & bit:
        return None
    sign = 1 - 2 * ((state & (bit - 1)).bit_count() % 2)
    created = _create(state ^ bit, target)
    if created is None:
        return None
    return created[0], sign * created[1]


def creation_operator(spatial_orbitals: int, spin_orbital: int) -> np.ndarray:
    """Return c^dagger in the full bit basis (row = final, column = initial).

    There are two spin orbitals per spatial orbital. Bit zero is d_up, bit
    one is d_down, and increasing integers label the Fock basis states.
    """
    if (
        isinstance(spatial_orbitals, (bool, np.bool_))
        or not isinstance(spatial_orbitals, Integral)
        or not 1 <= spatial_orbitals <= MAX_SPATIAL_ORBITALS
    ):
        raise ValueError(f"spatial_orbitals must be an integer in [1, {MAX_SPATIAL_ORBITALS}]")
    if (
        isinstance(spin_orbital, (bool, np.bool_))
        or not isinstance(spin_orbital, Integral)
        or not 0 <= spin_orbital < 2 * spatial_orbitals
    ):
        raise ValueError("spin_orbital must index one of the spin orbitals")
    dimension = 1 << (2 * int(spatial_orbitals))
    operator = np.zeros((dimension, dimension))
    for state in range(dimension):
        created = _create(state, int(spin_orbital))
        if created is not None:
            final, sign = created
            operator[final, state] = sign
    return operator


@dataclass
class AndersonSolution:
    """Complete eigenbasis, grouped by (N_up, N_down), not globally sorted.

    ``vectors[:, i]`` is the full-bit-basis eigenvector of ``energies[i]``
    in ``sectors[i]``. The ground energy is ``energies.min()``.
    """

    model: AndersonModel
    energies: np.ndarray
    sectors: list[tuple[int, int]]
    vectors: np.ndarray
    _bases: dict[tuple[int, int], np.ndarray] = field(repr=False)
    _sector_slices: dict[tuple[int, int], slice] = field(repr=False)
    _observables: np.ndarray = field(repr=False)
    _transitions: list[list[tuple[slice, slice, np.ndarray]]] | None = field(
        default=None, init=False, repr=False
    )

    @property
    def sector_energies(self) -> dict[tuple[int, int], np.ndarray]:
        return {sector: self.energies[part] for sector, part in self._sector_slices.items()}

    def _boltzmann(self, T: float) -> tuple[np.ndarray, np.ndarray, float]:
        T = _finite_real(T, "T")
        if T <= 0:
            raise ValueError("T must be strictly positive")
        with np.errstate(over="ignore", under="ignore"):
            scaled = (self.energies - self.energies.min()) / T
            weights = np.exp(-scaled)
        partition = weights.sum()
        return weights / partition, scaled, float(np.log(partition))

    def thermodynamics(self, T: float) -> dict[str, float]:
        """Return energy, free energy, entropy, and heat capacity per chain."""
        T = _finite_real(T, "T")
        probabilities, scaled, log_partition = self._boltzmann(T)
        # Ignore only floating-point underflow, not a physical thermal cutoff.
        # Scaled excitation moments avoid E-F cancellation and division by T^2.
        active = probabilities > 0
        x = scaled[active]
        p = probabilities[active]
        mean = float(p @ x)
        return {
            "E": float(probabilities @ self.energies),
            "F": float(self.energies.min() - T * log_partition),
            "S": log_partition + mean,
            "C": float(p @ ((x - mean) ** 2)),
        }

    def expectations(self, T: float) -> dict[str, float]:
        """Return identity, impurity occupancy, and impurity double occupancy."""
        probabilities, _, _ = self._boltzmann(T)
        occupancy, double_occupancy = self._observables @ probabilities
        return {
            "I": float(probabilities.sum()),
            "n_d": float(occupancy),
            "n_d_ud": float(double_occupancy),
        }

    def _creation_blocks(self) -> list[list[tuple[slice, slice, np.ndarray]]]:
        if self._transitions is None:
            transitions: list[list[tuple[slice, slice, np.ndarray]]] = [[], []]
            for spin in (0, 1):
                for sector, source_basis in self._bases.items():
                    target = list(sector)
                    target[spin] += 1
                    target_sector = tuple(target)
                    if target_sector not in self._bases:
                        continue
                    target_basis = self._bases[target_sector]
                    target_rows = {int(state): row for row, state in enumerate(target_basis)}
                    source_slice = self._sector_slices[sector]
                    target_slice = self._sector_slices[target_sector]
                    source_vectors = self.vectors[source_basis, source_slice]
                    target_vectors = self.vectors[target_basis, target_slice]
                    applied = np.zeros((len(target_basis), len(source_basis)))
                    for row, state in enumerate(source_basis):
                        created = _create(int(state), spin)
                        if created is not None:
                            final, sign = created
                            applied[target_rows[final], :] = sign * source_vectors[row, :]
                    matrix_elements = target_vectors.T @ applied
                    transitions[spin].append(
                        (source_slice, target_slice, matrix_elements ** 2)
                    )
            self._transitions = transitions
        return self._transitions

    def greens(self, T: float, frequencies: np.ndarray) -> np.ndarray:
        """Return G_up(z), G_down(z) from the full finite-temperature Lehmann sum.

        Frequencies are a one-dimensional array of finite real or complex z,
        away from poles. No broadening or factor of i is added implicitly.
        """
        frequencies = np.asarray(frequencies)
        if frequencies.ndim != 1 or frequencies.dtype.kind not in "iufc":
            raise ValueError("frequencies must be a one-dimensional numeric array")
        frequencies = frequencies.astype(np.complex128)
        if not np.all(np.isfinite(frequencies)):
            raise ValueError("frequencies must be finite")
        probabilities, _, _ = self._boltzmann(T)
        result = np.zeros((2, len(frequencies)), dtype=np.complex128)
        if not len(frequencies):
            return result
        for spin, blocks in enumerate(self._creation_blocks()):
            for source, target, matrix_elements_squared in blocks:
                differences = self.energies[target, None] - self.energies[None, source]
                weights = matrix_elements_squared * (
                    probabilities[target, None] + probabilities[None, source]
                )
                nonzero = weights != 0
                differences = differences[nonzero]
                weights = weights[nonzero]
                # One frequency at a time bounds memory independently of the grid size.
                for index, frequency in enumerate(frequencies):
                    result[spin, index] += np.sum(weights / (frequency - differences))
        return result


def solve(model: AndersonModel) -> AndersonSolution:
    """Enumerate bitstrings, construct sector Hamiltonians, and retain all states."""
    if not isinstance(model, AndersonModel):
        raise TypeError("model must be an AndersonModel")
    spatial_orbitals = len(model.zeta) + 1
    dimension = 1 << (2 * spatial_orbitals)
    up_mask = sum(1 << (2 * site) for site in range(spatial_orbitals))
    down_mask = up_mask << 1
    grouped = {
        (up, down): []
        for up in range(spatial_orbitals + 1)
        for down in range(spatial_orbitals + 1)
    }
    for state in range(dimension):
        sector = ((state & up_mask).bit_count(), (state & down_mask).bit_count())
        grouped[sector].append(state)

    energies = np.empty(dimension)
    vectors = np.zeros((dimension, dimension))
    observables = np.empty((2, dimension))
    sectors = []
    bases = {}
    sector_slices = {}
    hoppings = [
        (2 * site + spin, 2 * (site + 1) + spin, amplitude)
        for site, amplitude in enumerate((model.V, *model.t))
        for spin in (0, 1)
        if amplitude != 0
    ]
    start = 0
    for sector, states in grouped.items():
        basis = np.asarray(states, dtype=np.int64)
        rows = {state: row for row, state in enumerate(states)}
        up = basis & 1
        down = (basis >> 1) & 1
        diagonal = (
            model.epsilon_d * (up + down)
            + model.U * up * down
            + model.B / 2 * (up - down)
        )
        for site, onsite in enumerate(model.zeta, start=1):
            bath_up = (basis >> (2 * site)) & 1
            bath_down = (basis >> (2 * site + 1)) & 1
            diagonal += onsite * (bath_up + bath_down - 1)
        hamiltonian = np.diag(diagonal)
        for column, state in enumerate(states):
            for left, right, amplitude in hoppings:
                for target, source in ((left, right), (right, left)):
                    hopped = _hop(state, target, source)
                    if hopped is not None:
                        final, sign = hopped
                        hamiltonian[rows[final], column] += amplitude * sign
        eigenvalues, eigenvectors = np.linalg.eigh(hamiltonian)
        part = slice(start, start + len(states))
        energies[part] = eigenvalues
        vectors[basis, part] = eigenvectors
        observables[:, part] = np.asarray([up + down, up * down]) @ (eigenvectors ** 2)
        sectors.extend([sector] * len(states))
        bases[sector] = basis
        sector_slices[sector] = part
        start = part.stop
    return AndersonSolution(model, energies, sectors, vectors, bases, sector_slices, observables)


def wilson_hopping(n: int, Lambda: float, D: float = 1.0) -> float:
    """Physical z=1 flat-band hopping t_n, without an iteration-dependent rescaling."""
    if isinstance(n, (bool, np.bool_)) or not isinstance(n, Integral) or n < 0:
        raise ValueError("n must be a nonnegative integer")
    Lambda = _finite_real(Lambda, "Lambda")
    D = _finite_real(D, "D")
    if Lambda <= 1:
        raise ValueError("Lambda must be greater than one")
    if D <= 0:
        raise ValueError("D must be positive")
    n = int(n)
    alpha = log(Lambda)
    # expm1 preserves 1-q^k as Lambda approaches one.
    numerator = -expm1(-(n + 1) * alpha)
    denominator = sqrt((-expm1(-(2 * n + 1) * alpha)) * (-expm1(-(2 * n + 3) * alpha)))
    return D * ((-expm1(-alpha) / alpha) * exp(-n * alpha / 2) * numerator / denominator)


def flat_band_model(parameters: dict) -> AndersonModel:
    """Build the finite physical z=1 chain from explicit flat-band parameters."""
    required = ("epsilon_d", "U", "Gamma", "D", "Lambda", "z", "bath_sites")
    missing = [name for name in required if name not in parameters]
    if missing:
        raise ValueError(f"missing flat-band parameters: {', '.join(missing)}")
    z = _finite_real(parameters["z"], "z")
    if z != 1:
        raise ValueError("only z=1 is supported")
    Gamma = _finite_real(parameters["Gamma"], "Gamma")
    if Gamma < 0:
        raise ValueError("Gamma must be nonnegative")
    D = _finite_real(parameters["D"], "D")
    Lambda = _finite_real(parameters["Lambda"], "Lambda")
    # Validate D and Lambda even for a single bath site, which has no t entries.
    wilson_hopping(0, Lambda, D)
    bath_sites = parameters["bath_sites"]
    if (
        isinstance(bath_sites, (bool, np.bool_))
        or not isinstance(bath_sites, Integral)
        or not 1 <= bath_sites < MAX_SPATIAL_ORBITALS
    ):
        raise ValueError(f"bath_sites must be an integer in [1, {MAX_SPATIAL_ORBITALS - 1}]")
    bath_sites = int(bath_sites)
    return AndersonModel(
        epsilon_d=parameters["epsilon_d"],
        U=parameters["U"],
        V=sqrt(2 / pi) * sqrt(D) * sqrt(Gamma),
        zeta=(0.0,) * bath_sites,
        t=tuple(wilson_hopping(n, Lambda, D) for n in range(bath_sites - 1)),
        B=parameters.get("B", 0.0),
    )
