"""Independent finite-star physics; no solver inputs or candidate coefficients."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from ed_siam import AndersonModel


FIXTURE = Path(__file__).resolve().parent / "fixtures" / "asymmetric_star"


def load_star(fixture=FIXTURE):
    star = json.loads((Path(fixture) / "star.json").read_text())
    if set(star) != {"energies", "couplings", "epsilon_d", "U", "B"}:
        raise ValueError("unexpected physical star fields")
    e, v = np.asarray(star["energies"]), np.asarray(star["couplings"])
    if (e.shape != (4,) or v.shape != (4,) or not np.all(np.isfinite(e))
            or not np.all(np.isfinite(v)) or not np.all(v > 0)
            or not np.all(e[::2] > 0) or not np.all(e[1::2] < 0)
            or len(set(e)) != 4 or not np.isfinite(star["epsilon_d"])
            or star["U"] != 0 or star["B"] != 0):
        raise ValueError("expected four distinct paired poles, positive couplings and U=B=0")
    return star


def jacobi(star):
    """Orthogonalize the fixed physical multiplication operator, with two full passes."""
    e, v = np.asarray(star["energies"]), np.asarray(star["couplings"])
    q = v / np.linalg.norm(v)
    basis, diagonal, hopping = [], [], []
    for index in range(len(e)):
        basis.append(q)
        residual = e * q
        diagonal.append(float(q @ residual))
        for _ in range(2):
            for vector in basis:
                residual -= (vector @ residual) * vector
        if index + 1 < len(e):
            hopping.append(float(np.linalg.norm(residual)))
            q = residual / hopping[-1]
    return np.asarray(diagonal), np.asarray(hopping)


def reference_model(star, bath_sites):
    if bath_sites not in (2, 3, 4):
        raise ValueError("reference prefix must contain 2, 3, or 4 bath sites")
    zeta, t = jacobi(star)
    return AndersonModel(star["epsilon_d"], 0, np.linalg.norm(star["couplings"]),
                         tuple(zeta[:bath_sites]), tuple(t[:bath_sites - 1]))


def direct_delta(star, frequencies):
    s = np.asarray(frequencies, dtype=complex)
    return np.sum(np.asarray(star["couplings"]) ** 2 / (s[..., None] - star["energies"]), axis=-1)


def reference_greens(star, frequencies, bath_sites):
    """Full star uses its poles directly; a strict prefix is a DIFFERENT system."""
    s = np.asarray(frequencies, dtype=complex)
    if bath_sites == 4:
        delta = direct_delta(star, s)
    else:
        model = reference_model(star, bath_sites)
        bath = np.diag(model.zeta) + np.diag(model.t, 1) + np.diag(model.t, -1)
        poles, vectors = np.linalg.eigh(bath)
        delta = np.sum((model.V * vectors[0]) ** 2 / (s[..., None] - poles), axis=-1)
    return 1 / (s - star["epsilon_d"] - delta)


def star_tables(star, D):
    """Tool/T units: positive magnitudes E/D, normalized amplitudes, physical theta."""
    if not np.isfinite(D) or D <= 0:
        raise ValueError("D must be finite and positive")
    e, v = np.asarray(star["energies"]), np.asarray(star["couplings"])
    amplitude = v / np.linalg.norm(v)
    return {"de_pos.dat": e[::2] / D, "de_neg.dat": -e[1::2] / D,
            "du_pos.dat": amplitude[::2], "du_neg.dat": amplitude[1::2],
            "theta.dat": np.asarray([np.pi * (v @ v)])}
