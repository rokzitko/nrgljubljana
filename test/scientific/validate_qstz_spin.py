"""Compare untruncated QSTZ spin Matsubara and T*chi data with free-fermion Wick sums."""

import argparse
from collections import Counter
import configparser
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np


FIXTURE = Path(__file__).resolve().parent / "fixtures/qstz_spin"
MATS = 4


def require(condition, message):
    if not condition:
        raise ValueError(message)


def compare(actual, expected, label, atol=1e-11, rtol=1e-10):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape and actual.size > 0, f"{label}: wrong shape or empty data")
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), f"{label}: nonfinite values")
    error = np.abs(actual - expected)
    require(np.all(error <= atol + rtol * np.maximum(np.abs(actual), np.abs(expected))),
            f"{label}: max error={error.max():.6g}; actual={actual}; expected={expected}")
    return float(error.max())


def spin_matsubara(h, temperature):
    """One Cartesian component, -<T_tau Sz(tau) Sz(0)>, with Sz=(n_up-n_down)/2."""
    require(np.isfinite(temperature) and temperature > 0, "temperature must be positive and finite")
    energies, vectors = np.linalg.eigh(h)
    weights = np.abs(vectors[0]) ** 2
    f = np.exp(-np.logaddexp(0, energies / temperature))
    omega = 2 * np.pi * temperature * np.arange(MATS)
    answer = np.zeros(MATS, dtype=complex)
    for i, ei in enumerate(energies):
        for j, ej in enumerate(energies):
            weight = weights[i] * weights[j] / 2
            if abs(ei - ej) <= 1e-14:
                # Includes off-diagonal pairs in a degenerate eigenspace.
                answer[0] -= weight * f[i] * (1 - f[i]) / temperature
            else:
                answer += weight * (f[i] - f[j]) / (1j * omega + ei - ej)
    return answer


def dimensions(charges):
    """Count (Q,2S+1,Tz) multiplets from occupation states, not NRG tables."""
    sectors = Counter({(0, 0, 0): 1})
    for charge in charges:
        following = Counter()
        for (q, sz2, tz), count in sectors.items():
            for up, down in ((0, 0), (1, 0), (0, 1), (1, 1)):
                following[q + up + down - 1, sz2 + up - down, tz + charge * (up + down)] += count
        sectors = following
    return {(q, sz2 + 1, tz): count - sectors[q, sz2 + 2, tz]
            for (q, sz2, tz), count in sectors.items()
            if sz2 >= 0 and count > sectors[q, sz2 + 2, tz]}


def inspect_fixture(fixture):
    provenance = json.loads((fixture / "provenance.json").read_text())
    require(set(provenance["files"]) == {"param", "model.m", "data"}, "incomplete fixture manifest")
    for name, digest in provenance["files"].items():
        require(hashlib.sha256((fixture / name).read_bytes()).hexdigest() == digest,
                f"fixture checksum mismatch: {name}")
    parameters = (fixture / "param").read_text()
    config = configparser.ConfigParser(interpolation=None)
    config.read_string(parameters)
    eps, v = config.getfloat("extra", "eps"), config.getfloat("extra", "V")
    text = (fixture / "data").read_text()
    for header in (r"^#!9$", r"^# symtype\s+QSTZ$", r"^# SCALE\s+1$"):
        require(re.search(header, text, re.M), "unexpected fixture format, symmetry or scaling")
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.startswith("#")]
    channels, nmax, nsectors = map(int, lines[0].split())
    require((channels, nmax) == (3, 1), "expected exactly one three-site extension")
    counts, expanded, position = {}, [], 1
    for _ in range(nsectors):
        qn = tuple(map(int, lines[position].split()))
        count = int(lines[position + 1])
        energies = list(map(float, lines[position + 2].split()))
        require(len(qn) == 3 and qn not in counts and count == len(energies), "invalid seed sector")
        counts[qn] = count
        expanded.extend(energies * qn[1])
        position += 3
    require(counts == dimensions((0, 1, 0, -1)), "seed must contain all 256 physical states")
    require(lines[position] == "f 0 0" and lines.count("tsigma_d") == 1, "missing boundary/spin operator")
    require(lines.count("z") == lines.count("e") == 1, "missing/repeated chain or ground-energy block")
    trailer = lines[lines.index("z") + 1:]
    require(len(trailer) == 18 and trailer[::3] == ["1"] * 6, "wrong chain coefficient dimensions")
    coefficients = np.array([[float(x) for x in trailer[i + 1:i + 3]] for i in range(0, 18, 3)])
    # Flat Z discretization at Lambda=2,z=1,D=1, including the unused final hopping.
    hopping = [0.5 / np.log(2) * 2 ** (-n / 2) * (1 - 2 ** (-(n + 1))) /
               np.sqrt((1 - 2 ** (-(2 * n + 1))) * (1 - 2 ** (-(2 * n + 3)))) for n in (0, 1)]
    compare(coefficients, np.array([hopping] * 3 + [[0, 0]] * 3), "chain", atol=5e-15, rtol=5e-15)
    # Order: impurity, (+,0,-) seed sites, (+,0,-) extension sites.
    h = np.zeros((7, 7))
    h[0, 0], h[0, 2], h[2, 0] = eps, v, v
    for channel in range(3):
        h[1 + channel, 1 + channel], h[4 + channel, 4 + channel] = coefficients[3 + channel]
        h[1 + channel, 4 + channel] = h[4 + channel, 1 + channel] = coefficients[channel, 0]
    seed_energies = np.array([0.0])
    for energy in np.repeat(np.linalg.eigvalsh(h[:4, :4]), 2):
        seed_energies = np.concatenate((seed_energies, seed_energies + energy))
    ground = float(lines[lines.index("e") + 1])
    compare(np.sort(expanded) + ground, np.sort(seed_energies), "physical seed spectrum")
    return parameters, h


def check_retention(log):
    require("Number of states (multiplicity taken into account): 256" in log, "wrong initial space")
    rows = re.findall(r"nrkept=(\d+) nrkeptmult=(\d+) nrall=(\d+) nrallmult=(\d+)", log)
    counts = dimensions((0, 1, 0, -1, 1, 0, -1))
    multiplets = sum(counts.values())
    require(sum(ss * count for (_, ss, _), count in counts.items()) == 16384, "wrong reference dimension")
    require(rows and all(tuple(map(int, row)) == (multiplets, 16384, multiplets, 16384) for row in rows),
            f"truncation or incomplete diagonalization: {rows}")
    require(log.count("[NRG] Iteration 1/1") == log.count("[DM] Iteration 1/1") == 1,
            "expected exactly one NRG and one DM extension")


def check_chit(rows, h):
    """Check the complete CHIT curve: final 7-site shell, then 4-site seed, sorted by T."""
    rows = np.asarray(rows)
    require(rows.shape == (2, 2) and np.all(np.isfinite(rows)), "expected two finite CHIT (T,T*chi) rows")
    require(0 < rows[0, 0] < rows[1, 0], "CHIT shell temperatures must be positive and increasing")
    # CHIT reports T*chi at chitp*scale, not at the fixed FDM temperature.
    # The bosonic Green function has the opposite sign to the susceptibility.
    expected = [-temperature * spin_matsubara(shell, temperature)[0].real
                for shell, temperature in zip((h, h[:4, :4]), rows[:, 0])]
    return compare(rows[:, 1], expected, "spin T*chi")


def run(nrg, fixture, temperature, work_root):
    parameters, h = inspect_fixture(fixture)
    expected = spin_matsubara(h, temperature)
    parameters, count = re.subn(r"(?m)^T=.*$", f"T={temperature}\nspecchit=sigma_d-sigma_d", parameters)
    require(count == 1, "expected exactly one temperature parameter")
    work_root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="qstz-spin-", dir=work_root.resolve()))
    scratch = directory / "scratch"
    scratch.mkdir()
    shutil.copyfile(fixture / "data", directory / "data")
    (directory / "param").write_text(parameters)
    environment = os.environ.copy()
    environment.update({name: "1" for name in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                                              "BLIS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")})
    environment.update(OMP_DYNAMIC="FALSE", MKL_DYNAMIC="FALSE", LC_ALL="C",
                       TMPDIR=str(scratch), TMP=str(scratch), TEMP=str(scratch))
    report = {"temperature": temperature, "nrg": str(nrg.resolve()), "numpy": np.__version__,
              "expected_real": expected.real.tolist(), "expected_imag": expected.imag.tolist()}
    print(f"QSTZ spin T={temperature}; run directory: {directory}", flush=True)
    try:
        with (directory / "log").open("w") as log:
            result = subprocess.run([str(nrg.resolve()), "-w", str(scratch)], cwd=directory, env=environment,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=120)
        require(result.returncode == 0, f"NRG exited {result.returncode}; see {directory / 'log'}")
        require((directory / "DONE").is_file(), "missing fresh DONE")
        check_retention((directory / "log").read_text())
        rows = np.loadtxt(directory / "spin_FDMmats_dens_sigma_d-sigma_d.dat", ndmin=2)
        require(rows.shape == (MATS, 3), "expected four bosonic Matsubara rows, including zero")
        compare(rows[:, 0], 2 * np.pi * temperature * np.arange(MATS), "frequencies", atol=1e-14, rtol=1e-14)
        report["max_real_error"] = compare(rows[:, 1], expected.real, "Re spin correlator")
        report["max_imag_error"] = compare(rows[:, 2], expected.imag, "Im spin correlator")
        chit = np.loadtxt(directory / "chit_CHIT_dens_sigma_d-sigma_d.dat", ndmin=2)
        report["max_chit_error"] = check_chit(chit, h)
        report["chit_curve"] = chit.tolist()
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (directory / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: four bosonic frequencies and two CHIT shells; errors "
          f"{report['max_real_error']:.3g}, {report['max_imag_error']:.3g}, {report['max_chit_error']:.3g}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nrg", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, default=FIXTURE)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        run(args.nrg, args.fixture, args.temperature, args.work_root)
    except (ValueError, OSError, KeyError, IndexError, subprocess.SubprocessError) as error:
        parser.exit(1, f"QSTZ spin validation failed: {error}\n")


if __name__ == "__main__":
    main()
