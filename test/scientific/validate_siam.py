"""Run untruncated NRG and compare it with independently constructed finite-chain ED."""

from __future__ import annotations

import argparse
from dataclasses import replace
import hashlib
import json
from math import comb
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np

from ed_siam import flat_band_model, solve, wilson_hopping


MATS = 64
ATOL = 1e-11
RTOL = 1e-10


def require(condition, message):
    if not condition:
        raise ValueError(message)


def finite_values(tokens):
    values = np.asarray([float(token) for token in tokens])
    require(np.all(np.isfinite(values)), "nonfinite numerical output")
    return values


def case_digest(case):
    return hashlib.sha256(json.dumps(case, sort_keys=True, allow_nan=False).encode()).hexdigest()


def file_digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load_case(fixture):
    case = json.loads((Path(fixture) / "model.json").read_text())
    require(set(case) == {"model", "symmetry"}, "unexpected case specification fields")
    require(case["symmetry"] in {"QS", "QSZ"}, "only QS and QSZ fixtures are supported")
    require(set(case["model"]) == {
        "epsilon_d", "U", "Gamma", "D", "Lambda", "z", "bath_sites", "B"
    }, "unexpected physical model fields")
    model = flat_band_model(case["model"])
    require(len(model.zeta) >= 2, "NRG requires at least one extension beyond d+f0")
    require(model.B == 0 or case["symmetry"] == "QSZ", "QS requires zero magnetic field")
    return case


def parameter_text(case, temperature, mode, mmax=80, precision=1000):
    require(np.isfinite(temperature) and temperature > 0, "temperature must be finite and positive")
    require(mode in {"rescaled", "absolute"}, "unknown execution mode")
    p = case["model"]
    pair = "A_d-A_d" if case["symmetry"] == "QS" else "A_d-A_d-u A_d-A_d-d"
    extra = {"eps": p["epsilon_d"], "U": p["U"], "Gamma": p["Gamma"], "B": p["B"]}
    params = {
        "model": "SIAM", "variant": "EPS" if p["B"] == 0 else "MAGFIELDEPS",
        "symtype": case["symmetry"], "band": "flat", "bandrescale": p["D"],
        "discretization": "Z", "Lambda": p["Lambda"], "z": p["z"],
        "Ninit": 0, "Nmax": p["bath_sites"] - 1, "tri": "old",
        "wilsonchain": "legacy", "mMAX": mmax, "prec": precision,
        "data_has_rescaled_energies": "false", "polarized": "false", "substeps": "false",
        "ops": "I A_d n_d n_d_ud", "specd": pair,
        "keep": 4 ** (p["bath_sites"] + 1), "keepenergy": -1, "keepmin": 0,
        "safeguard": 0, "fixeps": 0, "diag": "dsyevd", "diagratio": 1,
        "diag_mode": "serial", "mult": "blas", "strategy": "all",
        "lastall": "true", "lastalloverride": "false", "calc0": "true",
        "T": temperature, "fdmmats": "true", "fdmexpv": "true", "fdmexpvn": 0,
        "fdm_cutoff": 0, "checkrho": "true", "mats": MATS,
        "finite": "false", "finitemats": "false", "dmnrg": "false", "dmnrgmats": "false",
        "cfs": "false", "cfsgt": "false", "cfsls": "false",
        "fdm": "false", "fdmgt": "false", "fdmls": "false",
        "broaden": "false", "savebins": "false", "clip_tol_imag": 0,
        "absolute": "true" if mode == "absolute" else "false", "silent": "false", "done": "true",
        "dumpsubspaces": "true", "reportdiagonal": 4 ** (p["bath_sites"] + 1),
        "reportdiagonallast": "false", "dumpscaled": "false", "dumpabs": "true", "dumpEscale": 1,
        "prec_xy": 17, "prec_custom": 17, "prec_td": 17,
        "width_custom": 24, "width_td": 24,
    }
    return "\n\n".join(
        f"[{section}]\n" + "\n".join(f"{key}={value}" for key, value in values.items())
        for section, values in (("extra", extra), ("param", params))
    ) + "\n"


def dimensions(orbitals, symmetry):
    qsz = {
        (up + down - orbitals, up - down + 1): comb(orbitals, up) * comb(orbitals, down)
        for up in range(orbitals + 1) for down in range(orbitals + 1)
    }
    if symmetry == "QSZ":
        return qsz
    return {
        (charge, spin): count
        for charge in range(-orbitals, orbitals + 1)
        for spin in range(1, orbitals + 2)
        if (count := qsz.get((charge, spin), 0) - qsz.get((charge, spin + 2), 0)) > 0
    }


def compare(label, actual, expected, checks, atol=ATOL, rtol=RTOL):
    actual, expected = np.asarray(actual), np.asarray(expected)
    require(actual.shape == expected.shape and actual.size > 0,
            f"{label}: shape mismatch or empty comparison: {actual.shape} vs {expected.shape}")
    require(np.all(np.isfinite(actual)) and np.all(np.isfinite(expected)), f"{label}: nonfinite values")
    errors = np.abs(actual - expected)
    limits = atol + rtol * np.maximum(np.abs(actual), np.abs(expected))
    worst = np.unravel_index(np.argmax(errors / limits), errors.shape)
    require(np.all(errors <= limits),
            f"{label} at {worst}: actual={actual[worst]}, expected={expected[worst]}, "
            f"error={errors[worst]:.6g}, limit={limits[worst]:.6g}")
    checks.append({"quantity": label, "count": actual.size, "max_absolute_error": float(errors.max())})


def compare_spectrum(sectors, solution, symmetry, label, checks):
    orbitals = len(solution.model.zeta) + 1
    counts = dimensions(orbitals, symmetry)
    require(set(sectors) == set(counts), f"{label}: missing or unexpected sectors")
    expanded = {}
    for (charge, spin), energies in sectors.items():
        require(len(energies) == counts[charge, spin], f"{label}: wrong dimension for {(charge, spin)}")
        projections = [spin] if symmetry == "QSZ" else range(2 - spin, spin + 1, 2)
        for projection in projections:
            expanded.setdefault((charge, projection), []).extend(energies)
    expected = {
        (up + down - orbitals, up - down + 1): values
        for (up, down), values in solution.sector_energies.items()
    }
    require(set(expanded) == set(expected), f"{label}: incorrect physical sectors")
    ground = min(min(values) for values in expanded.values())
    compare(f"{label}/ground", ground, solution.energies.min(), checks)
    for sector, values in expected.items():
        energies = np.sort(expanded[sector])
        compare(f"{label}/sector{sector}", energies, values, checks)
        compare(f"{label}/gaps{sector}", energies - ground, values - solution.energies.min(), checks)


def inspect_data(path, case, checks):
    """Read only the real legacy seed spectrum, e block, and final chain trailer."""
    text = Path(path).read_text()
    require(re.search(r"^#!9\s*$", text, re.M), "fixture must have format #!9")
    require(re.search(r"^# symtype\s+" + case["symmetry"] + r"\s*$", text, re.M), "wrong seed symmetry")
    require(re.search(r"^# SCALE\s+1(?:\.0*)?\s*$", text, re.M), "fixture seed must be in physical units")
    require("COMPLEX" not in text, "complex fixture format is not supported")
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.startswith("#")]
    channels, nmax, nsectors = map(int, lines[0].split())
    require(channels == 1 and nmax == case["model"]["bath_sites"] - 1, "wrong fixture chain length")
    sectors, position = {}, 1
    for _ in range(nsectors):
        sector = tuple(map(int, lines[position].split()))
        count = int(lines[position + 1])
        energies = finite_values(lines[position + 2].split())
        require(len(sector) == 2 and sector not in sectors and count == len(energies) and count > 0,
                "invalid seed sector")
        sectors[sector] = energies
        position += 3
    require(lines[position].startswith("f "), "missing seed boundary operators")
    require(lines.count("e") == 1 and lines.count("z") == 1, "missing or repeated e/z fixture markers")
    ground = finite_values([lines[lines.index("e") + 1]])[0]
    position = lines.index("z") + 1
    arrays = []
    for _ in range(2):
        require(int(lines[position]) == nmax, "wrong coefficient maximum index")
        arrays.append(finite_values(lines[position + 1:position + nmax + 2]))
        position += nmax + 2
    require(position == len(lines), "unexpected payload after legacy chain trailer")
    p = case["model"]
    compare("chain/hoppings", arrays[0], [wilson_hopping(n, p["Lambda"], p["D"]) for n in range(nmax + 1)],
            checks, atol=5e-15, rtol=5e-15)
    compare("chain/onsites", arrays[1], np.zeros(nmax + 1), checks, atol=5e-15, rtol=5e-15)
    model = flat_band_model(p)
    seed = solve(replace(model, zeta=model.zeta[:1], t=()))
    compare_spectrum({key: values + ground for key, values in sectors.items()}, seed,
                     case["symmetry"], "seed_input", checks)


def parse_report(text):
    blocks, current, sector = [], None, None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if match := re.fullmatch(r"=== Report N=(\d+)", line):
            current, sector = {}, None
            blocks.append((int(match[1]), current))
        elif match := re.fullmatch(r"Sector I=(-?\d+)\s+(-?\d+)", line):
            sector = (int(match[1]), int(match[2]))
            require(current is not None and sector not in current, "duplicate or misplaced report sector")
            current[sector] = []
        elif match := re.fullmatch(r"I=(-?\d+)\s+(-?\d+)\s+n=(\d+)\s+E=(\S+)(?:\s+(.*))?", line):
            require(current is not None and sector == (int(match[1]), int(match[2])), "misplaced report state")
            require(int(match[3]) == len(current[sector]), "noncontiguous report state indices")
            current[sector].append(finite_values([match[4]])[0])
            if match[5]:
                finite_values([field.split("=", 1)[1] for field in match[5].split()])
        else:
            raise ValueError(f"invalid report line: {line}")
    require(blocks and all(block and all(block.values()) for _, block in blocks), "empty spectrum report")
    return blocks


def parse_subspaces(text):
    blocks, current, expected = [], None, None
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if match := re.fullmatch(r"Iteration (\d+)", line):
            require(current is None or (current and expected == len(current)), "empty block or incorrect subspace count")
            current, expected = {}, None
            blocks.append((int(match[1]), current))
        elif match := re.fullmatch(r"len_dm=(\d+)", line):
            require(current is not None and expected is None, "misplaced subspace count")
            expected = int(match[1])
        elif match := re.fullmatch(r"I=(-?\d+)\s+(-?\d+)\s+kept=(\d+)\s+total=(\d+)", line):
            sector = (int(match[1]), int(match[2]))
            require(current is not None and expected is not None and sector not in current, "duplicate/misplaced subspace")
            current[sector] = (int(match[3]), int(match[4]))
        else:
            raise ValueError(f"invalid subspaces line: {line}")
    require(blocks and current and expected == len(current), "empty or truncated subspaces output")
    return blocks


def read_table(path, fields, headers):
    lines = [line.strip() for line in Path(path).read_text().splitlines() if line.strip()]
    require(len(lines) == headers + 1 and all(line.startswith("#") for line in lines[:headers]),
            f"{path}: expected headers and exactly one data row")
    names = lines[headers - 1].lstrip("#").split()
    require(names == fields, f"{path}: unexpected columns {names}")
    values = finite_values(lines[-1].split())
    require(len(values) == len(names), f"{path}: incorrect column count")
    return dict(zip(names, values))


def run_nrg(executable, directory, timeout=120):
    scratch = directory / "scratch"
    scratch.mkdir()
    environment = os.environ.copy()
    environment.update({name: "1" for name in (
        "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "BLIS_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"
    )})
    environment.update(OMP_DYNAMIC="FALSE", MKL_DYNAMIC="FALSE", LC_ALL="C",
                       TMPDIR=str(scratch), TMP=str(scratch), TEMP=str(scratch))
    with (directory / "log").open("w") as log:
        subprocess.run([str(executable), "-w", str(scratch)], cwd=directory, env=environment,
                       stdout=log, stderr=subprocess.STDOUT, check=True, timeout=timeout)
    require((directory / "DONE").is_file() and not (directory / "DONE").is_symlink(), "NRG exited without a fresh DONE")


def run_case(nrg, fixture, temperature, mode, work_root):
    fixture, nrg = Path(fixture).resolve(), Path(nrg).resolve()
    require(nrg.is_file() and os.access(nrg, os.X_OK), f"NRG executable is not runnable: {nrg}")
    case = load_case(fixture)
    provenance = json.loads((fixture / "provenance.json").read_text())
    require(provenance["case_sha256"] == case_digest(case), "model specification changed: regenerate the fixture")
    require(provenance["data_sha256"] == file_digest(fixture / "data"), "fixture data checksum mismatch")
    parameters = parameter_text(case, temperature, mode)
    work_root = Path(work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    directory = Path(tempfile.mkdtemp(prefix="run-", dir=work_root))
    checks = []
    report = {"fixture": fixture.name, "temperature": temperature, "mode": mode,
              "nrg": str(nrg), "numpy": np.__version__, "numpy_configuration": np.show_config(mode="dicts"),
              "checks": checks}
    print(f"Scientific validation: {fixture.name}, T={temperature}, {mode}\nRun directory: {directory}", flush=True)
    try:
        inspect_data(fixture / "data", case, checks)
        shutil.copyfile(fixture / "data", directory / "data")
        (directory / "param").write_text(parameters)
        run_nrg(nrg, directory)
        spectra = parse_report((directory / "report.nrg").read_text())
        subspaces = parse_subspaces((directory / "subspaces.dat").read_text())
        model = flat_band_model(case["model"])
        nbath = len(model.zeta)
        require([index for index, _ in spectra] == [0, *range(nbath - 1)], "wrong prefix report sequence")
        require([index for index, _ in subspaces] == list(range(nbath - 1)), "wrong subspace iteration sequence")
        for prefix, (_, sectors) in enumerate(spectra, start=1):
            solution = solve(replace(model, zeta=model.zeta[:prefix], t=model.t[:prefix - 1]))
            compare_spectrum(sectors, solution, case["symmetry"], f"prefix{prefix}", checks)
            if prefix > 1:
                kept = subspaces[prefix - 2][1]
                expected = dimensions(prefix + 1, case["symmetry"])
                require(set(kept) == set(expected), f"prefix{prefix}: missing/unexpected retained sectors")
                for sector, count in expected.items():
                    require(kept[sector] == (count, count) and len(sectors[sector]) == count,
                            f"prefix{prefix} sector{sector}: truncation or missing eigenpairs: {kept[sector]}, expected {count}")
        custom = read_table(directory / "customfdm", ["T", "I", "n_d", "n_d_ud"], 2)
        compare("customfdm/T", custom.pop("T"), temperature, checks)
        for name, expected in solution.expectations(temperature).items():
            compare(f"customfdm/{name}", custom[name], expected, checks)
        thermo = read_table(directory / "tdfdm", ["T", "E_fdm", "C_fdm", "F_fdm", "S_fdm"], 1)
        compare("tdfdm/T", thermo.pop("T"), temperature, checks)
        for name, expected in solution.thermodynamics(temperature).items():
            compare(f"tdfdm/{name}", thermo[f"{name}_fdm"], expected, checks)
        omega = (2 * np.arange(MATS) + 1) * np.pi * temperature
        greens = solution.greens(temperature, 1j * omega)
        suffixes = [""] if case["symmetry"] == "QS" else ["-u", "-d"]
        for spin, suffix in enumerate(suffixes):
            path = directory / f"spec_FDMmats_dens_A_d-A_d{suffix}.dat"
            rows = np.loadtxt(path, ndmin=2)
            require(rows.shape == (MATS, 3), f"{path.name}: wrong Matsubara table dimensions")
            compare(f"Matsubara{suffix}/omega", rows[:, 0], omega, checks, atol=1e-14, rtol=1e-14)
            compare(f"Matsubara{suffix}/real", rows[:, 1], greens[spin].real, checks)
            compare(f"Matsubara{suffix}/imag", rows[:, 2], greens[spin].imag, checks)
        require(checks and sum(check["count"] for check in checks) > 0, "no numerical comparisons performed")
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (directory / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS: {sum(check['count'] for check in checks)} values; "
          f"largest absolute error {max(check['max_absolute_error'] for check in checks):.3g}", flush=True)
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nrg", type=Path, required=True)
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--temperature", type=float, required=True)
    parser.add_argument("--mode", choices=("rescaled", "absolute"), required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        run_case(args.nrg, args.fixture, args.temperature, args.mode, args.work_root)
    except (ValueError, OSError, KeyError, IndexError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Scientific validation failed: {error}\n")


if __name__ == "__main__":
    main()
