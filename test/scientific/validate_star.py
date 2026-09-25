"""Qualify scalar frontends against independent asymmetric finite-star physics."""

from __future__ import annotations

import argparse
from configparser import ConfigParser, Error as ConfigError
import hashlib
from io import StringIO
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

import numpy as np

from ed_siam import solve
from finite_star import FIXTURE, jacobi, load_star, reference_greens, reference_model, star_tables
from validate_siam import (
    MATS, case_digest, compare, compare_spectrum, dimensions, file_digest, finite_values,
    parameter_text, parse_report, parse_subspaces, read_table, require, run_nrg,
)


SOURCE_HASH_SCOPE = "nrginit/**/*.m excluding paths with a mix* component; SHA256 of canonical source-hash map"


def bath_count(backend):
    require(backend in ("legacy", "rkpw"), "unknown star backend")
    return 3 if backend == "legacy" else 4


def parameters(star, backend, D, temperature, mode, tri=None, template=False):
    # Reuse the solver controls, not the flat-band reference or fixture schema.
    case = {"symmetry": "QS", "model": {
        "epsilon_d": star["epsilon_d"], "U": 0, "B": 0, "Gamma": 1 / D,
        "D": D, "Lambda": 2, "z": 1, "bath_sites": bath_count(backend),
    }}
    config = ConfigParser()
    config.optionxform = str
    config.read_string(parameter_text(case, temperature, mode, backend=backend))
    config["param"].update({
        "Ninit": "1", "keep": "1024", "mMAX": "1", "prec": "80", "preccpp": "2000",
        "ops": "I A_d n_d", "fdmexpvn": "1", "adapt": "false", "rescalexi": "false",
        "nrgchain_tables_load": "true",
    })
    if tri is not None:
        config["param"]["tri"] = tri
    if template:
        config["param"]["options"] = "GENERATE_TEMPLATE"
        # GENERATE_TEMPLATE disables reconstruction; avoid numeric RKPW export guards.
        config["param"]["tri"] = "old"
    else:
        config["param"]["hook_pre_lanczosinit"] = "import_star.m"
    output = StringIO()
    config.write(output, space_around_delimiters=False)
    return output.getvalue()


def table_texts(star, D):
    return {name: "".join(f"{value:.17g}\n" for value in values) for name, values in star_tables(star, D).items()}


def write_tables(directory, star, D):
    for name, text in table_texts(star, D).items():
        (directory / name).write_text(text)


def table_block(marker, arrays):
    return marker + "\n" + "".join(
        str(len(values) - 1) + "\n" + "".join(f"{value:.17g}\n" for value in values)
        for values in arrays)


def t_block(star, D):
    tables = star_tables(star, D)
    return table_block("T", [tables[name] for name in ("de_pos.dat", "de_neg.dat", "du_pos.dat", "du_neg.dat")])


def check_coefficients(xi, zeta, star, count, checks, label="chain"):
    diagonal, hopping = jacobi(star)
    compare(label + "/xi", xi, np.r_[hopping, 0][:count], checks, atol=5e-15, rtol=5e-15)
    compare(label + "/zeta", zeta, diagonal[:count], checks, atol=5e-15, rtol=5e-15)
    if count == 4:
        require(xi[-1] == 0, label + ": finite support must have exact terminal xi[3]=0")


def inspect_data(path, star, backend, D, checks, runtime=False):
    text = Path(path).read_text()
    require(re.search(r"^#!9\s*$", text, re.M), "expected real #!9 data")
    require("COMPLEX" not in text and re.search(r"^# symtype\s+QS\s*$", text, re.M), "expected real QS seed")
    require(re.search(r"^# SCALE\s+1(?:\.0*)?\s*$", text, re.M), "seed must have SCALE 1")
    lines = [line.strip() for line in text.splitlines() if line.strip() and not line.startswith("#")]
    channels, nmax, nsectors = map(int, lines[0].split())
    require(channels == 1 and nmax + 1 == bath_count(backend), "wrong data chain length")
    sectors, position = {}, 1
    for _ in range(nsectors):
        sector = tuple(map(int, lines[position].split()))
        count = int(lines[position + 1])
        values = finite_values(lines[position + 2].split())
        require(len(sector) == 2 and sector not in sectors and count == len(values), "invalid seed sector")
        sectors[sector] = values
        position += 3
    require(lines[position] == "f 0 0", "missing boundary operator")
    require(all(lines.count(marker) == 1 for marker in ("e", "z", "dA_d", "sI", "sn_d")),
            "missing or repeated data blocks")
    ground = finite_values([lines[lines.index("e") + 1]])[0]
    compare_spectrum({sector: values + ground for sector, values in sectors.items()},
                     solve(reference_model(star, 2)), "QS", "seed64", checks)
    position = lines.index("z") + 1
    arrays = []
    for _ in range(2):
        count = int(lines[position]) + 1
        require(count == (2 if runtime else nmax + 1), "wrong coefficient count")
        arrays.append(finite_values(lines[position + 1:position + 1 + count]))
        position += count + 1
    check_coefficients(*arrays, star, len(arrays[0]), checks, "data/z")
    if runtime:
        require(lines[position] == "T" and lines.count("T") == 1, "missing runtime T block")
        expected = t_block(star, D).splitlines()
        actual = lines[position:]
        compare("data/T", finite_values(actual[1:]), finite_values(expected[1:]), checks,
                atol=5e-15, rtol=5e-15)
        position = len(lines)
    require(position == len(lines), "unexpected data trailer")


def template_files(text):
    return {"data.in", *re.findall(r"^DIAG\s+(\S+)\s*$", text, re.M),
            *re.findall(r"^(op\.\S+)\s*$", text, re.M)}


def verify_fixture(fixture, backend, lane):
    star = load_star(fixture)
    asset = fixture / ("template" if lane == "instantiate" else backend)
    provenance = json.loads((asset / "provenance.json").read_text())
    fields = {
        "kernel_version", "source_revision", "initializer_sources_sha256", "initializer_source_files",
        "source_hash_scope", "preparation_script_sha256", "hook_sha256", "entrypoint_sha256",
        "import_check_sha256", "parameters", "precision", "mMAX", "physical_input_files", "kind",
        "backend", "star_sha256", "star_file_sha256", "files",
    }
    require(isinstance(provenance, dict) and fields <= provenance.keys(),
            "missing generation provenance fields")
    require(provenance["star_sha256"] == case_digest(star), "physical star provenance mismatch")
    require(provenance["star_file_sha256"] == file_digest(fixture / "star.json"), "star file checksum mismatch")
    required = {"data"}
    if lane == "instantiate":
        text = (asset / "data.in").read_text()
        required = template_files(text)
        require(len(required) > 1, "missing template Hamiltonians")
        require(provenance["kind"] == "symbolic-template", "incorrect template provenance")
    else:
        require(provenance["backend"] == backend and provenance["kind"] == "cpp-seed",
                "incorrect seed backend provenance")
    require(isinstance(provenance["files"], dict) and set(provenance["files"]) == required, "incomplete fixture manifest")
    for name in required:
        require(Path(name).name == name, "invalid fixture filename")
        require(provenance["files"][name] == file_digest(asset / name), f"fixture checksum mismatch: {name}")
    require({path.name for path in asset.iterdir()} == required | {"provenance.json"},
            "undeclared fixture asset")
    for name in ("kernel_version", "parameters"):
        require(isinstance(provenance[name], str) and bool(provenance[name].strip()),
                "missing generation provenance: " + name)
    require(isinstance(provenance["source_revision"], str)
            and re.fullmatch(r"(?:[0-9a-f]{40}|[0-9a-f]{64})", provenance["source_revision"]),
            "invalid source revision")
    for name in ("initializer_sources_sha256", "preparation_script_sha256", "hook_sha256",
                 "entrypoint_sha256", "import_check_sha256"):
        require(isinstance(provenance[name], str) and re.fullmatch(r"[0-9a-f]{64}", provenance[name]),
                "invalid provenance hash: " + name)
    require(provenance["source_hash_scope"] == SOURCE_HASH_SCOPE, "incorrect initializer source hash scope")
    sources = provenance["initializer_source_files"]
    require(isinstance(sources, dict) and sources, "missing initializer source map")
    for name, digest in sources.items():
        path = Path(name)
        require(not path.is_absolute() and str(path) == name and ".." not in path.parts
                and path.parts[0] == "nrginit" and path.suffix == ".m"
                and not any(part.lower().startswith("mix") for part in path.parts),
                "invalid initializer source path")
        require(isinstance(digest, str) and re.fullmatch(r"[0-9a-f]{64}", digest), "invalid initializer source hash")
    require({"nrginit/initial.m", "nrginit/sneg.m", "nrginit/wilson.m"} <= sources.keys(),
            "missing initializer entry points")
    require(provenance["initializer_sources_sha256"] == case_digest(sources), "initializer source manifest mismatch")
    expected_inputs = {name: hashlib.sha256(text.encode()).hexdigest() for name, text in table_texts(star, 2).items()}
    require(provenance["physical_input_files"] == expected_inputs, "physical input hash/unit contract mismatch")
    require(type(provenance["precision"]) is int and provenance["precision"] == 80
            and type(provenance["mMAX"]) is int and provenance["mMAX"] == 1, "incorrect preparation precision/cutoff")
    generation_backend = provenance["backend"]
    require(generation_backend in ("legacy", "rkpw"), "invalid generation backend")
    try:
        config = ConfigParser(interpolation=None)
        config.optionxform = str
        config.read_string(provenance["parameters"])
        require(set(config.sections()) == {"extra", "param"} and not config.defaults(), "incorrect parameter sections")
        p, extra = config["param"], config["extra"]
        strings = {"model": "SIAM", "variant": "EPS", "symtype": "QS", "band": "flat", "discretization": "Z",
                   "wilsonchain": "legacy", "tri": "old" if lane == "instantiate" else "cpp",
                   "tridiag_method": "lanczos" if generation_backend == "legacy" else "rkpw"}
        for name, expected in strings.items():
            require(p[name] == expected, "incorrect generation parameter: " + name)
        for name, expected in {"Ninit": 1, "Nmax": bath_count(generation_backend) - 1,
                               "mMAX": provenance["mMAX"], "prec": provenance["precision"], "preccpp": 2000}.items():
            require(p.getint(name) == expected, "incorrect generation parameter: " + name)
        for name, expected in {"bandrescale": 2, "Lambda": 2, "z": 1, "T": 0.05}.items():
            require(p.getfloat(name) == expected, "incorrect generation parameter: " + name)
        for name, expected in {"eps": star["epsilon_d"], "U": 0, "B": 0, "Gamma": 0.5}.items():
            require(extra.getfloat(name) == expected, "incorrect generation impurity parameter: " + name)
        for name in ("data_has_rescaled_energies", "polarized", "substeps", "rescalexi", "adapt", "absolute"):
            require(p.getboolean(name) is False, "incorrect generation parameter: " + name)
        require(p.getboolean("nrgchain_tables_load") is True, "generation must load the saved star")
        require(p.get("options", "").split() == (["GENERATE_TEMPLATE"] if lane == "instantiate" else []),
                "incorrect generation options")
        require(p.get("hook_pre_lanczosinit", "") == ("" if lane == "instantiate" else "import_star.m"),
                "incorrect generation hook")
    except (ConfigError, KeyError) as error:
        raise ValueError(f"invalid generation parameters: {error}") from error
    return star, asset


def run_tool(command, directory, name):
    with (directory / (name + ".log")).open("w") as log:
        subprocess.run([str(arg) for arg in command], cwd=directory, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=180)


def qualify_solver(nrg, directory, star, backend, D, temperature, checks, runtime=False):
    inspect_data(directory / "data", star, backend, D, checks, runtime)
    run_nrg(Path(nrg).resolve(), directory)
    nbath = bath_count(backend)
    if runtime:
        log = (directory / "log").read_text()
        arrays = []
        for name in ("xi", "zeta"):
            matches = re.findall(r"^\s*" + name + r"\((\d+)\)=(\S+)\s*$", log, re.M)
            require([int(i) for i, _ in matches] == list(range(nbath)), "missing runtime " + name)
            arrays.append(D * finite_values([value for _, value in matches]))
        check_coefficients(*arrays, star, nbath, checks, "runtime/physical")
    spectra = parse_report((directory / "report.nrg").read_text())
    subspaces = parse_subspaces((directory / "subspaces.dat").read_text())
    require([i for i, _ in spectra] == [0, *range(1, nbath - 1)], "wrong Ninit=1 spectrum sequence")
    require([i for i, _ in subspaces] == list(range(1, nbath - 1)), "wrong subspace sequence")
    for prefix, (_, sectors) in enumerate(spectra, start=2):
        solution = solve(reference_model(star, prefix))
        compare_spectrum(sectors, solution, "QS", f"prefix{prefix}", checks)
        require(sum(spin * len(values) for (_, spin), values in sectors.items()) == 4 ** (prefix + 1),
                "wrong physical state count")
        if prefix > 2:
            kept = subspaces[prefix - 3][1]
            expected = dimensions(prefix + 1, "QS")
            require(set(kept) == set(expected), "missing retained sectors")
            require(all(kept[key] == (count, count) for key, count in expected.items()),
                    "unexpected many-body truncation")
    custom = read_table(directory / "customfdm", ["T", "I", "n_d"], 2)
    compare("customfdm/T", custom.pop("T"), temperature, checks)
    for name in ("I", "n_d"):
        compare("customfdm/" + name, custom[name], solution.expectations(temperature)[name], checks)
    thermo = read_table(directory / "tdfdm", ["T", "E_fdm", "C_fdm", "F_fdm", "S_fdm"], 1)
    compare("tdfdm/T", thermo.pop("T"), temperature, checks)
    for name, value in solution.thermodynamics(temperature).items():
        compare("tdfdm/" + name, thermo[name + "_fdm"], value, checks)
    omega = (2 * np.arange(MATS) + 1) * np.pi * temperature
    expected = reference_greens(star, 1j * omega, nbath)
    rows = np.loadtxt(directory / "spec_FDMmats_dens_A_d-A_d.dat", ndmin=2)
    require(rows.shape == (MATS, 3), "wrong Matsubara shape")
    compare("G/omega", rows[:, 0], omega, checks, atol=1e-14, rtol=1e-14)
    compare("G/real", rows[:, 1], expected.real, checks)
    compare("G/imag", rows[:, 2], expected.imag, checks)


def run_case(nrg, fixture, backend, lane, work_root, D=2, temperature=0.05, mode="rescaled",
             nrgchain=None, instantiate=None, kernel=None, source=None):
    fixture = Path(fixture).resolve()
    require(lane in ("runtime", "nrgchain", "instantiate", "nrginit"), "unknown star frontend")
    bath_count(backend)
    if lane == "nrginit":
        star, asset = load_star(fixture), None
    else:
        star, asset = verify_fixture(fixture, backend, lane)
    work_root = Path(work_root).resolve()
    work_root.mkdir(parents=True, exist_ok=True)
    # Own the parent too, excluding initializer module lookup in caller directories.
    root = Path(tempfile.mkdtemp(prefix=f"{lane}-{backend}-", dir=work_root))
    directory = root / "run"
    directory.mkdir()
    checks = []
    nbath = bath_count(backend)
    report = {"backend": backend, "frontend": lane, "D": D, "temperature": temperature, "mode": mode,
              "star_sha256": case_digest(star),
              "reference": "direct four-pole star" if nbath == 4 else "independent three-site Jacobi prefix",
              "star_poles": 4, "bath_sites": nbath, "omitted_bath_sites": 4 - nbath,
              "seed_states": 64, "final_states": 4 ** (nbath + 1), "discarded_many_body_states": 0,
              "nrg": str(Path(nrg).resolve()), "numpy": np.__version__, "checks": checks}
    print(f"Star qualification: {directory}", flush=True)
    try:
        tri = "cpp" if lane == "runtime" else None
        (directory / "param").write_text(parameters(star, backend, D, temperature, mode, tri))
        write_tables(directory, star, D)
        if lane in ("runtime", "nrgchain"):
            inspect_data(asset / "data", star, backend, 2, checks, runtime=True)
            text = (asset / "data").read_text()
            if lane == "runtime":
                text = text.split("\nT\n")[0] + "\n" + t_block(star, D)
            else:
                run_tool([Path(nrgchain).resolve(), "l"], directory, "nrgchain")
                xi, zeta = [np.loadtxt(directory / name, ndmin=1) for name in ("xi.dat", "zeta.dat")]
                check_coefficients(xi, zeta, star, nbath, checks, "nrgchain")
                text = text.split("\nz\n")[0] + "\n" + table_block("z", [xi, zeta])
            (directory / "data").write_text(text)
        elif lane == "instantiate":
            shutil.copytree(asset, directory / "template")
            run_tool([Path(instantiate).resolve(), "--generate-temporaries"], directory, "instantiate")
            for suffix in (".dat", "1.dat"):
                xi, zeta = [np.loadtxt(directory / (name + suffix), ndmin=1) for name in ("xi", "zeta")]
                check_coefficients(xi, zeta, star, nbath, checks, "instantiate/" + suffix)
                compare("theta/" + suffix, np.loadtxt(directory / ("theta" + suffix), ndmin=1),
                        star_tables(star, D)["theta.dat"], checks, atol=5e-15, rtol=5e-15)
        else:
            from prepare_star import initialize
            initialize(Path(source), Path(kernel), directory)
        compare("theta", np.loadtxt(directory / "theta.dat", ndmin=1),
                star_tables(star, D)["theta.dat"], checks, atol=5e-15, rtol=5e-15)
        qualify_solver(nrg, directory, star, backend, D, temperature, checks, lane == "runtime")
        report["status"] = "passed"
    except Exception as error:
        report.update(status="failed", error=str(error))
        raise
    finally:
        (directory / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(f"PASS {lane}/{backend}: {sum(c['count'] for c in checks)} values; "
          f"{report['final_states']} states, bath prefix {nbath}/4", flush=True)
    return directory


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fixture", type=Path, default=FIXTURE)
    parser.add_argument("--backend", choices=("legacy", "rkpw"), required=True)
    parser.add_argument("--lane", choices=("runtime", "nrgchain", "instantiate", "nrginit"), required=True)
    parser.add_argument("--nrg", type=Path, required=True)
    parser.add_argument("--nrgchain", type=Path)
    parser.add_argument("--instantiate", type=Path)
    parser.add_argument("--kernel", type=Path)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--D", type=float, choices=(1, 2), default=2)
    parser.add_argument("--temperature", type=float, default=0.05)
    parser.add_argument("--mode", choices=("rescaled", "absolute"), default="rescaled")
    args = parser.parse_args()
    for lane, option in (("nrgchain", "nrgchain"), ("instantiate", "instantiate"), ("nrginit", "kernel")):
        if args.lane == lane and getattr(args, option) is None:
            parser.error(f"--{option} is required for {lane}")
    try:
        run_case(**vars(args))
    except (ValueError, OSError, KeyError, IndexError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Star qualification failed: {error}\n")


if __name__ == "__main__":
    main()
