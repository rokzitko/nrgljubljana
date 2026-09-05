"""Regenerate NRG inputs with the repository initializer, then qualify them with ED."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

from validate_siam import (
    case_digest, file_digest, inspect_data, load_case, parameter_text, require, run_case,
)


def prepare(source, fixture, kernel, nrg, work_root, mmax=80, precision=1000, write=False):
    source, fixture = source.resolve(), fixture.resolve()
    case = load_case(fixture)
    require(mmax >= 20 and precision >= 50, "insufficient discretization cutoff/precision")
    work_root.mkdir(parents=True, exist_ok=True)
    # Both working directory and parent are owned by this invocation: the
    # initializer searches '.' and '..' for optional model/operator modules.
    root = Path(tempfile.mkdtemp(prefix="prepare-", dir=work_root.resolve()))
    candidate = root / fixture.name
    candidate.mkdir()
    shutil.copyfile(fixture / "model.json", candidate / "model.json")
    parameters = parameter_text(case, 0.05, "rescaled", mmax, precision)
    (candidate / "param").write_text(parameters)
    nrgdir = source / "nrginit"
    # Use the initializer entry points directly to exclude ~/sneg.m and other
    # local overrides without changing HOME or Mathematica's license discovery.
    script = f'''
NRGDIR={json.dumps(str(nrgdir))};
PACKAGEPATH={{NRGDIR}};
SYMTYPE="runtime";
WriteString[$Output[[1]], "SCIENTIFIC_KERNEL_VERSION=", $Version, "\\n"];
If[Get[FileNameJoin[{{NRGDIR, "sneg.m"}}]] === $Failed, Exit[1]];
PARSED=False;
If[Get[FileNameJoin[{{NRGDIR, "initial.m"}}]] === $Failed, Exit[1]];
makedata["data"];
WriteString[$Output[[1]], "SCIENTIFIC_INITIALIZATION_SUCCESS\\n"];
Exit[0];
'''
    print(f"Preparing {fixture.name}: {candidate}", flush=True)
    with (candidate / "initialization.log").open("w") as log:
        subprocess.run([str(kernel.resolve()), "-noinit", "-noprompt", "-batchinput", "-batchoutput"],
                       input=script, text=True, cwd=candidate, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=600)
    log = (candidate / "initialization.log").read_text()
    require("SCIENTIFIC_INITIALIZATION_SUCCESS" in log and (candidate / "data").is_file(),
            f"initialization did not complete; see {candidate / 'initialization.log'}")
    data = candidate / "data"
    data.write_text("\n".join(line.rstrip() for line in data.read_text().splitlines()) + "\n")
    checks = []
    inspect_data(candidate / "data", case, checks)
    versions = [line.split("=", 1)[1] for line in log.splitlines() if line.startswith("SCIENTIFIC_KERNEL_VERSION=")]
    require(len(versions) == 1, "missing kernel version")
    sources = hashlib.sha256()
    for path in sorted(nrgdir.rglob("*.m")):
        sources.update(str(path.relative_to(source)).encode() + b"\0" + path.read_bytes())
    revision = subprocess.run(["git", "rev-parse", "HEAD"], cwd=source, text=True,
                              stdout=subprocess.PIPE, check=True).stdout.strip()
    provenance = {
        "case_sha256": case_digest(case), "data_sha256": file_digest(candidate / "data"),
        "source_revision": revision, "initializer_sources_sha256": sources.hexdigest(),
        "preparation_script_sha256": file_digest(Path(__file__)),
        "kernel_version": versions[0], "mMAX": mmax, "prec": precision,
        "parameters": parameters,
    }
    (candidate / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    # Regeneration is judged by physical results, not eigenvector phases or bytes.
    for temperature in (0.05, 0.2):
        for mode in ("rescaled", "absolute"):
            run_case(nrg, candidate, temperature, mode, root / "validation")
    if write:
        for name in ("data", "provenance.json"):
            shutil.copyfile(candidate / name, fixture / name)
        print(f"Updated qualified fixture: {fixture}", flush=True)
    else:
        require((fixture / "data").is_file(), "no committed fixture to check")
        inspect_data(fixture / "data", case, [])
        provenance = json.loads((fixture / "provenance.json").read_text())
        require(provenance["case_sha256"] == case_digest(case), "committed case provenance mismatch")
        require(provenance["data_sha256"] == file_digest(fixture / "data"), "committed fixture checksum mismatch")
        print(f"Regeneration qualified without modifying {fixture}", flush=True)
    return candidate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--nrg", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--mmax", type=int, default=80)
    parser.add_argument("--precision", type=int, default=1000)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true", help="replace fixtures only after validation passes")
    action.add_argument("--check", action="store_true", help="validate regeneration without changing source files")
    args = parser.parse_args()
    try:
        prepare(args.source, args.fixture, args.kernel, args.nrg, args.work_root,
                args.mmax, args.precision, args.write)
    except (ValueError, OSError, KeyError, IndexError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Fixture preparation failed: {error}\n")


if __name__ == "__main__":
    main()
