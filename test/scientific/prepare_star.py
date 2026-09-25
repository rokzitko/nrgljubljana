"""Generate new Ninit=1 finite-star assets using source-only Mathematica entry points."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
import shutil
import subprocess
import tempfile

from finite_star import FIXTURE, load_star
from validate_siam import case_digest, file_digest, require
from validate_star import SOURCE_HASH_SCOPE, parameters, run_case, template_files, write_tables


def initialize(source, kernel, directory, template=False):
    source, kernel = source.resolve(), kernel.resolve()
    hook = Path(__file__).with_name("import_star.m")
    import_check = Path(__file__).with_name("check_star_import.m")
    if not template:
        shutil.copyfile(hook, directory / "import_star.m")
    script = f'''
NRGDIR={json.dumps(str(source / "nrginit"))};
PACKAGEPATH={{NRGDIR}};
SYMTYPE="runtime";
WriteString[$Output[[1]], "SCIENTIFIC_KERNEL_VERSION=", $Version, "\\n"];
If[Get[{json.dumps(str(import_check))}] === $Failed, Exit[1]];
If[Get[FileNameJoin[{{NRGDIR, "sneg.m"}}]] === $Failed, Exit[1]];
PARSED=False;
If[Get[FileNameJoin[{{NRGDIR, "initial.m"}}]] === $Failed, Exit[1]];
makedata["data"];
WriteString[$Output[[1]], "SCIENTIFIC_INITIALIZATION_SUCCESS\\n"];
Exit[0];
'''
    (directory / "entrypoint.m").write_text(script)
    with (directory / "initialization.log").open("w") as log:
        subprocess.run([str(kernel), "-noinit", "-noprompt", "-batchinput", "-batchoutput"],
                       input=script, text=True, cwd=directory, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=600)
    log = (directory / "initialization.log").read_text()
    versions = re.findall(r"^SCIENTIFIC_KERNEL_VERSION=(.+)$", log, re.M)
    output = directory / ("data.in" if template else "data")
    require("SCIENTIFIC_INITIALIZATION_SUCCESS" in log and output.is_file() and len(versions) == 1,
            f"incomplete initializer; see {directory / 'initialization.log'}")
    require("SCIENTIFIC_STAR_IMPORT_EXACT" in log, "finite-star import precision check did not run")
    output.write_text("\n".join(line.rstrip() for line in output.read_text().splitlines()) + "\n")
    source_files = {}
    for path in sorted((source / "nrginit").rglob("*.m")):
        if any(part.lower().startswith("mix") for part in path.relative_to(source).parts):
            continue
        name = str(path.relative_to(source))
        source_files[name] = file_digest(path)
    provenance = {
        "kernel_version": versions[0],
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip(),
        "initializer_sources_sha256": case_digest(source_files),
        "initializer_source_files": source_files,
        "source_hash_scope": SOURCE_HASH_SCOPE,
        "preparation_script_sha256": file_digest(Path(__file__)),
        "hook_sha256": file_digest(hook), "entrypoint_sha256": file_digest(directory / "entrypoint.m"),
        "import_check_sha256": file_digest(import_check),
        "parameters": (directory / "param").read_text(), "precision": 80, "mMAX": 1,
        "physical_input_files": {name: file_digest(directory / name) for name in
                                 ("de_pos.dat", "de_neg.dat", "du_pos.dat", "du_neg.dat", "theta.dat")},
    }
    (directory / "generation.json").write_text(json.dumps(provenance, indent=2) + "\n")
    return provenance


def prepare(source, fixture, kernel, nrg, instantiate, backend, work_root, template=False, write=False):
    star = load_star(fixture)
    work_root.mkdir(parents=True, exist_ok=True)
    root = Path(tempfile.mkdtemp(prefix="prepare-star-", dir=work_root.resolve()))
    generation = root / "generation"
    generation.mkdir()
    print(f"Generating finite-star {'template' if template else backend + ' seed'}: {generation}", flush=True)
    (generation / "param").write_text(parameters(star, backend, 2, 0.05, "rescaled",
                                                tri=None if template else "cpp", template=template))
    write_tables(generation, star, 2)
    provenance = initialize(source, kernel, generation, template)
    candidate = root / "fixture"
    candidate.mkdir()
    shutil.copyfile(fixture / "star.json", candidate / "star.json")
    asset_name = "template" if template else backend
    asset = candidate / asset_name
    asset.mkdir()
    names = {"data"}
    if template:
        text = (generation / "data.in").read_text()
        names = template_files(text)
    for name in sorted(names):
        require(Path(name).name == name and not name.lower().startswith("mix"), "unexpected generated asset")
        text = (generation / name).read_text()
        (asset / name).write_text("\n".join(line.rstrip() for line in text.splitlines()) + "\n")
    provenance.update(
        kind="symbolic-template" if template else "cpp-seed", backend=backend,
        star_sha256=case_digest(star), star_file_sha256=file_digest(candidate / "star.json"),
        files={name: file_digest(asset / name) for name in sorted(names)},
    )
    (asset / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    for temperature in (0.05, 0.2):
        for mode in ("rescaled", "absolute"):
            run_case(nrg, candidate, backend, "instantiate" if template else "runtime", root / "qualification",
                     temperature=temperature, mode=mode, instantiate=instantiate)
    if write:
        destination = fixture / asset_name
        if destination.exists():
            require({path.name for path in destination.iterdir()} == names | {"provenance.json"},
                    "refusing to overwrite an unexpected asset set")
        destination.mkdir(exist_ok=True)
        for name in sorted(names | {"provenance.json"}):
            shutil.copyfile(asset / name, destination / name)
        print(f"Published {len(names)} qualified text assets plus provenance in {destination}", flush=True)
    return candidate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--fixture", type=Path, default=FIXTURE)
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--nrg", type=Path, required=True)
    parser.add_argument("--instantiate", type=Path, required=True)
    parser.add_argument("--backend", choices=("legacy", "rkpw"), required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--template", action="store_true")
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--write", action="store_true", help="publish star assets only after qualification")
    action.add_argument("--check", action="store_true", help="generate and qualify without modifying fixtures")
    args = vars(parser.parse_args())
    args.pop("check")
    try:
        prepare(**args)
    except (ValueError, OSError, KeyError, IndexError, subprocess.SubprocessError) as error:
        parser.exit(1, f"Star preparation failed: {error}\n")


if __name__ == "__main__":
    main()
