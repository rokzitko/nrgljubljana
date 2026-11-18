"""Prepare the small QSTZ spin input; qualification is a separate solver run."""

import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import tempfile


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--basis", type=Path, required=True, help="external QSTZbasis.m (regeneration only)")
    parser.add_argument("--kernel", type=Path, required=True)
    parser.add_argument("--work-root", type=Path, required=True)
    parser.add_argument("--write", action="store_true", help="publish generated data/provenance, not solver outputs")
    args = parser.parse_args()
    source = args.source.resolve()
    fixture = Path(__file__).resolve().parent / "fixtures/qstz_spin"
    args.work_root.mkdir(parents=True, exist_ok=True)
    # The initializer searches '.' and '..'; both must be fresh and owned.
    root = Path(tempfile.mkdtemp(prefix="qstz-prepare-", dir=args.work_root.resolve()))
    candidate = root / "fixture"
    candidate.mkdir()
    for name in ("model.m", "param"):
        shutil.copyfile(fixture / name, candidate / name)
    shutil.copyfile(args.basis, candidate / "QSTZbasis.m")
    script = f'''
NRGDIR={json.dumps(str(source / "nrginit"))};
PACKAGEPATH={{NRGDIR}};
SYMTYPE="runtime";
WriteString[$Output[[1]], "QSTZ_KERNEL_VERSION=", $Version, "\\n"];
If[Get[FileNameJoin[{{NRGDIR, "sneg.m"}}]] === $Failed, Exit[1]];
PARSED=False;
If[Get[FileNameJoin[{{NRGDIR, "initial.m"}}]] === $Failed, Exit[1]];
makedata["data"];
WriteString[$Output[[1]], "QSTZ_INITIALIZATION_SUCCESS\\n"];
Exit[0];
'''
    (candidate / "initialize.m").write_text(script)
    print(f"Preparing: {candidate}", flush=True)
    with (candidate / "initialization.log").open("w") as log:
        subprocess.run([str(args.kernel.resolve()), "-noinit", "-noprompt", "-batchinput", "-batchoutput"],
                       input=script, text=True, cwd=candidate, stdout=log, stderr=subprocess.STDOUT,
                       check=True, timeout=600)
    log = (candidate / "initialization.log").read_text()
    if "QSTZ_INITIALIZATION_SUCCESS" not in log or not (candidate / "data").is_file():
        raise RuntimeError(f"initializer failed; see {candidate / 'initialization.log'}")
    data = candidate / "data"
    data.write_text("\n".join(line.rstrip() for line in data.read_text().splitlines()) + "\n")
    sources = hashlib.sha256()
    for path in sorted((source / "nrginit").rglob("*.m")):
        sources.update(str(path.relative_to(source)).encode() + b"\0" + path.read_bytes())
    provenance = {
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip(),
        "initializer_sources_sha256": sources.hexdigest(),
        "kernel_version": next(line.split("=", 1)[1] for line in log.splitlines()
                               if line.startswith("QSTZ_KERNEL_VERSION=")),
        "basis_sha256": digest(candidate / "QSTZbasis.m"),
        "preparation_script_sha256": digest(Path(__file__)),
        "files": {name: digest(candidate / name) for name in ("data", "param", "model.m")},
    }
    (candidate / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if args.write:
        for name in ("data", "provenance.json"):
            shutil.copyfile(candidate / name, fixture / name)
    print("Qualification required: run validate_qstz_spin.py at T=0.2 and T=0.7.", flush=True)


if __name__ == "__main__":
    main()
