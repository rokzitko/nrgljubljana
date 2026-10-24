"""CLI/native comparisons run only on disposable copies of templates and inputs."""

import argparse
import difflib
import math
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile


def run(command, directory):
    result = subprocess.run(command, cwd=directory, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{command} failed ({result.returncode}):\n{result.stdout}")
    return result.stdout


def compare_files(native, reference, temporary_parent):
    native_bytes, reference_bytes = native.read_bytes(), reference.read_bytes()
    if native_bytes == reference_bytes:
        return
    # Keep evidence outside the disposable fixture tree, even after unwinding.
    artifacts = Path(tempfile.mkdtemp(prefix="workflow-mismatch-", dir=temporary_parent))
    (artifacts / "native").write_bytes(native_bytes)
    (artifacts / "reference").write_bytes(reference_bytes)
    diff = "".join(difflib.unified_diff(
        reference_bytes.decode(errors="replace").splitlines(keepends=True),
        native_bytes.decode(errors="replace").splitlines(keepends=True),
        fromfile=str(reference), tofile=str(native)))
    (artifacts / "diff.log").write_text(diff)
    preview = "".join(diff.splitlines(keepends=True)[:80])
    raise AssertionError(f"Native/reference mismatch: {native}\nSaved files and full diff: {artifacts}\n{preview}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--driver", required=True)
    parser.add_argument("--tools", required=True)
    parser.add_argument("--template", required=True)
    parser.add_argument("--temporary-parent", required=True)
    parser.add_argument("--consumer", help="Optional independently built installed-package consumer")
    parser.add_argument("--current-input", help="Also copy param/Delta.dat from this calculation; requires --consumer and --nmax")
    parser.add_argument("--nmax", type=int, help="Explicit already-resolved chain length for --current-input")
    args = parser.parse_args()
    if args.current_input and (not args.consumer or not args.nmax):
        parser.error("--current-input requires --consumer and explicit --nmax")
    tools = Path(args.tools)
    with tempfile.TemporaryDirectory(prefix="workflow-", dir=args.temporary_parent) as temporary:
        root = Path(temporary)
        host = root / "host"
        host.mkdir()
        (root / "scratch").mkdir()
        (root / "empty").mkdir()
        parameter = """[extra]
eps1=-0.13
U1=0.4
B1=0
[param]
symtype=QS
Lambda=2
band=asymode
adapt=false
polarized=false
density_interpolation=steffen
xmax=4
outputstep=0.125
f_method=ode
Nmax=0
preccpp=384
tridiag_method=lanczos
z=0.7
bandrescale=1
dos=Delta.dat
keep=16
keepmin=1
T=0.01
done=true
data_has_rescaled_energies=false
"""
        for name in ("A", "B"):
            directory = root / name
            (directory / "inputs").mkdir(parents=True)
            shutil.copytree(args.template, directory / "template")
            text = parameter if name == "A" else parameter.replace("Lambda=2", "Lambda=3").replace(
                "z=0.7", "z=0.35").replace("tridiag_method=lanczos", "tridiag_method=rkpw").replace("eps1=-0.13", "eps1=-0.25")
            (directory / "inputs/solver.param").write_text(text)
            (directory / "param").write_text(text)
            (directory / "Delta.dat").write_text("-1 0.15\n-0.5 0.18\n-0.01 0.2\n0.01 0.21\n0.5 0.23\n1 0.25\n")
        a = root / "A"
        shutil.copytree(a / "template", a / "broken-template")
        (a / "broken-template/ham_0.1").write_text("{{unknown_workflow_variable}}\n")
        shutil.copytree(a / "template", a / "complex-template")
        data = a / "complex-template/data.in"
        data.write_text("# COMPLEX\n" + data.read_text())
        (a / "bad-adapt").write_text(parameter.replace("adapt=false", "adapt=false\nhardgap=true\nboundary=0.1").replace("xmax=4", "xmax=80"))
        (a / "unsupported").write_text(parameter.replace("symtype=QS", "symtype=U1"))
        (a / "integral-param").write_text(parameter.replace("f_method=ode", "f_method=integral"))
        (a / "ode-param").write_text(parameter.replace("density_interpolation=steffen", "density_interpolation=linear"))
        sentinel = root / "outside-sentinel"
        sentinel.write_text("Must not be overwritten through an output symlink.\n")
        for filename in ("data", "FSOL.dat", "FSOLNEG.dat"):
            (a / filename).symlink_to(sentinel)

        for name, nmax in (("A", 6), ("B", 5)):
            directory = root / ("cli-" + name)
            shutil.copytree(root / name, directory)
            param = directory / "inputs/solver.param"
            param.write_text(param.read_text().replace("Nmax=0", f"Nmax={nmax}"))
            for branch in ("P", "N"):
                run([str(tools / "adapt"), "--integral", branch, "inputs/solver.param"], directory)
            command = [str(tools / "instantiate"), "--param", "inputs/solver.param", "--template-dir", "template"]
            run(command + ["--generate-temporaries"], directory)
            vectors = {path.name: path.read_bytes() for path in directory.glob("vec.*")}
            assert vectors, "Missing seed eigenvectors"
            for data in vectors.values():
                rows, cols = struct.unpack_from("=II", data)
                values = struct.unpack_from(f"={rows * cols}d", data, 8)
                for row in range(rows):
                    vector = values[row * cols:(row + 1) * cols]
                    assert max(vector, key=abs) > 0, "Noncanonical seed eigenvector sign"
                    assert all(math.copysign(1.0, value) > 0 for value in vector if value == 0)
            for filename in vectors:
                (directory / filename).unlink()
            run(command + ["--diag-seed-only"], directory)
            assert vectors == {path.name: path.read_bytes() for path in directory.glob("vec.*")}, "Seed modes disagree"

        before = {path.relative_to(root) for path in root.rglob("*")}
        run([args.driver, str(root)], host)
        for name in ("A", "B"):
            for filename in ("FSOL.dat", "FSOLNEG.dat", "data"):
                native, cli = root / name / filename, root / ("cli-" + name) / filename
                compare_files(native, cli, args.temporary_parent)
            # Preparation must not leak intermediate chains or staging directories.
            assert not list((root / name).glob("xi*.dat"))
            assert not list((root / name).glob("theta*.dat"))
            for path in (root / name).iterdir():
                if path.is_dir():
                    assert path.relative_to(root) in before, f"Staging directory leaked: {path}"
        run([str(tools / "adapt"), "P", "ode-param"], root / "cli-A")
        compare_files(a / "FSOL-ode.dat", root / "cli-A/FSOL.dat", args.temporary_parent)
        lines = (a / "data").read_text().splitlines()
        energy = float(lines[lines.index("e") + 1])
        assert energy == float((a / "result").read_text()), "Result does not match generated ground energy"
        assert list(host.iterdir()) == [], "Artifacts outside calculation directories"
        assert list((root / "empty").iterdir()) == [], "Failure leaked solver artifacts"
        assert list((root / "scratch").iterdir()) == [], "Scratch not cleaned"
        assert sentinel.read_text() == "Must not be overwritten through an output symlink.\n"
        after = {path.relative_to(root) for path in root.rglob("*")}
        assert all(path.parts[0] in ("A", "B") for path in after - before), "Unexpected output location"
        if args.consumer:
            installed = root / "installed-consumer"
            installed.mkdir()
            shutil.copytree(args.template, installed / "template")
            (installed / "param").write_text(parameter)
            shutil.copyfile(a / "Delta.dat", installed / "Delta.dat")
            run([args.consumer, str(installed), str(root / "scratch")], host)
            assert (installed / "DONE").is_file()
            compare_files(installed / "data", a / "data", args.temporary_parent)
            assert list((root / "scratch").iterdir()) == []
            assert list(host.iterdir()) == []
        if args.current_input:
            current = root / "current-input"
            current.mkdir()
            shutil.copytree(args.template, current / "template")
            for filename in ("param", "Delta.dat"):
                shutil.copyfile(Path(args.current_input) / filename, current / filename)
            original_param = (current / "param").read_bytes()
            reference = root / "current-reference"
            shutil.copytree(current, reference)
            # The native caller supplies nmax; only the disposable CLI reference is edited.
            lines = (reference / "param").read_text().splitlines()
            lines = [line for line in lines if not line.strip().startswith("Nmax=")]
            index = lines.index("[param]") + 1
            lines.insert(index, f"Nmax={args.nmax}")
            (reference / "param").write_text("\n".join(lines) + "\n")
            for branch in ("P", "N"):
                run([str(tools / "adapt"), "--integral", branch], reference)
            run([str(tools / "instantiate")], reference)
            run([args.consumer, str(current), str(root / "scratch"), "--prepare-only", str(args.nmax)], host)
            for filename in ("FSOL.dat", "FSOLNEG.dat", "data"):
                compare_files(current / filename, reference / filename, args.temporary_parent)
            assert (current / "param").read_bytes() == original_param
            assert sorted(path.name for path in current.iterdir()) == ["Delta.dat", "FSOL.dat", "FSOLNEG.dat", "data", "param", "template"]
            assert list(host.iterdir()) == []
            print(f"Current input preparation matches CLI with explicit nmax={args.nmax}")
        print("Workflow contract passed: CLI equivalence, reader validation, A/B/A, failures, state, solve and cleanup")


if __name__ == "__main__":
    main()
