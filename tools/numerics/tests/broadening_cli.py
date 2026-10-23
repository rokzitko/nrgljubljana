"""Optional independent CLI parity oracle. Subprocesses belong only in tests."""

import math
from pathlib import Path
import subprocess
import sys
import tempfile


def table(path):
    return [tuple(map(float, line.split())) for line in path.read_text().splitlines()]


def main():
    native, cli = map(lambda value: str(Path(value).resolve()), sys.argv[1:])
    with tempfile.TemporaryDirectory(prefix="nrg-broadening-parity-") as temporary:
        root = Path(temporary)
        subprocess.run([native, "cli_fixtures", str(root)], check=True)
        for case in sorted(root.iterdir()):
            arguments = (case / "arguments.txt").read_text().splitlines()
            expected = table(case / "native.dat")
            for jobs in (1, 3, 20):
                run = subprocess.run([cli, "--jobs", str(jobs), *arguments], cwd=case, text=True, capture_output=True)
                if run.returncode:
                    raise AssertionError(f"CLI failed in {case.name}: {run.stdout}\n{run.stderr}")
                actual = table(case / "spec.dat")
                if len(actual) != len(expected):
                    raise AssertionError(f"Output length mismatch in {case.name}")
                for (x, y), (ex, ey) in zip(actual, expected):
                    if x != ex or not math.isclose(y, ey, rel_tol=2e-13, abs_tol=2e-14):
                        raise AssertionError(f"Case {case.name}, jobs={jobs}, x={x}: CLI={y}, native={ey}")
        print("10 signed broadening cases matched the CLI at 1, 3 and 20 workers.")


if __name__ == "__main__":
    main()
