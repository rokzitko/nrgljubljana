"""Generate the same three-channel bath in two spin symmetries and extend it once."""

import math
import os
from pathlib import Path
import re
import subprocess
import sys
import tempfile


def require(ok, message):
    if not ok:
        raise ValueError(message)


def run_case(case, symmetry, package, kernel, nrg):
    case.mkdir(parents=True)
    operators = ""
    if symmetry == "QSZTZ":
        # One central (Tz=0) spin component and a single spin-up fermion.
        operators = "ops=SZf0 A0\nspecchit=SZf0-SZf0\nspecgt=A0-A0\ngtp=1\nchitp=1\nprec_xy=16\n"
        (case / "modeloperators.m").write_text('mtDoubletOp["A0", If[#2 == 0 && #3 == UP, f[#1, 0, UP], 0]&]\n')
    (case / "param").write_text(f"""[extra]
Gamma=0.1
[param]
symtype={symmetry}
model=model.m
Ninit=0
Nmax=1
Lambda=2
band=flat
discretization=Z
tri=old
mMAX=80
prec=100
diag=dsyevd
diag_mode=serial
mult=blas
keep=4096
keepenergy=-1
calc0=true
silent=false
done=true
prec_td=16
width_td=24
{operators}""")
    (case / "model.m").write_text(f"""def3ch[0];
{symmetry.lower()}ImpurityCharges = {{}};
H = H0 /. f[ca_, ch_, sp_] :> f[ca, 1 - ch, sp];
""")
    # Run the installed initializer directly, pinning PACKAGEPATH instead of
    # the launcher's ~/ and ~/nrg overrides. basis.m still does real module
    # dispatch; neither the fresh case nor its parent contains a basis file.
    script = '''NRGDIR = Environment["NRG_INIT_DIR"];
PACKAGEPATH = {NRGDIR};
SYMTYPE = "runtime";
Get["sneg.m", Path -> PACKAGEPATH];
Get["initial.m", Path -> PACKAGEPATH];
makedata["data"];
Print["INITIALIZER COMPLETE"];
Exit[0];
'''

    def run(command, **kwargs):
        result = subprocess.run(command, cwd=case, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, **kwargs)
        if result.returncode:
            raise RuntimeError(f"{command[0]} exited {result.returncode}:\n{result.stdout}")
        return result.stdout

    log = run([kernel, "-noinit", "-noprompt", "-batchinput", "-batchoutput"], input=script,
              env=dict(os.environ, NRG_INIT_DIR=str(package)))
    require("INITIALIZER COMPLETE" in log and (case / "data").is_file() and
            f"Loading module {symmetry}basis.m" in log, f"incomplete initializer: {log}")
    log = run([nrg])
    require((case / "DONE").is_file(), log)
    require(re.search(r"Number of states \(multiplicity taken into account\): 64\b", log), log)
    # Six spinful orbitals have 4096 states; the number of spin multiplets
    # equals the number of states with Sz=0 or Sz=1/2, binomial(13,6).
    kept = 1716 if symmetry == "QSTZ" else 4096
    counts = re.findall(r"nrkept=(\d+) nrkeptmult=(\d+) nrall=(\d+) nrallmult=(\d+)", log)
    require(counts == [(str(kept), "4096", str(kept), "4096")], f"{symmetry} retention: {counts}")
    print(f"{symmetry}: 64 -> 4096 physical states ({kept} retained eigenvectors)")
    header, seed, step = (case / "td").read_text().splitlines()
    fields = header.lstrip("#").split()
    rows = []
    for line in (seed, step):
        values = list(map(float, line.split()))
        require(len(fields) == len(values) and all(map(math.isfinite, values)), f"invalid td row: {line}")
        rows.append(dict(zip(fields, values)))
    if symmetry == "QSZTZ":
        # The central orbital has E=0 in the seed, then equal weights at +/-xi0
        # in the dimer. GT is -f'(E) for one spin; T*chi_z also has interlevel
        # matrix elements of the local Sz. Output order is increasing T (step, seed).
        xi0 = 1 / (math.sqrt(7) * math.log(2))  # flat Z discretization, Lambda=2
        for filename, spin in (("chit_CHIT_dens_SZf0-SZf0.dat", True), ("gt_GT_dens_A0-A0.dat", False)):
            curve = (case / filename).read_text().splitlines()
            require(len(curve) == 2, f"{filename}: expected seed and one extension")
            for line, td, energy in zip(curve, reversed(rows), (xi0, 0)):
                temp, value = map(float, line.split())
                require(abs(temp - td["T"]) < 1e-10, f"{filename}: wrong temperature {temp}")
                f = 1 / (1 + math.exp(energy / temp))
                if spin:
                    expected = f * (1 - f) / 4 + temp * (1 - 2 * f) / (8 * energy) if energy else 1 / 8
                else:
                    expected = f * (1 - f) / temp
                require(abs(value - expected) < 1e-10, f"{filename} at T={temp}: {value} != {expected}")
                print(f"{filename} at T={temp:.6g}: {value:.12g} (absolute error {abs(value - expected):.3g})")
    return rows


def main():
    cmake, install_tree, destination, kernel, nrg = sys.argv[1:]
    with tempfile.TemporaryDirectory(prefix="qsztz-", dir=Path.cwd()) as tmp:
        root = Path(tmp)
        # DESTDIR also confines builds configured with an absolute destination.
        subprocess.run([cmake, "--install", install_tree, "--prefix", "/stage"], check=True,
                       env=dict(os.environ, DESTDIR=tmp))
        package = root / (Path("/stage") / destination).relative_to("/")
        require((package / "QSZTZbasis.m").is_file(), "QSZTZbasis.m missing from CMake installation")
        reference, additive = [run_case(root / "runs" / sym, sym, package, kernel, nrg)
                               for sym in ("QSTZ", "QSZTZ")]
        for i, (ref, row) in enumerate(zip(reference, additive)):
            errors = {k: abs(row[k] - ref[k]) for k in ("T", "<Sz^2>", "<Tz^2>", "<Q>", "<Q^2>",
                                                         "<E>", "<E^2>", "C", "F", "S")}
            errors.update({k: abs(row[k]) for k in ("<Sz>", "<Tz>")})
            require(max(errors.values()) < 1e-10, f"TD row {i} absolute differences: {errors}")
            print(f"TD row {i}: max absolute difference {max(errors.values()):.3g}")
        require(abs(additive[0]["S"] - math.log(64)) < 1e-10, "degenerate seed entropy != ln(64)")


if __name__ == "__main__":
    main()
