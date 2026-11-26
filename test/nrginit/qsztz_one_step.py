"""Generate the same three-channel bath in two spin symmetries and extend it once."""

import math
import os
from pathlib import Path
import re
import struct
import subprocess
import sys
import tempfile


XI0 = 1 / (math.sqrt(7) * math.log(2))  # first flat-band Z hopping, Lambda=2
FDM_T = 0.2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def close(label, value, expected):
    require(abs(value - expected) < 1e-10, f"{label}: {value} != {expected}")


def signed_weights(case, prefix, op):
    path = case / f"{prefix}_FDM_dens_{op}-{op}.bin"
    weights = [0., 0.]  # negative and positive frequency, without taking abs(weight)
    for omega, weight in struct.iter_unpack("dd", path.read_bytes()):
        require(math.isfinite(omega) and omega != 0 and math.isfinite(weight), f"invalid bin in {path}")
        require(math.copysign(1, omega) * weight >= -1e-12, f"{op}: wrong bin sign at {omega}: {weight}")
        weights[omega > 0] += weight
    print(f"{case.name} {op}: W(-), W(+) = {weights[0]:.12g}, {weights[1]:.12g}")
    return weights


def check_field(case, field, components):
    # Each spin sees the free dimer [[v, xi0], [xi0, 0]]. The other two
    # channels are spectators. Differentiate its energies w.r.t. v to get
    # the probabilities on the field site, and fill the modes at FDM_T.
    def modes(v):
        gap = math.hypot(v, 2 * XI0)
        for sign in (-1, 1):
            e = (v + sign * gap) / 2
            yield e, (1 + sign * v / gap) / 2, 1 / (1 + math.exp(e / FDM_T))

    up, down = list(modes(-field / 2)), list(modes(field / 2))
    mag = (sum(p * f for _, p, f in up) - sum(p * f for _, p, f in down)) / 2
    temp, measured = map(float, (case / "customfdm").read_text().splitlines()[-1].split())
    close("FDM temperature", temp, FDM_T)
    close("local <Sz>", measured, mag)
    print(f"b={field:g}: <Sz>={measured:.12g} (free dimer {mag:.12g})")
    for op in components:
        initial, final = (down, up) if op == "SpLoc" else (up, down)
        expected = [0., 0.]
        for ei, pi, fi in initial:
            for ej, pj, fj in final:
                expected[ej > ei] += pi * pj * (fi - fj)
        weights = signed_weights(case, "spin", op)
        for value, ref in zip(weights, expected):
            close(f"b={field:g} {op} signed weight", value, ref)
        # NRG autocorrelates O^dagger with O: SpLoc integrates to -2<Sz>.
        close(f"{op} commutator", sum(weights), (-2 if op == "SpLoc" else 2) * measured)


def run_case(case, symmetry, package, kernel, nrg, field=0):
    case.mkdir(parents=True)
    if symmetry == "QSTZ":
        operators = "ops=sigma_f\nspect=sigma_f-sigma_f\nspecchit=sigma_f-sigma_f\n"
    else:
        # Request only the lowering component in the reversed field.
        components = ("SpLoc", "SmLoc") if field >= 0 else ("SmLoc",)
        pairs = " ".join(f"{op}-{op}" for op in components)
        operators = f"spect={pairs}\nops=SZf0 {' '.join(components)}"
        if field == 0:
            operators += f" A0\nspecs=SZf0-SZf0\nspecchit=SZf0-SZf0 {pairs}\nspecgt=A0-A0\ngtp=1\n"
        else:
            operators += "\nfdmexpv=true\nprec_custom=16\n"
        (case / "modeloperators.m").write_text('''Join[
  mtQSZTZLadderOps["SpLoc", "SmLoc", spinplus[f[0]]],
  mtDoubletOp["A0", If[#2 == 0 && #3 == UP, f[#1, 0, UP], 0]&]]
''')
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
prec_xy=16
chitp=1
T={FDM_T}
fdm=true
savebins=true
broaden=false
{operators}""")
    (case / "model.m").write_text(f"""def3ch[0];
{symmetry.lower()}ImpurityCharges = {{}};
H = H0 /. f[ca_, ch_, sp_] :> f[ca, 1 - ch, sp];
H = H - ({field}) spinz[f[0]];
""")
    # Run the installed initializer directly, pinning PACKAGEPATH instead of
    # the launcher's ~/ and ~/nrg overrides. basis.m still does real module
    # dispatch; neither the fresh case nor its parent contains a basis file.
    bootstrap = '''NRGDIR = Environment["NRG_INIT_DIR"];
PACKAGEPATH = {NRGDIR};
SYMTYPE = "runtime";
Get["sneg.m", Path -> PACKAGEPATH];
Get["initial.m", Path -> PACKAGEPATH];
'''

    def run(command, status=0, **kwargs):
        result = subprocess.run(command, cwd=case, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                text=True, **kwargs)
        if result.returncode != status:
            raise RuntimeError(f"{case.name}: {command[0]} exited {result.returncode}, expected {status}:\n{result.stdout}")
        return result.stdout

    def initialize(body, status=0):
        return run([kernel, "-noinit", "-noprompt", "-batchinput", "-batchoutput"], status,
                   input=bootstrap + body + "\nExit[0];\n", env=dict(os.environ, NRG_INIT_DIR=str(package)))

    log = initialize('makedata["data"]; Print["INITIALIZER COMPLETE"];')
    require("INITIALIZER COMPLETE" in log and (case / "data").is_file() and
            f"Loading module {symmetry}basis.m" in log, f"incomplete initializer: {log}")
    if symmetry == "QSZTZ":
        headers = re.findall(r"^t(SpLoc|SmLoc)$", (case / "data").read_text(), re.M)
        require(sorted(headers) == sorted(components), f"requested {components}, generated {headers}")
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
    if field:
        check_field(case, field, components)
    else:
        # The central orbital has E=0 in the seed, then equal weights at +/-xi0
        # in the dimer. GT is -f'(E) for one spin; T*chi_z also has interlevel
        # matrix elements of the local Sz. Output order is increasing T (step, seed).
        for op in (("sigma_f",) if symmetry == "QSTZ" else ("SZf0", *components, "A0")):
            prefix = "gt_GT" if op == "A0" else "chit_CHIT"
            filename = f"{prefix}_dens_{op}-{op}.dat"
            curve = (case / filename).read_text().splitlines()
            require(len(curve) == 2, f"{filename}: expected seed and one extension")
            for line, td, energy in zip(curve, reversed(rows), (XI0, 0)):
                temp, value = map(float, line.split())
                require(abs(temp - td["T"]) < 1e-10, f"{filename}: wrong temperature {temp}")
                f = 1 / (1 + math.exp(energy / temp))
                if op == "A0":
                    expected = f * (1 - f) / temp
                else:
                    expected = f * (1 - f) / 4 + temp * (1 - 2 * f) / (8 * energy) if energy else 1 / 8
                    if op in ("SpLoc", "SmLoc"):
                        expected *= 2
                require(abs(value - expected) < 1e-10, f"{filename} at T={temp}: {value} != {expected}")
                print(f"{filename} at T={temp:.6g}: {value:.12g} (absolute error {abs(value - expected):.3g})")
        # Request an active name so the guard cannot pass via operator filtering.
        names = '"sigma_f", "SmLoc"' if symmetry == "QSTZ" else '"SpLoc", "SpLoc"'
        error = initialize(f'mtQSZTZLadderOps[{names}, spinplus[f[0]]];', status=1)
        require("mtQSZTZLadderOps" in error and "Aborting." in error, error)
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
        runs = root / "runs"
        reference, additive = [run_case(runs / sym, sym, package, kernel, nrg)
                               for sym in ("QSTZ", "QSZTZ")]
        for i, (ref, row) in enumerate(zip(reference, additive)):
            errors = {k: abs(row[k] - ref[k]) for k in ("T", "<Sz^2>", "<Tz^2>", "<Q>", "<Q^2>",
                                                         "<E>", "<E^2>", "C", "F", "S")}
            errors.update({k: abs(row[k]) for k in ("<Sz>", "<Tz>")})
            require(max(errors.values()) < 1e-10, f"TD row {i} absolute differences: {errors}")
            print(f"TD row {i}: max absolute difference {max(errors.values()):.3g}")
        require(abs(additive[0]["S"] - math.log(64)) < 1e-10, "degenerate seed entropy != ln(64)")

        # At zero field the interlevel transition has local matrix element
        # 1/2 per spin flip; the Cartesian component carries half its weight.
        w = math.tanh(XI0 / (2 * FDM_T)) / 4
        for sym, prefix, op, factor in (("QSTZ", "spin", "sigma_f", 0.5), ("QSZTZ", "corr", "SZf0", 0.5),
                                        ("QSZTZ", "spin", "SpLoc", 1), ("QSZTZ", "spin", "SmLoc", 1)):
            neg, pos = signed_weights(runs / sym, prefix, op)
            close(op + " negative weight", neg, -factor * w)
            close(op + " positive weight", pos, factor * w)
        for field in (0.4, -0.4):
            run_case(runs / f"field{field}", "QSZTZ", package, kernel, nrg, field)


if __name__ == "__main__":
    main()
