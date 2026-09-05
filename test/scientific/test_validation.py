"""Focused unittest checks for scientific validation parsing and failure contracts."""

from collections import Counter
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
import json
from math import comb
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

from ed_siam import AndersonModel, flat_band_model, solve
import validate_siam as validation


QSZ_FIXTURE = Path(__file__).resolve().parent / "fixtures" / "siam_qsz"


class SpectrumTests(unittest.TestCase):
    def test_qsz_dimensions_match_fock_counts_and_spin_labels(self):
        for orbitals in range(1, 7):
            with self.subTest(orbitals=orbitals):
                expected = Counter(
                    (up.bit_count() + down.bit_count() - orbitals,
                     up.bit_count() - down.bit_count() + 1)
                    for up in range(1 << orbitals) for down in range(1 << orbitals)
                )
                counts = validation.dimensions(orbitals, "QSZ")
                self.assertEqual(counts, expected)
                self.assertEqual(sum(counts.values()), 4 ** orbitals)
                self.assertEqual(counts[0, 1 - orbitals], 1)
                self.assertEqual(counts[0, 1 + orbitals], 1)
                for up in range(orbitals + 1):
                    for down in range(orbitals + 1):
                        self.assertEqual(counts[up + down - orbitals, up - down + 1],
                                         comb(orbitals, up) * comb(orbitals, down))

    def test_qs_dimensions_count_multiplets_not_spin_states(self):
        self.assertEqual(validation.dimensions(2, "QS"), {
            (-2, 1): 1, (-1, 2): 2, (0, 1): 3, (0, 3): 1, (1, 2): 2, (2, 1): 1,
        })
        for orbitals in range(1, 7):
            with self.subTest(orbitals=orbitals):
                counts = validation.dimensions(orbitals, "QS")
                expanded = Counter()
                for (charge, multiplicity), count in counts.items():
                    self.assertGreater(count, 0)
                    for twice_sz in range(1 - multiplicity, multiplicity, 2):
                        expanded[charge, twice_sz + 1] += count
                self.assertEqual(expanded, validation.dimensions(orbitals, "QSZ"))
                self.assertEqual(sum(spin * count for (_, spin), count in counts.items()),
                                 4 ** orbitals)
                for charge in range(-orbitals, orbitals + 1):
                    self.assertEqual(sum(spin * count for (q, spin), count in counts.items()
                                         if q == charge), comb(2 * orbitals, charge + orbitals))

    def test_qs_spectrum_expands_doublets_and_triplets(self):
        solution = solve(AndersonModel(-0.4, 1.0, 0.0, (0.0,), ()))
        sectors = {(-2, 1): [0.0], (-1, 2): [0.0, -0.4],
                   (0, 1): [0.2, 0.0, -0.4], (0, 3): [-0.4],
                   (1, 2): [0.2, -0.4], (2, 1): [0.2]}
        checks = []
        validation.compare_spectrum(sectors, solution, "QS", "multiplets", checks)
        self.assertEqual(len(checks), 19)
        self.assertEqual(sum(check["count"] for check in checks), 33)
        self.assertIn("multiplets/sector(0, -1)", [check["quantity"] for check in checks])
        self.assertIn("multiplets/sector(0, 3)", [check["quantity"] for check in checks])

    def test_qsz_spectrum_maps_up_minus_down_plus_one(self):
        solution = solve(AndersonModel(-0.3, 0.7, 0.0, (0.0,), (), B=0.2))
        sectors = {(up + down - 2, up - down + 1): energies[::-1]
                   for (up, down), energies in solution.sector_energies.items()}
        checks = []
        validation.compare_spectrum(sectors, solution, "QSZ", "spin", checks)
        self.assertEqual(sum(check["count"] for check in checks), 33)
        reflected = {(charge, 2 - spin): energies for (charge, spin), energies in sectors.items()}
        with self.assertRaisesRegex(ValueError, r"spin/sector.*actual=.*expected=.*error=.*limit="):
            validation.compare_spectrum(reflected, solution, "QSZ", "spin", [])

    def test_spectrum_rejects_missing_extra_and_wrong_dimension_sectors(self):
        solution = solve(AndersonModel(-0.3, 0.7, 0.2, (0.0,), ()))
        sectors = {(up + down - 2, up - down + 1): energies
                   for (up, down), energies in solution.sector_energies.items()}
        for invalid, diagnostic in (
            ({}, "missing or unexpected sectors"),
            ({key: values for key, values in sectors.items() if key != (-2, 1)},
             "missing or unexpected sectors"),
            ({**sectors, (99, 1): [0.0]}, "missing or unexpected sectors"),
            ({**sectors, (0, 1): sectors[0, 1][:-1]}, r"wrong dimension for \(0, 1\)"),
        ):
            with self.subTest(diagnostic=diagnostic, sectors=list(invalid)):
                with self.assertRaisesRegex(ValueError, "spectrum: " + diagnostic):
                    validation.compare_spectrum(invalid, solution, "QSZ", "spectrum", [])


class ReportParsingTests(unittest.TestCase):
    def test_repeated_iteration_labels_preserve_distinct_blocks(self):
        blocks = validation.parse_report("""
            === Report N=0
            Sector I=-1 0
            I=-1 0 n=0 E=-1.25 I=1 n_d=0.5 n_d_ud=1e-3
            I=-1 0 n=1 E=2.5e-2

            === Report N=0
            Sector I=-1 0
            I=-1 0 n=0 E=-2
            Sector I=0 -1
            I=0 -1 n=0 E=0
            === Report N=1
            Sector I=1 2
            I=1 2 n=0 E=3
        """)
        self.assertEqual(blocks, [(0, {(-1, 0): [-1.25, 0.025]}),
                                  (0, {(-1, 0): [-2.0], (0, -1): [0.0]}),
                                  (1, {(1, 2): [3.0]})])
        self.assertIsNot(blocks[0][1], blocks[1][1])

    def test_report_rejects_empty_blocks(self):
        valid = "=== Report N=0\nSector I=0 1\nI=0 1 n=0 E=0\n"
        for text in ("", " \n", "=== Report N=0", "=== Report N=0\nSector I=0 1",
                     "=== Report N=0\n" + valid, valid + "Sector I=1 2\n"):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, "empty spectrum report"):
                validation.parse_report(text)

    def test_report_rejects_duplicate_or_misplaced_sectors_and_states(self):
        prefix = "=== Report N=0\nSector I=0 1\nI=0 1 n=0 E=0\n"
        for text, diagnostic in (
            (prefix + "Sector I=0 1", "duplicate or misplaced report sector"),
            ("Sector I=0 1", "duplicate or misplaced report sector"),
            ("I=0 1 n=0 E=0", "misplaced report state"),
            ("=== Report N=0\nI=0 1 n=0 E=0", "misplaced report state"),
            (prefix + "I=1 2 n=1 E=0", "misplaced report state"),
        ):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, diagnostic):
                validation.parse_report(text)

    def test_report_requires_contiguous_zero_based_state_indices(self):
        prefix = "=== Report N=0\nSector I=0 1\n"
        for states in ("I=0 1 n=1 E=0", "I=0 1 n=0 E=0\nI=0 1 n=0 E=1",
                       "I=0 1 n=0 E=0\nI=0 1 n=2 E=1"):
            with self.subTest(states=states):
                with self.assertRaisesRegex(ValueError, "noncontiguous report state indices"):
                    validation.parse_report(prefix + states)

    def test_report_rejects_malformed_lines(self):
        prefix = "=== Report N=0\nSector I=0 1\n"
        for line in ("=== Report N=-1", "Sector I=0", "Sector I=0 1 2", "Sector I=x 1",
                     "I=0 1 n=-1 E=0", "I=0 1 n=0", "unrecognized output"):
            with self.subTest(line=line), self.assertRaisesRegex(ValueError, "invalid report line:"):
                validation.parse_report(prefix + line)

    def test_report_rejects_nonfinite_energies_and_observables(self):
        prefix = "=== Report N=0\nSector I=0 1\nI=0 1 n=0 "
        for value in ("nan", "inf", "-inf", "1e400"):
            for fields in (f"E={value}", f"E=0 I=1 n_d={value}"):
                with self.subTest(fields=fields):
                    with self.assertRaisesRegex(ValueError, "nonfinite numerical output"):
                        validation.parse_report(prefix + fields)


class SubspaceParsingTests(unittest.TestCase):
    def test_subspaces_preserve_iteration_order_and_kept_total_counts(self):
        self.assertEqual(validation.parse_subspaces("""
            Iteration 0
            len_dm=2
            I=-1 0 kept=1 total=2
            I=0 -1 kept=1 total=1

            Iteration 1
            len_dm=1
            I=1 2 kept=4 total=4
        """), [(0, {(-1, 0): (1, 2), (0, -1): (1, 1)}), (1, {(1, 2): (4, 4)})])

    def test_subspaces_reject_empty_and_truncated_blocks(self):
        valid = "Iteration 1\nlen_dm=1\nI=0 1 kept=1 total=1\n"
        for text in ("", " \n", "Iteration 0", "Iteration 0\nlen_dm=0",
                     "Iteration 0\nlen_dm=2\nI=0 1 kept=1 total=1",
                     "Iteration 0\nlen_dm=2\nI=0 1 kept=1 total=1\n" + valid,
                     "Iteration 0\nlen_dm=0\n" + valid):
            with self.subTest(text=text):
                with self.assertRaisesRegex(ValueError, "empty|truncated|incorrect subspace count"):
                    validation.parse_subspaces(text)

    def test_subspaces_reject_duplicate_or_misplaced_records(self):
        for text, diagnostic in (
            ("len_dm=1", "misplaced subspace count"),
            ("Iteration 0\nlen_dm=1\nlen_dm=1", "misplaced subspace count"),
            ("I=0 1 kept=1 total=1", "duplicate/misplaced subspace"),
            ("Iteration 0\nI=0 1 kept=1 total=1", "duplicate/misplaced subspace"),
            ("Iteration 0\nlen_dm=2\nI=0 1 kept=1 total=1\nI=0 1 kept=1 total=1",
             "duplicate/misplaced subspace"),
        ):
            with self.subTest(text=text), self.assertRaisesRegex(ValueError, diagnostic):
                validation.parse_subspaces(text)

    def test_subspaces_reject_malformed_lines(self):
        for line in ("Iteration -1", "len_dm=-1", "I=0 kept=1 total=1",
                     "I=0 1 kept=-1 total=1", "I=0 1 kept=1 total=nan",
                     "I=0 1 kept=1 total=1 extra", "unrecognized output"):
            with self.subTest(line=line):
                with self.assertRaisesRegex(ValueError, "invalid subspaces line:"):
                    validation.parse_subspaces("Iteration 0\nlen_dm=1\n" + line)


class TableTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / "table"

    def test_table_accepts_one_or_two_headers_and_one_finite_row(self):
        for headers in (1, 2):
            with self.subTest(headers=headers):
                self.path.write_text(("# expectation values\n" if headers == 2 else "")
                                     + "# T I n_d\n\n0.05 1 1.25e+0\n\n")
                self.assertEqual(validation.read_table(self.path, ["T", "I", "n_d"], headers),
                                 {"T": 0.05, "I": 1.0, "n_d": 1.25})

    def test_table_rejects_missing_headers_or_wrong_row_count(self):
        for text in ("", "# T I\n", "T I\n0.05 1\n", "# T I\n0.05 1\n0.2 1\n"):
            with self.subTest(text=text):
                self.path.write_text(text)
                with self.assertRaisesRegex(ValueError, "expected headers and exactly one data row") as error:
                    validation.read_table(self.path, ["T", "I"], 1)
                self.assertIn(str(self.path), str(error.exception))
        self.path.write_text("not a comment\n# T I\n0.05 1\n")
        with self.assertRaisesRegex(ValueError, "expected headers and exactly one data row"):
            validation.read_table(self.path, ["T", "I"], 2)

    def test_table_rejects_wrong_columns_and_column_counts(self):
        for text, diagnostic in (("# I T\n1 0.05\n", "unexpected columns"),
                                 ("# T other\n0.05 1\n", "unexpected columns"),
                                 ("# T I\n0.05\n", "incorrect column count"),
                                 ("# T I\n0.05 1 2\n", "incorrect column count")):
            with self.subTest(text=text):
                self.path.write_text(text)
                with self.assertRaisesRegex(ValueError, diagnostic) as error:
                    validation.read_table(self.path, ["T", "I"], 1)
                self.assertIn(str(self.path), str(error.exception))

    def test_table_rejects_nonfinite_and_nonnumeric_values(self):
        for value, diagnostic in (("nan", "nonfinite numerical output"),
                                  ("inf", "nonfinite numerical output"),
                                  ("-inf", "nonfinite numerical output"),
                                  ("1e400", "nonfinite numerical output"),
                                  ("bad", "could not convert string to float")):
            with self.subTest(value=value):
                self.path.write_text(f"# T I\n0.05 {value}\n")
                with self.assertRaisesRegex(ValueError, diagnostic):
                    validation.read_table(self.path, ["T", "I"], 1)


class ComparisonTests(unittest.TestCase):
    def test_compare_accepts_scalars_and_values_within_tolerance(self):
        checks = []
        validation.compare("scalar", 0.25, 0.25, checks)
        validation.compare("default", [1.0 + 5e-11], [1.0], checks)
        validation.compare("tolerant", [0.0005, 101.0], [0.0, 100.0], checks, atol=1e-3, rtol=1e-2)
        validation.compare("boundary", 0.5, 0.0, checks, atol=0.5, rtol=0)
        self.assertEqual([check["count"] for check in checks], [1, 1, 2, 1])
        self.assertEqual(checks[0], {"quantity": "scalar", "count": 1, "max_absolute_error": 0.0})
        self.assertEqual(checks[2], {"quantity": "tolerant", "count": 2, "max_absolute_error": 1.0})

    def test_compare_rejects_empty_or_broadcastable_mismatched_shapes(self):
        for actual, expected in (([], []), ([[]], [[]]), ([], [1]), ([1], [[1]]), ([1, 2], [1])):
            with self.subTest(actual=actual, expected=expected):
                checks = []
                with self.assertRaisesRegex(ValueError, "shape: shape mismatch or empty comparison:"):
                    validation.compare("shape", actual, expected, checks)
                self.assertEqual(checks, [])

    def test_compare_rejects_nonfinite_actual_or_expected(self):
        for value in (float("nan"), float("inf"), -float("inf")):
            for actual, expected in (([value], [0.0]), ([0.0], [value])):
                with self.subTest(actual=actual, expected=expected):
                    checks = []
                    with self.assertRaisesRegex(ValueError, "observable: nonfinite values"):
                        validation.compare("observable", actual, expected, checks)
                    self.assertEqual(checks, [])

    def test_compare_mismatch_identifies_quantity_values_error_and_limit(self):
        checks = [{"quantity": "earlier", "count": 1, "max_absolute_error": 0.0}]
        before = checks.copy()
        with self.assertRaisesRegex(ValueError, r"density at \(") as error:
            validation.compare("density", [[0, 2], [30, 0]], [[0, 0], [0, 0]], checks, atol=0.1, rtol=0)
        for diagnostic in ("actual=30", "expected=0", "error=30", "limit=0.1"):
            self.assertIn(diagnostic, str(error.exception))
        self.assertEqual(checks, before)


class InspectDataTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.path = Path(temporary.name) / "data"
        shutil.copyfile(QSZ_FIXTURE / "data", self.path)
        self.text = self.path.read_text()
        self.case = validation.load_case(QSZ_FIXTURE)

    def test_qsz_fixture_has_physical_seed_and_chain_coefficients(self):
        checks = []
        validation.inspect_data(self.path, self.case, checks)
        self.assertEqual([check["quantity"] for check in checks[:3]],
                         ["chain/hoppings", "chain/onsites", "seed_input/ground"])
        self.assertEqual(sum(check["count"] for check in checks), 2 * self.case["model"]["bath_sites"] + 33)
        self.assertLess(max(check["max_absolute_error"] for check in checks), validation.ATOL)
        self.assertEqual(self.path.read_text(), self.text)

    def test_inspect_rejects_wrong_format_symmetry_scaling_and_boundary(self):
        for old, new, diagnostic in (
            ("#!9", "#!8", "fixture must have format #!9"),
            ("# symtype  QSZ", "# symtype  QS", "wrong seed symmetry"),
            ("# SCALE  1\n", "# SCALE  2\n", "fixture seed must be in physical units"),
            ("#!9", "#!9\n# COMPLEX", "complex fixture format is not supported"),
            ("1 2 9\n", "2 2 9\n", "wrong fixture chain length"),
            ("1 2 9\n", "1 1 9\n", "wrong fixture chain length"),
            ("f 0 0", "g 0 0", "missing seed boundary operators"),
        ):
            with self.subTest(new=new):
                self.assertIn(old, self.text)
                self.path.write_text(self.text.replace(old, new, 1))
                with self.assertRaisesRegex(ValueError, diagnostic):
                    validation.inspect_data(self.path, self.case, [])

    def test_inspect_rejects_invalid_seed_sectors_counts_and_energies(self):
        lines = self.text.splitlines()
        first = lines.index("-2 1")
        ground = lines.index("e") + 1
        for index, value, diagnostic in (
            (first + 3, "-2 1", "invalid seed sector"),
            (first, "-2 1 0", "invalid seed sector"),
            (first + 1, "0", "invalid seed sector"),
            (first + 1, "2", "invalid seed sector"),
            (first, "-3 1", "seed_input: missing or unexpected sectors"),
            (first + 2, "nan", "nonfinite numerical output"),
            (first + 2, str(float(lines[first + 2]) + 0.1), r"seed_input/sector.*actual=.*expected="),
            (ground, "inf", "nonfinite numerical output"),
            (ground, str(float(lines[ground]) + 0.1), r"seed_input/ground.*actual=.*expected="),
        ):
            with self.subTest(index=index, value=value):
                changed = lines.copy()
                changed[index] = value
                self.path.write_text("\n".join(changed) + "\n")
                with self.assertRaisesRegex(ValueError, diagnostic):
                    validation.inspect_data(self.path, self.case, [])

    def test_inspect_rejects_rescaled_seed_even_with_physical_scale_header(self):
        lines = self.text.splitlines()
        first = lines.index("-2 1")
        for sector in range(9):
            index = first + 3 * sector + 2
            lines[index] = " ".join(str(2 * float(value)) for value in lines[index].split())
        self.path.write_text("\n".join(lines) + "\n")
        with self.assertRaisesRegex(ValueError, r"seed_input/sector.*actual=.*expected=.*error="):
            validation.inspect_data(self.path, self.case, [])

    def test_inspect_rejects_missing_or_repeated_trailer_markers(self):
        for marker in ("e", "z"):
            for replacement in ("\n", f"\n{marker}\n{marker}\n"):
                with self.subTest(marker=marker, replacement=replacement):
                    self.path.write_text(self.text.replace(f"\n{marker}\n", replacement, 1))
                    with self.assertRaisesRegex(ValueError, "missing or repeated e/z fixture markers"):
                        validation.inspect_data(self.path, self.case, [])

    def test_inspect_rejects_wrong_xi_onsites_lengths_and_nonfinite_coefficients(self):
        seed, trailer = self.text.split("\nz\n")
        coefficients = trailer.splitlines()
        nmax = int(coefficients[0])
        for index, value, diagnostic in (
            (0, str(nmax - 1), "wrong coefficient maximum index"),
            (nmax + 2, str(nmax + 1), "wrong coefficient maximum index"),
            (1, str(float(coefficients[1]) + 0.01), r"chain/hoppings.*actual=.*expected="),
            (nmax + 1, str(float(coefficients[nmax + 1]) + 0.01), r"chain/hoppings.*actual=.*expected="),
            (nmax + 3, "0.01", r"chain/onsites.*actual=.*expected="),
            (1, "nan", "nonfinite numerical output"),
            (nmax + 3, "-inf", "nonfinite numerical output"),
        ):
            with self.subTest(index=index, value=value):
                changed = coefficients.copy()
                changed[index] = value
                self.path.write_text(seed + "\nz\n" + "\n".join(changed) + "\n")
                with self.assertRaisesRegex(ValueError, diagnostic):
                    validation.inspect_data(self.path, self.case, [])
        for changed in (coefficients[:-1], coefficients + ["0"]):
            with self.subTest(length=len(changed)):
                self.path.write_text(seed + "\nz\n" + "\n".join(changed) + "\n")
                with self.assertRaisesRegex(ValueError, "unexpected payload after legacy chain trailer"):
                    validation.inspect_data(self.path, self.case, [])


class RunNrgTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.directory = self.root / "run"
        self.directory.mkdir()
        self.executable = self.root / "fake-nrg"

    def test_nonzero_exit_fails_even_when_done_exists(self):
        self.executable.write_text(f"#!{sys.executable}\nfrom pathlib import Path\nimport sys\n"
                                   "Path('DONE').write_text('not success')\n"
                                   "print('intentional NRG failure', flush=True)\nsys.exit(7)\n")
        self.executable.chmod(0o755)
        with self.assertRaises(subprocess.CalledProcessError) as error:
            validation.run_nrg(self.executable, self.directory)
        self.assertEqual(error.exception.returncode, 7)
        self.assertEqual(error.exception.cmd, [str(self.executable), "-w", str(self.directory / "scratch")])
        self.assertTrue((self.directory / "DONE").is_file())
        self.assertIn("intentional NRG failure", (self.directory / "log").read_text())

    def test_zero_exit_without_done_is_not_success(self):
        self.executable.write_text(f"#!{sys.executable}\nprint('exited without marker')\n")
        self.executable.chmod(0o755)
        with self.assertRaisesRegex(ValueError, "NRG exited without a fresh DONE"):
            validation.run_nrg(self.executable, self.directory)
        self.assertIn("exited without marker", (self.directory / "log").read_text())

    def test_timeout_is_propagated(self):
        self.executable.write_text(f"#!{sys.executable}\nimport time\ntime.sleep(10)\n")
        self.executable.chmod(0o755)
        with self.assertRaises(subprocess.TimeoutExpired) as error:
            validation.run_nrg(self.executable, self.directory, timeout=0.05)
        self.assertGreaterEqual(error.exception.timeout, 0)
        self.assertLessEqual(error.exception.timeout, 0.05)
        self.assertIn("timed out", str(error.exception))
        self.assertEqual(error.exception.cmd, [str(self.executable), "-w", str(self.directory / "scratch")])
        self.assertFalse((self.directory / "DONE").exists())

    def test_success_requires_regular_done_and_isolated_single_thread_environment(self):
        self.executable.write_text(f"#!{sys.executable}\nfrom pathlib import Path\nimport json, os, sys\n"
                                   "Path('DONE').write_text('fresh marker')\n"
                                   "print(json.dumps({'args': sys.argv[1:], 'cwd': os.getcwd(), "
                                   "'env': dict(os.environ)}))\n")
        self.executable.chmod(0o755)
        threads = ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS", "BLIS_NUM_THREADS",
                   "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS")
        with mock.patch.dict(os.environ, {name: "8" for name in threads}):
            validation.run_nrg(self.executable, self.directory)
        result = json.loads((self.directory / "log").read_text())
        scratch = self.directory / "scratch"
        self.assertTrue(scratch.is_dir())
        self.assertEqual(result["args"], ["-w", str(scratch)])
        self.assertEqual(result["cwd"], str(self.directory))
        for name in threads:
            self.assertEqual(result["env"][name], "1")
        for name in ("TMPDIR", "TMP", "TEMP"):
            self.assertEqual(result["env"][name], str(scratch))
        for name in ("OMP_DYNAMIC", "MKL_DYNAMIC"):
            self.assertEqual(result["env"][name], "FALSE")
        self.assertEqual(result["env"]["LC_ALL"], "C")
        self.assertEqual((self.directory / "DONE").read_text(), "fresh marker")
        self.assertFalse((self.directory / "DONE").is_symlink())

    def test_symlink_done_is_not_success(self):
        self.executable.write_text(f"#!{sys.executable}\nfrom pathlib import Path\n"
                                   "Path('target').write_text('old marker')\n"
                                   "Path('DONE').symlink_to('target')\n")
        self.executable.chmod(0o755)
        with self.assertRaisesRegex(ValueError, "NRG exited without a fresh DONE"):
            validation.run_nrg(self.executable, self.directory)


class RunCaseTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fixture = self.root / "fixture"
        self.fixture.mkdir()
        for name in ("model.json", "data", "provenance.json"):
            shutil.copyfile(QSZ_FIXTURE / name, self.fixture / name)
        self.work_root = self.root / "runs"

    def test_repeated_failed_runs_are_fresh_and_leave_source_unchanged(self):
        for name in ("DONE", "param", "report.nrg", "validation.json"):
            (self.fixture / name).write_text(f"stale source {name}\n")
        before = {path.name: path.read_bytes() for path in self.fixture.iterdir()}
        directories = []

        def fail_run(executable, directory):
            self.assertEqual(executable, Path(sys.executable).resolve())
            self.assertEqual({path.name for path in directory.iterdir()}, {"data", "param"})
            self.assertEqual((directory / "data").read_bytes(), before["data"])
            self.assertIn("data_has_rescaled_energies=false\n", (directory / "param").read_text())
            directories.append(directory)
            (directory / "DONE").write_text("partial run")
            raise RuntimeError("injected NRG failure")

        with mock.patch.object(validation, "run_nrg", side_effect=fail_run) as run, redirect_stdout(StringIO()):
            for _ in range(2):
                with self.assertRaisesRegex(RuntimeError, "injected NRG failure"):
                    validation.run_case(sys.executable, self.fixture, 0.05, "rescaled", self.work_root)
        self.assertEqual(run.call_count, 2)
        self.assertEqual(len(set(directories)), 2)
        self.assertEqual(set(self.work_root.iterdir()), set(directories))
        for directory in directories:
            report = json.loads((directory / "validation.json").read_text())
            self.assertEqual(report["status"], "failed")
            self.assertEqual(report["error"], "injected NRG failure")
            self.assertEqual(report["temperature"], 0.05)
            self.assertEqual(report["mode"], "rescaled")
            self.assertGreater(sum(check["count"] for check in report["checks"]), 0)
        self.assertEqual({path.name: path.read_bytes() for path in self.fixture.iterdir()}, before)

    def test_run_case_rejects_unrunnable_executable_before_launch(self):
        nonexecutable = self.root / "not-executable"
        nonexecutable.write_text("not a program")
        nonexecutable.chmod(0o600)
        for executable in (self.root / "missing", nonexecutable):
            with self.subTest(executable=executable), mock.patch.object(validation, "run_nrg") as run:
                with self.assertRaisesRegex(ValueError, "NRG executable is not runnable:"):
                    validation.run_case(executable, self.fixture, 0.05, "rescaled", self.work_root)
                run.assert_not_called()
                self.assertFalse(self.work_root.exists())

    def test_run_case_rejects_changed_model_provenance_before_launch(self):
        case = json.loads((self.fixture / "model.json").read_text())
        case["model"]["Gamma"] += 0.01
        (self.fixture / "model.json").write_text(json.dumps(case))
        with mock.patch.object(validation, "run_nrg") as run:
            with self.assertRaisesRegex(ValueError, "model specification changed: regenerate the fixture"):
                validation.run_case(sys.executable, self.fixture, 0.05, "rescaled", self.work_root)
            run.assert_not_called()
        self.assertFalse(self.work_root.exists())

    def test_run_case_rejects_changed_data_checksum_before_launch(self):
        data = self.fixture / "data"
        data.write_bytes(data.read_bytes() + b"\n")
        with mock.patch.object(validation, "run_nrg") as run:
            with self.assertRaisesRegex(ValueError, "fixture data checksum mismatch"):
                validation.run_case(sys.executable, self.fixture, 0.05, "rescaled", self.work_root)
            run.assert_not_called()
        self.assertFalse(self.work_root.exists())

    def test_run_case_rejects_truncation_counts_and_wrong_iteration_sequences(self):
        model = flat_band_model(validation.load_case(self.fixture)["model"])
        reports, subspaces = [], []
        for prefix in range(1, len(model.zeta) + 1):
            solution = solve(replace(model, zeta=model.zeta[:prefix], t=model.t[:prefix - 1]))
            report = [f"=== Report N={max(0, prefix - 2)}"]
            retained = [f"Iteration {prefix - 2}", f"len_dm={len(solution.sector_energies)}"]
            for (up, down), energies in solution.sector_energies.items():
                sector = f"{up + down - prefix - 1} {up - down + 1}"
                report.append(f"Sector I={sector}")
                report.extend(f"I={sector} n={index} E={energy:.17g}" for index, energy in enumerate(energies))
                retained.append(f"I={sector} kept={len(energies)} total={len(energies)}")
            reports.append("\n".join(report) + "\n")
            if prefix > 1:
                subspaces.append("\n".join(retained) + "\n")
        report_text, subspace_text = "".join(reports), "".join(subspaces)
        for spectrum, retained, diagnostic in (
            (report_text, subspace_text.replace("kept=3 total=3", "kept=2 total=3", 1),
             r"prefix2 sector.*truncation or missing eigenpairs: \(2, 3\), expected 3"),
            (report_text, subspace_text.replace("kept=3 total=3", "kept=3 total=4", 1),
             r"prefix2 sector.*truncation or missing eigenpairs: \(3, 4\), expected 3"),
            (reports[0] + reports[1].replace("N=0", "N=1", 1) + reports[2], subspace_text,
             "wrong prefix report sequence"),
            (report_text, subspace_text.replace("Iteration 0", "Iteration 1", 1),
             "wrong subspace iteration sequence"),
        ):
            def write_outputs(executable, directory):
                (directory / "report.nrg").write_text(spectrum)
                (directory / "subspaces.dat").write_text(retained)

            with self.subTest(diagnostic=diagnostic), redirect_stdout(StringIO()):
                with mock.patch.object(validation, "run_nrg", side_effect=write_outputs) as run:
                    with self.assertRaisesRegex(ValueError, diagnostic) as error:
                        validation.run_case(sys.executable, self.fixture, 0.05, "rescaled", self.work_root)
                directory = run.call_args.args[1]
                failure = json.loads((directory / "validation.json").read_text())
                self.assertEqual(failure["status"], "failed")
                self.assertEqual(failure["error"], str(error.exception))


if __name__ == "__main__":
    unittest.main()
