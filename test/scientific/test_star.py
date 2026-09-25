"""Independent star oracle checks and failure injection at the real suite entry point."""

from configparser import ConfigParser
from contextlib import redirect_stdout
from dataclasses import replace
from io import StringIO
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

from ed_siam import solve
from finite_star import FIXTURE, direct_delta, jacobi, load_star, reference_greens, reference_model, star_tables
import validate_star as validation


ENABLED_BACKENDS = tuple(filter(None, os.environ.get("NRG_TEST_CHAIN_BACKENDS", "legacy,rkpw").split(",")))
if set(ENABLED_BACKENDS) - {"legacy", "rkpw"} or len(set(ENABLED_BACKENDS)) != len(ENABLED_BACKENDS):
    raise ValueError("invalid NRG_TEST_CHAIN_BACKENDS")


class OracleTests(unittest.TestCase):
    def test_known_mass_mean_width_and_complete_spectral_measure(self):
        star = load_star()
        zeta, t = jacobi(star)
        self.assertAlmostEqual(np.linalg.norm(star["couplings"]), 0.3, places=15)
        self.assertAlmostEqual(zeta[0], -0.14, places=15)
        self.assertAlmostEqual(t[0], np.sqrt(0.1424), places=15)
        self.assertNotAlmostEqual(t[0], 0.3)
        bath = np.diag(zeta) + np.diag(t, 1) + np.diag(t, -1)
        poles, vectors = np.linalg.eigh(bath)
        order = np.argsort(star["energies"])
        np.testing.assert_allclose(poles, np.asarray(star["energies"])[order], atol=1e-15, rtol=1e-15)
        np.testing.assert_allclose(0.09 * vectors[0] ** 2,
                                   np.asarray(star["couplings"])[order] ** 2, atol=1e-15, rtol=1e-15)

    def test_direct_star_and_full_jacobi_resolvents_agree(self):
        star = load_star()
        s = np.asarray([0.02j, 0.3j, 2j, 0.12 + 0.4j, -0.8 + 0.1j])
        hstar = np.diag([star["epsilon_d"], *star["energies"]])
        hstar[0, 1:] = hstar[1:, 0] = star["couplings"]
        model = reference_model(star, 4)
        hchain = np.diag([model.epsilon_d, *model.zeta])
        hchain += np.diag([model.V, *model.t], 1) + np.diag([model.V, *model.t], -1)
        expected = reference_greens(star, s, 4)
        for h in (hstar, hchain):
            np.testing.assert_allclose([np.linalg.inv(z * np.eye(5) - h)[0, 0] for z in s],
                                       expected, atol=2e-14, rtol=2e-14)
        self.assertGreater(np.max(np.abs(expected.real)), 0.1)
        np.testing.assert_allclose(1 / expected, s - star["epsilon_d"] - direct_delta(star, s))

    def test_prefix_and_full_green_functions_match_many_body_lehmann(self):
        star = load_star()
        s = np.asarray([0.05j, 0.4j, 2j])
        for count in (2, 3, 4):
            solution = solve(reference_model(star, count))
            for temperature in (0.05, 0.2):
                with self.subTest(count=count, temperature=temperature):
                    expected = reference_greens(star, s, count)
                    for spin in solution.greens(temperature, s):
                        np.testing.assert_allclose(spin, expected, atol=2e-13, rtol=2e-13)

    def test_full_many_body_energies_include_centered_bath_constant(self):
        star = load_star()
        h = np.diag([star["epsilon_d"], *star["energies"]])
        h[0, 1:] = h[1:, 0] = star["couplings"]
        one_particle = np.repeat(np.linalg.eigvalsh(h), 2)
        occupations = (np.arange(1024)[:, None] >> np.arange(10)) & 1
        expected = occupations @ one_particle - sum(star["energies"])
        np.testing.assert_allclose(np.sort(solve(reference_model(star, 4)).energies),
                                   np.sort(expected), atol=3e-14, rtol=3e-14)

    def test_legacy_prefix_is_not_the_full_star_or_a_three_pole_subset(self):
        star = load_star()
        s = 1j * np.asarray([0.03, 0.15, 0.6])
        prefix = reference_greens(star, s, 3)
        full = reference_greens(star, s, 4)
        subset = {**star, "energies": star["energies"][:3], "couplings": star["couplings"][:3]}
        wrong = 1 / (s - star["epsilon_d"] - direct_delta(subset, s))
        for candidate in (full, wrong):
            with self.assertRaisesRegex(ValueError, "prefix G"):
                validation.compare("prefix G", candidate, prefix, [])

    def test_seed_spectrum_detects_hybridization_sign_and_unit_errors(self):
        seed = reference_model(load_star(), 2)
        expected = np.sort(solve(seed).energies)
        # A hopping sign alone is gauge-equivalent, so mutate ONSITE signs instead.
        for changed in (replace(seed, V=seed.V * np.sqrt(np.pi)),
                        replace(seed, V=seed.V * np.sqrt(2)),
                        replace(seed, zeta=tuple(-z for z in seed.zeta)),
                        replace(seed, zeta=tuple(z / 2 for z in seed.zeta), t=(seed.t[0] / 2,))):
            with self.assertRaisesRegex(ValueError, "seed"):
                validation.compare("seed", np.sort(solve(changed).energies), expected, [])


class AdapterTests(unittest.TestCase):
    def test_tables_and_mathematica_gamma_preserve_physical_units(self):
        star = load_star()
        for D in (1, 2):
            tables = star_tables(star, D)
            np.testing.assert_allclose(tables["de_pos.dat"] * D, [0.9, 0.4])
            np.testing.assert_allclose(tables["de_neg.dat"] * D, [0.7, 0.2])
            np.testing.assert_allclose(tables["du_pos.dat"], [0.2, 0.4])
            np.testing.assert_allclose(tables["du_neg.dat"], [0.4, 0.8])
            self.assertAlmostEqual(tables["theta.dat"][0], np.pi * 0.09)
            for backend in ("legacy", "rkpw"):
                config = ConfigParser()
                config.read_string(validation.parameters(star, backend, D, 0.2, "absolute"))
                gamma = config["extra"].getfloat("Gamma")
                self.assertAlmostEqual(np.sqrt(D * gamma * tables["theta.dat"][0] / np.pi), 0.3)
                p = config["param"]
                self.assertEqual(p["tri"], "old" if backend == "legacy" else "rkpw")
                self.assertEqual(p["tridiag_method"], "lanczos" if backend == "legacy" else "rkpw")
                self.assertEqual(p.getint("Ninit"), 1)
                self.assertEqual(p.getint("Nmax") + 1, validation.bath_count(backend))
                self.assertEqual(p.getint("keep"), 1024)
                self.assertFalse(p.getboolean("rescalexi"))
                self.assertFalse(p.getboolean("data_has_rescaled_energies"))
                self.assertTrue(p.getboolean("nrgchain_tables_load"))
                self.assertEqual(p["hook_pre_lanczosinit"], "import_star.m")

    def test_extra_hopping_and_terminal_zero_are_required(self):
        star = load_star()
        zeta, t = jacobi(star)
        validation.check_coefficients(t, zeta[:3], star, 3, [])
        validation.check_coefficients(np.r_[t, 0], zeta, star, 4, [])
        for xi, onsite, count in ((t[:2], zeta[:3], 3), (np.r_[t, 1e-16], zeta, 4),
                                  (t, -zeta[:3], 3), (t / 2, zeta[:3] / 2, 3)):
            with self.subTest(count=count), self.assertRaises(ValueError):
                validation.check_coefficients(xi, onsite, star, count, [])

    @unittest.skipUnless(ENABLED_BACKENDS, "no enabled chain backends")
    def test_prepared_assets_are_complete_and_physically_qualified(self):
        for backend in ENABLED_BACKENDS:
            star, asset = validation.verify_fixture(FIXTURE, backend, "runtime")
            checks = []
            validation.inspect_data(asset / "data", star, backend, 2, checks, runtime=True)
            self.assertGreater(sum(check["count"] for check in checks), 128)
            validation.verify_fixture(FIXTURE, backend, "instantiate")


class SuiteFailureContracts:
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        self.fixture = self.root / "fixture"
        self.fixture.mkdir()
        shutil.copyfile(FIXTURE / "star.json", self.fixture / "star.json")
        for name in (self.backend, "template"):
            shutil.copytree(FIXTURE / name, self.fixture / name)

    def run_case(self, lane="runtime"):
        return validation.run_case(Path(sys.executable), self.fixture, self.backend, lane, self.root / "runs",
                                   nrgchain=Path(sys.executable), instantiate=Path(sys.executable))

    def test_missing_data_and_template_asset_fail_before_launch(self):
        for lane, relative in (("runtime", f"{self.backend}/data"), ("instantiate", "template/ham_-3.1"),
                               ("instantiate", "template/op.A_d_3.1_2.2")):
            path = self.fixture / relative
            contents = path.read_bytes()
            path.unlink()
            with self.subTest(lane=lane, file=relative), mock.patch.object(validation, "run_nrg") as nrg:
                with mock.patch.object(validation, "run_tool") as tool:
                    with self.assertRaises(FileNotFoundError):
                        self.run_case(lane=lane)
                    nrg.assert_not_called()
                    tool.assert_not_called()
                    self.assertFalse((self.root / "runs").exists())
            path.write_bytes(contents)

    def test_tampered_provenance_fails_on_actual_suite_path(self):
        for field, value, diagnostic in (("backend", "incorrect", "backend provenance"),
                                         ("star_sha256", "bad", "star provenance"),
                                         ("star_file_sha256", "bad", "star file checksum"),
                                         ("files", {}, "incomplete fixture manifest"),
                                         ("kernel_version", "", "missing generation provenance")):
            path = self.fixture / self.backend / "provenance.json"
            original = path.read_text()
            provenance = json.loads(original)
            provenance[field] = value
            path.write_text(json.dumps(provenance))
            with self.subTest(field=field), mock.patch.object(validation, "run_nrg") as nrg:
                with self.assertRaisesRegex(ValueError, diagnostic):
                    self.run_case()
                nrg.assert_not_called()
            path.write_text(original)
        self.assertFalse((self.root / "runs").exists())

    def test_changed_data_or_operator_is_not_silently_accepted(self):
        for lane, name in (("runtime", f"{self.backend}/data"), ("instantiate", "template/op.A_d_3.1_2.2")):
            path = self.fixture / name
            path.write_text(path.read_text() + "\n")
            with self.subTest(lane=lane), mock.patch.object(validation, "run_tool") as tool:
                with self.assertRaisesRegex(ValueError, "fixture checksum mismatch"):
                    self.run_case(lane=lane)
                tool.assert_not_called()

    def test_undeclared_template_assets_are_not_copied_to_the_frontend(self):
        (self.fixture / "template/unexpected").write_text("not a qualified asset\n")
        with mock.patch.object(validation, "run_tool") as tool:
            with self.assertRaisesRegex(ValueError, "undeclared fixture asset"):
                self.run_case(lane="instantiate")
            tool.assert_not_called()

    def test_seed_rescaling_is_detected_even_with_updated_checksum_and_scale_one(self):
        path = self.fixture / self.backend / "data"
        lines = path.read_text().splitlines()
        first = lines.index("-3 1")
        for sector in range(10):
            index = first + 3 * sector + 2
            lines[index] = " ".join(str(2 * float(value)) for value in lines[index].split())
        path.write_text("\n".join(lines) + "\n")
        manifest = self.fixture / self.backend / "provenance.json"
        provenance = json.loads(manifest.read_text())
        provenance["files"]["data"] = validation.file_digest(path)
        manifest.write_text(json.dumps(provenance))
        with mock.patch.object(validation, "run_nrg") as nrg, redirect_stdout(StringIO()):
            with self.assertRaisesRegex(ValueError, "seed64/sector"):
                self.run_case()
            nrg.assert_not_called()
        reports = list((self.root / "runs").glob("*/run/validation.json"))
        self.assertEqual(len(reports), 1)
        self.assertEqual(json.loads(reports[0].read_text())["status"], "failed")

    def test_wrong_theta_sign_and_scaling_from_tool_fail_before_solver(self):
        star = load_star()
        zeta, t = jacobi(star)
        count = validation.bath_count(self.backend)
        for mutation, diagnostic in (("theta", "theta"), ("sign", "zeta"), ("scale", "xi")):
            def tool(command, directory, name):
                np.savetxt(directory / "xi.dat", np.r_[t, 0][:count] / (2 if mutation == "scale" else 1))
                np.savetxt(directory / "zeta.dat", (-zeta if mutation == "sign" else zeta)[:count])
                if mutation == "theta":
                    np.savetxt(directory / "theta.dat", [0.09])  # Missing pi.
            with self.subTest(mutation=mutation), redirect_stdout(StringIO()):
                with mock.patch.object(validation, "run_tool", side_effect=tool):
                    with mock.patch.object(validation, "run_nrg") as nrg:
                        with self.assertRaisesRegex(ValueError, diagnostic):
                            self.run_case(lane="nrgchain")
                        nrg.assert_not_called()

    def test_backends_do_not_consume_each_others_prepared_results(self):
        other = "legacy" if self.backend == "rkpw" else "rkpw"
        self.assertFalse((self.fixture / other).exists())
        with redirect_stdout(StringIO()), mock.patch.object(validation, "run_nrg", side_effect=RuntimeError("stop")) as nrg:
            with self.assertRaisesRegex(RuntimeError, "stop"):
                self.run_case()
            nrg.assert_called_once()

    def test_generation_parameter_tampering_fails_before_launch(self):
        mutations = [
            ("param", "tri", "rkpw"), ("param", "tridiag_method", "invalid"),
            ("param", "Ninit", "0"), ("param", "Nmax", "99"), ("param", "mMAX", "2"),
            ("param", "prec", "81"), ("param", "bandrescale", "1"),
            ("param", "data_has_rescaled_energies", "true"), ("param", "rescalexi", "true"),
            ("param", "nrgchain_tables_load", "false"), ("param", "options", "CHOP"),
            ("param", "hook_pre_lanczosinit", ""), ("extra", "Gamma", "1"),
            ("extra", "eps", "0.17"), ("extra", "U", "0.1"),
        ]
        for lane in ("runtime", "instantiate"):
            path = self.fixture / ("template" if lane == "instantiate" else self.backend) / "provenance.json"
            original = path.read_text()
            for section, name, value in mutations:
                # The template has no hook; an injected hook must also be rejected.
                if lane == "instantiate" and name == "hook_pre_lanczosinit":
                    value = "import_star.m"
                provenance = json.loads(original)
                config = ConfigParser()
                config.optionxform = str
                config.read_string(provenance["parameters"])
                if name == "tridiag_method":
                    value = "rkpw" if config[section][name] == "lanczos" else "lanczos"
                config[section][name] = value
                output = StringIO()
                config.write(output)
                provenance["parameters"] = output.getvalue()
                path.write_text(json.dumps(provenance))
                with self.subTest(lane=lane, name=name), mock.patch.object(validation, "run_tool") as tool:
                    with mock.patch.object(validation, "run_nrg") as nrg:
                        with self.assertRaisesRegex(ValueError, "generation"):
                            self.run_case(lane=lane)
                        tool.assert_not_called()
                        nrg.assert_not_called()
                path.write_text(original)
        self.assertFalse((self.root / "runs").exists())

    def test_missing_claimed_metadata_fails_before_launch(self):
        for lane in ("runtime", "instantiate"):
            path = self.fixture / ("template" if lane == "instantiate" else self.backend) / "provenance.json"
            original = path.read_text()
            for name in json.loads(original):
                provenance = json.loads(original)
                del provenance[name]
                path.write_text(json.dumps(provenance))
                with self.subTest(lane=lane, field=name), mock.patch.object(validation, "run_tool") as tool:
                    with mock.patch.object(validation, "run_nrg") as nrg:
                        with self.assertRaisesRegex(ValueError, "missing generation provenance"):
                            self.run_case(lane=lane)
                        tool.assert_not_called()
                        nrg.assert_not_called()
                path.write_text(original)
        self.assertFalse((self.root / "runs").exists())

    def test_incoherent_metadata_maps_and_hashes_are_rejected(self):
        for lane in ("runtime", "instantiate"):
            path = self.fixture / ("template" if lane == "instantiate" else self.backend) / "provenance.json"
            original = path.read_text()
            record = json.loads(original)
            altered_sources = {**record["initializer_source_files"], "nrginit/wilson.m": "0" * 64}
            altered_inputs = {**record["physical_input_files"], "theta.dat": "0" * 64}
            mutations = [("precision", 81), ("mMAX", 2), ("source_revision", "unknown"),
                         ("initializer_source_files", {}), ("initializer_source_files", altered_sources),
                         ("initializer_source_files", {"../initial.m": "0" * 64}),
                         ("initializer_source_files", {"nrginit/initial.m": "bad"}),
                         ("source_hash_scope", "only loaded files"), ("physical_input_files", {}),
                         ("physical_input_files", altered_inputs)]
            mutations += [(name, "bad") for name in ("initializer_sources_sha256", "hook_sha256",
                          "entrypoint_sha256", "preparation_script_sha256", "import_check_sha256")]
            for name, value in mutations:
                path.write_text(json.dumps({**record, name: value}))
                with self.subTest(lane=lane, field=name), mock.patch.object(validation, "run_nrg") as nrg:
                    with mock.patch.object(validation, "run_tool") as tool:
                        with self.assertRaises(ValueError):
                            self.run_case(lane=lane)
                        nrg.assert_not_called()
                        tool.assert_not_called()
                path.write_text(original)

    def test_valid_historical_source_hashes_need_not_match_current_sources(self):
        path = self.fixture / self.backend / "provenance.json"
        record = json.loads(path.read_text())
        record["initializer_source_files"]["nrginit/wilson.m"] = "0" * 64
        record["initializer_sources_sha256"] = validation.case_digest(record["initializer_source_files"])
        for name in ("hook_sha256", "entrypoint_sha256", "preparation_script_sha256", "import_check_sha256"):
            record[name] = "1" * 64
        path.write_text(json.dumps(record))
        validation.verify_fixture(self.fixture, self.backend, "runtime")


@unittest.skipUnless("legacy" in ENABLED_BACKENDS, "legacy backend disabled")
class LegacySuiteFailureTests(SuiteFailureContracts, unittest.TestCase):
    backend = "legacy"


@unittest.skipUnless("rkpw" in ENABLED_BACKENDS, "rkpw backend disabled")
class RkpwSuiteFailureTests(SuiteFailureContracts, unittest.TestCase):
    backend = "rkpw"


@unittest.skipIf(os.environ.get("NRG_TEST_DISCOVERY_CHILD") == "1", "isolated discovery child")
class BackendDiscoveryTests(unittest.TestCase):
    def test_discovery_without_disabled_backend_assets(self):
        source = Path(__file__).resolve().parent
        for enabled in (*ENABLED_BACKENDS, ""):
            with self.subTest(enabled=enabled), tempfile.TemporaryDirectory() as temporary:
                isolated = Path(temporary) / "scientific"

                def ignored(directory, names):
                    if Path(directory) == FIXTURE:
                        return {"legacy", "rkpw"} - {enabled}
                    return {"__pycache__"}

                shutil.copytree(source, isolated, ignore=ignored)
                for disabled in {"legacy", "rkpw"} - {enabled}:
                    self.assertFalse((isolated / "fixtures/asymmetric_star" / disabled).exists())
                result = subprocess.run(
                    [sys.executable, "-B", "-m", "unittest", "discover", "-s", str(isolated), "-p", "test_*.py"],
                    cwd=isolated, env={**os.environ, "NRG_TEST_CHAIN_BACKENDS": enabled, "NRG_TEST_DISCOVERY_CHILD": "1"},
                    stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
                self.assertEqual(result.returncode, 0, result.stdout)


if __name__ == "__main__":
    unittest.main()
