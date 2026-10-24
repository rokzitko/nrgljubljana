import tempfile
import unittest
from pathlib import Path

from contract import compare_files


class ComparisonDiagnostics(unittest.TestCase):
    def test_equal_files_leave_no_artifacts(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            native, reference = root / "native", root / "reference"
            native.write_bytes(b"1.25\n")
            reference.write_bytes(native.read_bytes())
            compare_files(native, reference, root)
            self.assertEqual(sorted(path.name for path in root.iterdir()), ["native", "reference"])

    def test_mismatch_survives_fixture_cleanup_with_bounded_preview(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            with tempfile.TemporaryDirectory(dir=root) as fixture:
                native, reference = Path(fixture) / "native", Path(fixture) / "reference"
                native_bytes = b"1.2500000000000001\n" * 100
                reference_bytes = b"1.25\n" * 100
                native.write_bytes(native_bytes)
                reference.write_bytes(reference_bytes)
                with self.assertRaises(AssertionError) as failure:
                    compare_files(native, reference, root)
            artifacts, = root.glob("workflow-mismatch-*")
            self.assertEqual((artifacts / "native").read_bytes(), native_bytes)
            self.assertEqual((artifacts / "reference").read_bytes(), reference_bytes)
            diff = (artifacts / "diff.log").read_text()
            self.assertIn("-1.25\n", diff)
            self.assertIn("+1.2500000000000001\n", diff)
            self.assertGreater(len(diff.splitlines()), 80)
            self.assertIn(str(artifacts), str(failure.exception))
            self.assertLessEqual(len(str(failure.exception).splitlines()), 82)


if __name__ == "__main__":
    unittest.main()
