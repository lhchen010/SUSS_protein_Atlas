import io
import importlib.util
import os
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from runtime_utils import analysis_status, novel_call, resolve_executable, symmetric_tm


class RuntimeUtilsTests(unittest.TestCase):
    def test_resolve_executable_accepts_absolute_and_path_names(self):
        self.assertEqual(resolve_executable(sys.executable, "python"), str(Path(sys.executable).resolve()))
        self.assertTrue(Path(resolve_executable("sh", "shell")).is_file())

    def test_resolve_executable_rejects_missing_required_tool(self):
        with self.assertRaises(FileNotFoundError):
            resolve_executable("definitely-not-a-suss-tool", "missing")

    def test_symmetric_tm_modes(self):
        self.assertEqual(symmetric_tm(0.4, 0.8, "min"), 0.4)
        self.assertEqual(symmetric_tm(0.4, 0.8, "max"), 0.8)
        self.assertAlmostEqual(symmetric_tm(0.4, 0.8, "mean"), 0.6)
        with self.assertRaises(ValueError):
            symmetric_tm(0.4, 0.8, "median")

    def test_annotation_status_and_novel_are_tri_state(self):
        statuses = {"domain": "complete", "pdb": "complete", "afdb": "complete"}
        self.assertEqual(analysis_status(False, statuses, tuple(statuses)), "not_run")
        self.assertEqual(analysis_status(True, statuses, tuple(statuses)), "complete")
        statuses["afdb"] = "not_run"
        self.assertEqual(analysis_status(True, statuses, tuple(statuses)), "partial")
        self.assertIsNone(novel_call(False, False, False))
        self.assertTrue(novel_call(False, False, True))
        self.assertFalse(novel_call(True, False, True))


class PortalArchiveTests(unittest.TestCase):
    def test_pdb_name_normalization_replaces_existing_prefix(self):
        self.assertEqual(
            self.portal._normalized_pdb_name("old_TDZ13209.1.pdb", "new"),
            "new_TDZ13209.1.pdb",
        )

    @classmethod
    def setUpClass(cls):
        cls.runs = tempfile.TemporaryDirectory()
        os.environ["SUSS_RUNS_DIR"] = cls.runs.name
        spec = importlib.util.spec_from_file_location("suss_portal_test", ROOT / "portal" / "suss_portal.py")
        cls.portal = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.portal)

    @classmethod
    def tearDownClass(cls):
        cls.runs.cleanup()

    def test_engine_extraction_rejects_path_traversal(self):
        payload = io.BytesIO()
        with tarfile.open(fileobj=payload, mode="w:gz") as archive:
            info = tarfile.TarInfo("../escape.txt")
            data = b"unsafe"
            info.size = len(data)
            archive.addfile(info, io.BytesIO(data))
        payload.seek(0)
        with tempfile.TemporaryDirectory() as destination, tarfile.open(fileobj=payload, mode="r:gz") as archive:
            with self.assertRaises(ValueError):
                self.portal._safe_extract_engine(archive, destination)


if __name__ == "__main__":
    unittest.main()
