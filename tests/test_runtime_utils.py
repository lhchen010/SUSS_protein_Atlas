import io
import importlib.util
import os
import sys
import tarfile
import tempfile
import unittest
from pathlib import Path

import yaml


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from runtime_utils import (
    analysis_status,
    novel_call,
    resolve_executable,
    stage_qc_structures,
    symmetric_tm,
)


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

    def test_stage_qc_structures_excludes_failed_rows(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            pdb_dir = root / "pdb"
            output_dir = root / "qc_pdb"
            pdb_dir.mkdir()
            (pdb_dir / "strain_A.pdb").write_text("PASS")
            (pdb_dir / "strain_B.pdb").write_text("FAIL")
            qc = root / "qc.csv"
            qc.write_text(
                "acc,fn,length,mean_plddt,pass\n"
                "A,strain_A.pdb,25,83.0,True\n"
                "B,strain_B.pdb,20,45.0,False\n"
            )

            staged = stage_qc_structures(qc, pdb_dir, output_dir)

            self.assertEqual(staged, ["strain_A.pdb"])
            self.assertEqual((output_dir / "strain_A.pdb").read_text(), "PASS")
            self.assertFalse((output_dir / "strain_B.pdb").exists())

    def test_stage_qc_structures_fails_clearly_when_nothing_passes(self):
        with tempfile.TemporaryDirectory() as root:
            root = Path(root)
            pdb_dir = root / "pdb"
            pdb_dir.mkdir()
            (pdb_dir / "strain_A.pdb").write_text("FAIL")
            qc = root / "qc.csv"
            qc.write_text("acc,fn,pass\nA,strain_A.pdb,False\n")

            with self.assertRaisesRegex(ValueError, "No structures passed QC"):
                stage_qc_structures(qc, pdb_dir, root / "qc_pdb")


class PortalArchiveTests(unittest.TestCase):
    def test_pdb_name_normalization_replaces_existing_prefix(self):
        self.assertEqual(
            self.portal._normalized_pdb_name("old_TDZ13209.1.pdb", "new"),
            "new_TDZ13209.1.pdb",
        )

    def test_upload_filename_does_not_truth_test_field_storage(self):
        class Upload:
            filename = "query.pdb"

            def __bool__(self):
                raise TypeError("Cannot be converted to bool.")

        self.assertEqual(self.portal._uploaded_filename(Upload()), "query.pdb")
        self.assertEqual(self.portal._uploaded_filename(None), "")

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

    def test_missing_job_config_is_initialized_from_4070_reference(self):
        with tempfile.TemporaryDirectory() as engine:
            config = Path(engine) / "config"
            config.mkdir()
            reference = config / "config.yaml.4070.example"
            reference.write_text("signals:\n  esm_scope: family_representatives\n")

            path, initialized = self.portal._ensure_engine_config(engine)

            self.assertTrue(initialized)
            self.assertEqual(
                Path(path).read_text(),
                reference.read_text(),
            )

    def test_qc_settings_are_user_configurable_and_validated(self):
        self.assertEqual(
            self.portal._parse_qc_settings(
                {"min_length": "20", "max_length": "1500", "min_plddt": "45.5"}.get
            ),
            {"min_length": 20, "max_length": 1500, "min_plddt": 45.5},
        )
        with self.assertRaisesRegex(ValueError, "cannot exceed"):
            self.portal._parse_qc_settings(
                {"min_length": "100", "max_length": "20", "min_plddt": "50"}.get
            )

    def test_write_config_persists_user_qc_settings(self):
        with tempfile.TemporaryDirectory() as engine:
            config_dir = Path(engine) / "config"
            config_dir.mkdir()
            config_path = config_dir / "config.yaml"
            config_path.write_text("{}\n")
            meta = {
                "code": "bc",
                "species": "Botrytis cinerea",
                "host_range": "plants",
                "is_colleto": False,
                "rnaseq_xlsx": "",
                "tm": 0.5,
                "coverage": 0.5,
                "res": 1.0,
                "minfam": 2,
                "domain_evalue": 0.001,
                "domain_prob": 0.5,
                "domain_min_aln": 40,
                "domain_lddt": 0.5,
                "domain_alntm": 0.0,
                "domain_coverage": 0.0,
                "evalue": 0.001,
                "blast_coverage": 0.5,
                "min_length": 20,
                "max_length": 1500,
                "min_plddt": 45.0,
                "p2rank_profile": "alphafold",
                "steps": {},
                "rnaseq_mode": "false",
                "project_title": "",
                "analysis_scope": "both",
            }

            self.portal._write_config(engine, meta)

            self.assertEqual(
                yaml.safe_load(config_path.read_text())["qc"],
                {"min_plddt": 45.0, "min_length": 20, "max_length": 1500},
            )


if __name__ == "__main__":
    unittest.main()
