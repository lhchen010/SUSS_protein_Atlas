"""Tests for centralized flexible protein ID parsing and preflight validation."""

import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import id_utils
import validate_inputs


def test_extract_accession_supports_genbank_refseq_uniprot_afdb_and_custom():
    # GenBank (legacy 100% compatibility)
    assert id_utils.extract_accession("cor_TDZ13209.1.pdb") == "TDZ13209.1"
    assert id_utils.extract_accession("/tmp/cor_TDZ13209.1.pdb") == "TDZ13209.1"
    # RefSeq
    assert id_utils.extract_accession("bc_XP_001552341.1.pdb") == "XP_001552341.1"
    assert id_utils.extract_accession(">NP_001234.2 secreted protein") == "NP_001234.2"
    # AlphaFold DB
    assert id_utils.extract_accession("AF-P12345-F1-model_v4.cif") == "P12345"
    assert id_utils.extract_accession("AF-A0A0K2H345-F1-model_v4.pdb.gz") == "A0A0K2H345"
    # UniProt piped & bare
    assert id_utils.extract_accession(">tr|A0A0K2H345|A0A0K2H345_COLGL effector") == "A0A0K2H345"
    assert id_utils.extract_accession("strain_Q9Y261.cif") == "Q9Y261"
    # Custom fallback stem
    assert id_utils.extract_accession("myco_gene0042.pdb", strain="myco") == "gene0042"
    assert id_utils.extract_accession("myco_gene0042.pdb", strip_strain_prefix=True) == "gene0042"
    # Custom user regex
    assert (
        id_utils.extract_accession(
            "jgi_Colsp1_99812_model.cif",
            custom_pattern=r"Colsp1_(\d+)",
        )
        == "99812"
    )


def test_normalize_structure_filename_preserves_extension_and_replaces_prefix():
    assert (
        id_utils.normalize_structure_filename("old_TDZ13209.1.pdb", "new")
        == "new_TDZ13209.1.pdb"
    )
    assert (
        id_utils.normalize_structure_filename("AF-P12345-F1-model_v4.cif", "bc")
        == "bc_P12345.cif"
    )
    assert (
        id_utils.normalize_structure_filename("bc_XP_001552341.1.mmcif.gz", "bc2")
        == "bc2_XP_001552341.1.mmcif.gz"
    )


def test_validate_inputs_accepts_non_genbank_ids_and_detects_duplicates(tmp_path):
    pdb_dir = tmp_path / "pdb"
    pdb_dir.mkdir()
    sample_pdb = (
        "ATOM      1  CA  ALA A   1       0.000   0.000   0.000  1.00 88.50           C\n"
        "ATOM      2  CA  GLY A   2       3.800   0.000   0.000  1.00 91.20           C\n"
        "END\n"
    )
    (pdb_dir / "AF-P12345-F1-model_v4.pdb").write_text(sample_pdb, encoding="utf-8")
    (pdb_dir / "bc_XP_001552341.1.pdb").write_text(sample_pdb, encoding="utf-8")

    cfg = {
        "strain": {"code": "bc"},
        "input": {"pdb_dir": str(pdb_dir)},
        "qc": {"min_length": 2, "max_length": 1000, "min_plddt": 50},
    }
    errors, warns, info = validate_inputs.validate(cfg)
    assert errors == []
    assert info["n_structures"] == 2

    # Duplicate accession across two filenames should fail preflight
    (pdb_dir / "bc_P12345.pdb").write_text(sample_pdb, encoding="utf-8")
    errors2, _, _ = validate_inputs.validate(cfg)
    assert any("duplicate protein accession" in e for e in errors2)
