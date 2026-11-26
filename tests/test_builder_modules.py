"""Tests for modularized builder packages and importable domain_workbench helpers."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "builders"))
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import export_bundles
import html_builder
import pdb_io
import svg_charts
import domain_workbench


def test_html_builder_reexports_modular_helpers():
    assert html_builder._superpose_pdb is pdb_io._superpose_pdb
    assert html_builder._structure_path is pdb_io._structure_path
    assert html_builder._newick_to_svg is svg_charts._newick_to_svg
    assert html_builder._svg_matrix is svg_charts._svg_matrix
    assert html_builder._xlsx_b64 is export_bundles._xlsx_b64
    assert html_builder._domain_xlsx_b64 is export_bundles._domain_xlsx_b64


def test_domain_workbench_is_importable_without_snakemake_global():
    assert callable(domain_workbench.main)
    assert callable(domain_workbench.aligned_identity_matrix)
    mat = domain_workbench.aligned_identity_matrix(
        {"seg1": "ACDE", "seg2": "ACDF"},
        ["seg1", "seg2"],
    )
    assert mat[0][0] == 1.0
    assert mat[0][1] == 0.75


def test_annotate_has_no_stale_accre_reference():
    src = (ROOT / "workflow" / "scripts" / "annotate.py").read_text(encoding="utf-8")
    assert "ACCRE" not in src


def test_matrices_filters_blast_with_configured_thresholds():
    src = (ROOT / "workflow" / "scripts" / "matrices.py").read_text(encoding="utf-8")
    sf = (ROOT / "workflow" / "Snakefile").read_text(encoding="utf-8")
    assert "select_blast_relationships" in src
    assert "evalue=config[\"classification\"][\"blast_evalue\"]" in sf


def test_html_builder_counts_cysteines_without_requiring_rate4site_pdb():
    src = (ROOT / "workflow" / "builders" / "html_builder.py").read_text(encoding="utf-8")
    assert "cys_candidates" in src
    assert "_structural_conservation.pdb" in src
