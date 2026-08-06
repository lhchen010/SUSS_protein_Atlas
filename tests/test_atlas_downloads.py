import base64
import io
import json
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "builders"))

import html_builder


def test_matrix_pair_stats_separate_all_pairs_from_detected_pairs():
    matrix = pd.DataFrame(
        {
            "member": ["A", "B", "C"],
            "A": [1.0, 0.6, 0.0],
            "B": [0.6, 1.0, 0.3],
            "C": [0.0, 0.3, 1.0],
        }
    )

    stats = html_builder._matrix_pair_stats(matrix)

    assert stats == {
        "mean_all": 0.3,
        "mean_detected": 0.45,
        "maximum": 0.6,
        "n_pairs": 3,
        "n_detected": 2,
        "n_possible": 3,
        "n_missing": 0,
    }


def test_matrix_stats_and_hub_do_not_zero_fill_missing_pairs():
    matrix = pd.DataFrame(
        {
            "member": ["A", "B", "C"],
            "A": [1.0, 0.8, np.nan],
            "B": [0.8, 1.0, 0.7],
            "C": [np.nan, 0.7, 1.0],
        }
    )

    stats = html_builder._matrix_pair_stats(matrix)
    hub, mean_tm, measured, expected = html_builder._hub_from_tm(
        matrix, ["A", "B", "C"]
    )

    assert stats["mean_all"] == 0.75
    assert stats["n_pairs"] == 2
    assert stats["n_missing"] == 1
    assert (hub, mean_tm, measured, expected) == ("B", 0.75, 2, 2)


def test_matrix_agreement_uses_only_mutually_measured_pairs():
    foldseek = pd.DataFrame(
        {
            "member": ["A", "B", "C", "D"],
            "A": [1.0, 0.8, np.nan, 0.5],
            "B": [0.8, 1.0, 0.7, np.nan],
            "C": [np.nan, 0.7, 1.0, 0.6],
            "D": [0.5, np.nan, 0.6, 1.0],
        }
    )
    usalign = pd.DataFrame(
        {
            "member": ["A", "B", "C", "D"],
            "A": [1.0, 0.82, 0.2, 0.52],
            "B": [0.82, 1.0, 0.72, 0.3],
            "C": [0.2, 0.72, 1.0, 0.62],
            "D": [0.52, 0.3, 0.62, 1.0],
        }
    )

    stats = html_builder._matrix_agreement_stats(foldseek, usalign)

    assert stats["n_compared"] == 4
    assert stats["pearson_r"] == 1.0
    assert stats["max_abs_diff"] == 0.02
    assert stats["n_disagree"] == 0
    assert stats["usalign_mean"] == pytest.approx(0.53)


def test_site_correlations_join_on_residue_number():
    signature = pd.DataFrame(
        {
            "resi": [30, 10, 20],
            "conservation": [3.0, 1.0, 2.0],
            "rel_sasa": [0.1, 0.3, 0.2],
        }
    )

    stats = html_builder._site_correlations(
        signature, {10: -1.0, 20: -2.0, 30: -3.0}
    )

    assert stats["conservation"] == pytest.approx(-1.0)
    assert stats["sasa"] == pytest.approx(1.0)
    assert stats["n_conservation"] == 3


def _pdb(coords):
    lines = []
    for serial, (x, y, z) in enumerate(coords, 1):
        lines.append(
            f"ATOM  {serial:5d}  CA  ALA A{serial:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 90.00           C"
        )
    return "\n".join(lines) + "\nEND\n"


def test_superpose_pdb_recovers_known_rigid_transform():
    reference = np.asarray(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
    )
    rotation = np.asarray([[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]])
    mobile = reference @ rotation + np.asarray([12.0, -7.0, 3.0])

    aligned_pdb, stats = html_builder._superpose_pdb(_pdb(mobile), _pdb(reference))

    assert stats["method"] == "ca_order"
    assert stats["n_ca"] == 4
    assert stats["rmsd"] < 1e-6
    np.testing.assert_allclose(html_builder._ca_coordinates(aligned_pdb), reference, atol=1e-3)


def test_superpose_uses_foldmason_gap_correspondence():
    reference = np.asarray(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]]
    )
    mobile_core = reference[[0, 2, 3]] + np.asarray([5.0, 2.0, -1.0])
    mobile = np.vstack([mobile_core[0], [99.0, 99.0, 99.0], mobile_core[1:]])

    _, stats = html_builder._superpose_pdb(
        _pdb(mobile), _pdb(reference), mobile_aln="AA-AA", ref_aln="A-AAA"
    )

    assert stats["method"] == "foldmason"
    assert stats["n_ca"] == 3
    assert stats["rmsd"] < 1e-6


def test_foldmason_alignment_records_preserve_gaps_and_map_members(tmp_path):
    alignment = tmp_path / "F0.aln"
    alignment.write_text(
        ">cor_A1.pdb\nACD-EF\n"
        ">cor_A2.pdb\nA-DGEF\n"
    )

    records = html_builder._read_fasta_records(alignment)
    mapped = html_builder._records_by_member(records, ["A1", "A2"])

    assert mapped == {"A1": "ACD-EF", "A2": "A-DGEF"}


def test_structure_bundle_contains_individual_pdbs_and_manifest():
    encoded = html_builder._structures_zip_b64(
        "F0", {"A1": _pdb([[0, 0, 0], [1, 0, 0], [0, 1, 0]]),
               "A2": _pdb([[2, 2, 2], [3, 2, 2], [2, 3, 2]])}
    )

    with zipfile.ZipFile(io.BytesIO(base64.b64decode(encoded))) as archive:
        assert set(archive.namelist()) == {
            "F0_structures/A1.pdb", "F0_structures/A2.pdb", "F0_structures/manifest.tsv"
        }
        assert archive.read("F0_structures/A1.pdb").startswith(b"ATOM")


def test_structure_path_accepts_a_unique_non_strain_prefix(tmp_path):
    structure = tmp_path / "legacy_A1.pdb"
    structure.write_text(_pdb([[0, 0, 0], [1, 0, 0], [0, 1, 0]]))

    resolved = html_builder._structure_path(
        "A1", input_dir=str(tmp_path), strain="current"
    )

    assert resolved == str(structure)


def test_structure_path_rejects_ambiguous_prefixes(tmp_path):
    (tmp_path / "strain1_A1.pdb").write_text(_pdb([[0, 0, 0]]))
    (tmp_path / "strain2_A1.pdb").write_text(_pdb([[1, 1, 1]]))

    with pytest.raises(RuntimeError, match="Ambiguous structures for A1"):
        html_builder._structure_path("A1", input_dir=str(tmp_path))


def test_family_workbook_contains_complete_evidence_sheets():
    matrix = pd.DataFrame({"member": ["A1", "A2"], "A1": [1.0, 0.6], "A2": [0.6, 1.0]})
    encoded = html_builder._xlsx_b64(
        fam="F0",
        members=["A1", "A2"],
        annotation=pd.DataFrame([
            {"acc": "A1", "family": "F0", "annotation_status": "complete",
             "interpro_status": "complete", "foldseek_pdb_status": "complete",
             "foldseek_afdb_status": "complete", "effectorp_status": "complete",
             "deeptmhmm_status": "complete",
             "pfam_domains": "PF00001", "pdb_hit": "1ABC", "afdbsp_name": "Protein alpha",
             "effectorp": "effector", "n_TMR": 0, "novel": False},
            {"acc": "A2", "family": "F0", "annotation_status": "partial",
             "interpro_status": "complete", "foldseek_pdb_status": "complete",
             "foldseek_afdb_status": "complete", "effectorp_status": "complete",
             "deeptmhmm_status": "failed",
             "pfam_domains": "", "pdb_hit": "", "afdbsp_name": "",
             "effectorp": "non-effector", "n_TMR": 1, "novel": False},
        ]),
        tm=matrix,
        usm=matrix,
        idm=matrix,
        blast_pairs=pd.DataFrame([{"q": "A1", "t": "A2", "pident": 18.0, "class": "core_SUSS"}]),
        sig=pd.DataFrame([{"resi": 1, "conservation": 0.8}]),
        exp=pd.DataFrame([{"acc": "A1", "control": 1.0, "infection": 4.0}]),
        pocket_entry={
            "ref": "A1",
            "fpocket_status": "complete",
            "p2rank_status": "complete",
            "fpocket": {"top_score": 2.1, "n_pockets": 1, "lining_residues": [1, 2],
                        "pockets": [{"pocket_id": 1, "score": 2.1, "lining_residues": [1, 2]}]},
            "p2rank_profile": "alphafold",
            "p2rank": {"top_score": 0.8, "top_probability": 0.73, "n_pockets": 1,
                       "lining_residues": [2, 3],
                       "pockets": [{"pocket_id": 1, "score": 0.8, "probability": 0.73,
                                    "lining_residues": [2, 3]}]},
        },
        pocket_raw={
            "fpocket_pockets": pd.DataFrame([{"pocket_id": 1, "score": 2.1, "volume": 42.0}]),
            "p2rank_pockets": pd.DataFrame([{"rank": 1, "score": 0.8, "residue_ids": "A_2 A_3"}]),
        },
        trees={"foldtree": "(A1,A2);", "lddt": "(A2,A1);"},
        tree_status={
            "metrics": {
                "foldtree": {
                    "status": "complete_with_fallback",
                    "rooting_method": "midpoint",
                    "source_stage": "pre_root",
                    "reason": "mad_unavailable_or_invalid",
                }
            }
        },
        fit_stats={"A1": {"reference": "A1", "method": "reference", "n_ca": 3, "rmsd": 0.0}},
    )

    workbook = pd.ExcelFile(io.BytesIO(base64.b64decode(encoded)))
    assert {
        "README", "members", "annotation", "foldseek_TM", "usalign_TM", "blast_identity", "blast_pairs",
        "pocket_summary", "pocket_predictions", "pocket_residues", "fpocket_pockets",
        "p2rank_pockets", "foldtree", "RNAseq", "per_site",
        "superposition",
    }.issubset(workbook.sheet_names)
    pockets = workbook.parse("pocket_predictions")
    assert set(pockets["method"]) == {"fpocket", "p2rank"}
    assert pockets.loc[pockets.method == "p2rank", "probability"].iloc[0] == 0.73
    pocket_summary = workbook.parse("pocket_summary").set_index("method")
    assert pocket_summary.loc["p2rank", "profile"] == "alphafold"
    annotation = workbook.parse("annotation")
    assert list(annotation["acc"]) == ["A1", "A2"]
    assert {"annotation_status", "pfam_domains", "pdb_hit", "afdbsp_name",
            "effectorp", "n_TMR", "novel", "interpro_status", "foldseek_pdb_status",
            "foldseek_afdb_status", "effectorp_status", "deeptmhmm_status"}.issubset(
                annotation.columns)
    foldtree = workbook.parse("foldtree").set_index("metric")
    assert foldtree.loc["foldtree", "rooting_method"] == "midpoint"
    assert foldtree.loc["foldtree", "status"] == "complete_with_fallback"


def test_singleton_workbook_keeps_direct_evidence_and_omits_family_analyses():
    annotation = pd.DataFrame([{
        "acc": "PROT1",
        "family": "singleton",
        "annotation_status": "complete",
        "interpro_status": "complete",
        "foldseek_pdb_status": "complete",
        "foldseek_afdb_status": "complete",
        "effectorp_status": "complete",
        "deeptmhmm_status": "complete",
        "pfam_domains": "PF12345",
        "interpro_entries": "IPR012345",
        "pdb_hit": "4XYZ",
        "pdb_tm": 0.63,
        "afdbsp_hit": "AF-Q9TEST-F1",
        "afdbsp_name": "Secreted test protein",
        "afdbsp_tm": 0.71,
        "effectorp": "effector",
        "n_TMR": 0,
        "novel": False,
    }])
    encoded = html_builder._xlsx_b64(
        fam="PROT1",
        members=["PROT1"],
        annotation=annotation,
        tm=None,
        usm=None,
        idm=None,
        blast_pairs=None,
        sig=None,
        exp=pd.DataFrame([{"acc": "PROT1", "control": 1.0, "infection": 8.0}]),
        pocket_entry={
            "ref": "PROT1",
            "fpocket_status": "complete",
            "fpocket": {
                "top_score": 2.4,
                "n_pockets": 1,
                "lining_residues": [2, 3],
                "pockets": [{"pocket_id": 1, "score": 2.4, "lining_residues": [2, 3]}],
            },
        },
        pocket_raw={},
        trees={},
        tree_status={},
        fit_stats={},
        analysis_kind="singleton",
    )

    workbook = pd.ExcelFile(io.BytesIO(base64.b64decode(encoded)))
    assert {"README", "members", "annotation", "pocket_summary",
            "pocket_predictions", "pocket_residues", "RNAseq"}.issubset(
                workbook.sheet_names)
    assert {"foldseek_TM", "usalign_TM", "blast_identity", "blast_pairs",
            "foldtree", "superposition", "per_site"}.isdisjoint(workbook.sheet_names)
    row = workbook.parse("annotation").iloc[0]
    assert row["pdb_hit"] == "4XYZ"
    assert row["pdb_tm"] == 0.63
    assert row["afdbsp_name"] == "Secreted test protein"
    assert row["afdbsp_tm"] == 0.71


def test_singleton_annotation_payload_preserves_foldseek_scores_and_statuses():
    payload = html_builder._annotation_payload(pd.DataFrame([{
        "acc": "PROT1",
        "family": "singleton",
        "annotation_status": "complete",
        "foldseek_pdb_status": "complete",
        "foldseek_afdb_status": "complete",
        "pdb_hit": "4XYZ",
        "pdb_tm": 0.63,
        "afdbsp_hit": "AF-Q9TEST-F1",
        "afdbsp_name": "Secreted test protein",
        "afdbsp_tm": 0.71,
        "pfam_domains": "",
        "interpro_entries": "",
        "effectorp": "non-effector",
        "n_TMR": 1,
        "novel": pd.NA,
    }]))

    member = payload["members"][0]
    assert payload["label"] == "Secreted test protein"
    assert member["pdb_tm"] == 0.63
    assert member["afdb_tm"] == 0.71
    assert member["novel"] is None
    assert member["annotation_status"] == "complete"
    assert member["foldseek_pdb_status"] == "complete"
    assert member["foldseek_afdb_status"] == "complete"


def test_domain_search_gene_aliases_include_accessions_and_segments():
    renderer = (
        ROOT / "workflow" / "builders" / "template" / "renderer.js"
    ).read_text(encoding="utf-8")

    assert (
        'field==="gene"||field==="acc"||field==="accession"||field==="protein"'
        in renderer
    )
    assert "[m.acc,m.segment_id,m.gene]" in renderer
    assert "segments overlap this label" in renderer


def test_singleton_p2rank_display_uses_probability_not_ranking_score():
    renderer = (
        ROOT / "workflow" / "builders" / "template" / "renderer.js"
    ).read_text(encoding="utf-8")
    builder = (
        ROOT / "workflow" / "builders" / "html_builder.py"
    ).read_text(encoding="utf-8")

    assert 's.pocket_metric==="probability"?"prob ":"score "' in renderer
    assert '"pocket_metric": (' in builder
    assert '"pocket_value": (' in builder


def test_build_atlas_embeds_singleton_as_independent_payload(tmp_path):
    results = tmp_path / "results"
    structures = tmp_path / "pdb"
    results.mkdir()
    structures.mkdir()
    pd.DataFrame(columns=[
        "family", "n_members", "mean_TM", "mean_identity", "suss_pct",
        "mean_pLDDT", "mean_len", "max_identity",
    ]).to_csv(results / "master.csv", index=False)
    pd.DataFrame([{
        "acc": "PROT1", "family": "singleton", "community": -1,
        "deg": 0, "plddt": 91.2, "length": 3,
    }]).to_csv(results / "members.csv", index=False)
    pd.DataFrame([{
        "acc": "PROT1", "family": "singleton", "annotation_status": "complete",
        "foldseek_pdb_status": "complete", "foldseek_afdb_status": "complete",
        "pdb_hit": "4XYZ", "pdb_tm": 0.63, "afdbsp_hit": "AF-Q9TEST-F1",
        "afdbsp_name": "Secreted test protein", "afdbsp_tm": 0.71,
        "pfam_domains": "PF12345(1-3)", "interpro_entries": "IPR012345(1-2)",
        "effectorp": "effector", "n_TMR": 0, "novel": False,
    }]).to_csv(results / "annotation.csv", index=False)
    pd.DataFrame([{
        "acc": "PROT1", "control": 1.0, "infection": 8.0,
    }]).to_csv(results / "rnaseq_expression.csv", index=False)
    (results / "seqs.fasta").write_text(">PROT1\nAAA\n")
    (structures / "cor_PROT1.pdb").write_text(
        _pdb([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
    )
    (results / "pockets.json").write_text(
        '{"PROT1":{"ref":"PROT1","p2rank":{"top_score":2.4,'
        '"top_probability":0.72,"n_pockets":1,"lining_residues":[2]}}}'
    )
    pd.DataFrame([{
        "domain_family": "D0", "n_segments": 1, "n_proteins": 1,
        "n_edges": 0, "mean_probability": 0.9, "mean_lddt": 0.7,
        "mean_aligned_residues": 3,
    }]).to_csv(results / "domain_families.csv", index=False)
    pd.DataFrame([{
        "domain_family": "D0", "segment_id": "PROT1:1-3", "acc": "PROT1",
        "start": 1, "end": 3, "length": 3, "community": 0,
    }]).to_csv(results / "domain_members.csv", index=False)
    pd.DataFrame(columns=[
        "domain_family", "source", "target", "evalue", "prob", "bits",
        "lddt", "fident", "alnlen", "shorter_coverage",
    ]).to_csv(results / "domain_edges.csv", index=False)
    pd.DataFrame(columns=[
        "source_family", "target_family", "n_edges", "mean_probability",
        "max_probability", "mean_lddt", "max_lddt", "mean_aligned_residues",
    ]).to_csv(results / "domain_cross_edges.csv", index=False)
    (results / "domain_workbench.json").write_text(json.dumps({
        "schema_version": 2,
        "families": {
            "D0": {
                "hub": "PROT1:1-3",
                "members": ["PROT1:1-3"],
                "transforms": {},
                "fit_stats": {},
                "structural_msa": {},
                "three_di_msa": {},
                "foldmason_guide_newick": "(PROT1__1-3:0.1);",
                "tree_label_map": {"PROT1__1-3": "PROT1:1-3"},
                "sequence_msa": {},
                "sequence_newick": "",
                "sequence_subgroups": [],
                "usalign_labels": ["PROT1:1-3"],
                "usalign_matrix": [[1.0]],
                "foldtree_trees": {},
                "foldtree_status": {},
                "status": {},
            }
        },
    }))
    out_html = results / "atlas.html"

    built = html_builder.build_atlas(
        master_csv=str(results / "master.csv"),
        cards_dir=str(results / "cards"),
        composition_xlsx=str(results / "composition.xlsx"),
        annotation_csv=str(results / "annotation.csv"),
        results_dir=str(results),
        out_html=str(out_html),
        mode="single",
        atlas_name="test",
        config={
            "strain": {"code": "cor", "species": "Test species"},
            "input": {"pdb_dir": str(structures)},
            "output": {"project_title": "Test singleton atlas"},
        },
    )

    html = out_html.read_text()
    assert built["families"] == 0
    assert built["singletons"] == 1
    assert '"kind": "singleton"' in html
    assert '"SINGLETONS": [{"id": "PROT1"' in html
    assert "Secreted test protein" in html
    assert '"pdb_tm": 0.63' in html
    assert '"afdb_tm": 0.71' in html
    assert '"NET": {"nodes": [], "edges": []}' in html
    assert '"REFPDB": {}' in html
    assert '"msa": {}' in html
    assert '"structures_zip_b64":' not in html
    assert '"DNET": {"nodes": [{"domain_family": "D0"' in html
    assert '"segment_id": "PROT1:1-3"' in html
    assert '"label": "PF12345"' in html
    assert '"pocket_residues": [2]' in html
    assert '"pocket_metric": "probability"' in html
    assert '"pocket_value": 0.72' in html
    assert '"package_zip_b64": "' in html
    assert '"parent_structures_zip_b64": "' in html
    assert "FoldMason structural guide tree" in html
    assert html.count("ATOM      1") == 1


def test_old_pocket_results_are_enriched_from_raw_outputs():
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        p2dir = root / "p2rank" / "F0" / "out"
        fpdir = root / "fpocket" / "F0" / "A1_out" / "pockets"
        p2dir.mkdir(parents=True)
        fpdir.mkdir(parents=True)
        pd.DataFrame([
            {"rank": 1, "score": 0.9, "residue_ids": "A_2 A_3"},
            {"rank": 2, "score": 0.5, "residue_ids": "A_8"},
        ]).to_csv(p2dir / "cor_A1.pdb_predictions.csv", index=False)
        pd.DataFrame([
            {"rank": 1, "score": 99.0, "residue_ids": "A_99"},
        ]).to_csv(p2dir / "cor_OLD.pdb_predictions.csv", index=False)
        (fpdir.parent / "A1_info.txt").write_text(
            "Pocket 1 :\n Score : 3.2\nPocket 2 :\n Score : 1.1\n"
        )
        (fpdir / "pocket1_atm.pdb").write_text(_pdb([[0, 0, 0], [1, 0, 0]]))

        enriched = html_builder._enrich_pocket_entry(
            tmp, "F0", {"ref": "A1", "fpocket": {"top_score": 3.2},
                        "p2rank": {"top_score": 0.9}}
        )

        assert len(enriched["p2rank"]["pockets"]) == 2
        assert enriched["p2rank"]["top_score"] == 0.9
        assert enriched["p2rank"]["lining_residues"] == [2, 3]
        assert len(enriched["fpocket"]["pockets"]) == 2
        assert enriched["fpocket"]["pockets"][0]["lining_residues"] == [1, 2]
        raw = html_builder._pocket_raw_tables(tmp, "F0", enriched)
        assert list(raw["p2rank_pockets"].columns) == ["rank", "score", "residue_ids"]
        assert raw["p2rank_pockets"].iloc[0]["score"] == 0.9
        assert set(raw["fpocket_pockets"]["pocket_id"]) == {1, 2}


def test_full_length_family_rnaseq_uses_run_level_expression_table():
    expression = pd.DataFrame([
        {"acc": "A1", "control": 1.0, "infection": 4.0},
        {"acc": "A2", "control": 2.0, "infection": 8.0},
        {"acc": "OTHER", "control": 9.0, "infection": 9.0},
    ])

    subset = html_builder._family_expression(expression, ["A2", "A1"])

    assert set(subset["acc"]) == {"A1", "A2"}
    assert "OTHER" not in set(subset["acc"])


def test_renderer_uses_aligned_payload_and_zip_download():
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()
    assert "PAY[fam].transforms" in renderer
    assert "alignedPdb(curFam,m)" in renderer
    assert "structures_zip_b64" in renderer
    assert "All structures (ZIP)" in renderer
    assert "All structures (multi-PDB)" not in renderer


def test_backend_atlas_uses_lazy_artifacts_and_streaming_routes():
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()
    builder = (ROOT / "workflow" / "builders" / "html_builder.py").read_text()
    portal = (ROOT / "portal" / "suss_portal.py").read_text()

    assert "function fetchStructure" in renderer
    assert "function fetchReference" in renderer
    assert "function hasReference" in renderer
    assert 'artifactUrl("structure",acc)' in renderer
    assert 'artifactUrl("reference",key)' in renderer
    assert 'artifactUrl("xlsx",fam)' in renderer
    assert 'artifactUrl("structures",fam)' in renderer
    assert "hasStruct=BACKEND.enabled||" in renderer
    assert "if(!pdb&&BACKEND.enabled&&!structureFetchFailed[m.acc])" in renderer
    assert 'if mode == "backend"' in builder
    assert "store_download" in builder
    assert "store_reference" in builder
    assert "REFAVAIL=REFAVAIL" in builder
    assert 'elif u.path == "/artifact"' in portal
    assert "def _stream_file" in portal
    assert 'html_mode="backend"' in portal
    assert 'kind == "reference"' in portal


def test_network_search_supports_annotation_fields_and_highlighting():
    prefix = (ROOT / "workflow" / "builders" / "template" / "prefix.html").read_text()
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()

    assert 'id="searchinput"' in prefix
    assert 'id="searchstatus"' in prefix
    assert 'id="clearsearch"' in prefix
    assert "function familyMatches" in renderer
    assert "function applyNetworkSearch" in renderer
    assert 'field==="gene"' in renderer
    assert 'field==="annotation"' in renderer
    assert 'field==="effector"' in renderer
    assert 'field==="tmr"' in renderer
    assert 'field==="structtm"' in renderer
    assert "networkNodes.update" in renderer
    assert "network.focus" in renderer


def test_singleton_workbench_is_separate_from_family_network():
    prefix = (ROOT / "workflow" / "builders" / "template" / "prefix.html").read_text()
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()

    assert 'id="modesingletons"' in prefix
    assert 'id="singletons"' in prefix
    assert "function showSingleton" in renderer
    assert "function singletonMatches" in renderer
    assert "function buildSingletonStructPane" in renderer
    assert "function singletonDlbtn" in renderer
    assert "function basePdb" in renderer
    assert "function esmPdb" in renderer
    assert "function pdbWithValues" in renderer
    assert "Foldseek PDB100" in renderer
    assert "Foldseek AFDB / Swiss-Prot" in renderer
    assert 'setMode(\\\'quality\\\')' in renderer
    assert "function singletonTab" in renderer


def test_sequence_viewer_and_alignment_downloads_are_exposed_without_singleton_msa():
    prefix = (ROOT / "workflow" / "builders" / "template" / "prefix.html").read_text()
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()
    builder = (ROOT / "workflow" / "builders" / "html_builder.py").read_text()

    assert "Sequence + MSA" in renderer
    assert "function dlSeqs" in renderer
    assert "function dlMsa" in renderer
    assert "function buildSequencePane" in renderer
    assert "function renderSequenceViewer" in renderer
    assert "Sequence MSA" in renderer
    assert "Structural MSA (AA)" in renderer
    assert "Structural MSA (3Di)" in renderer
    assert "function dlAlignment" in renderer
    assert "MAFFT_sequence_MSA" in renderer
    assert "FoldMason_AA_MSA" in renderer
    assert "singletonTab(3)" in renderer
    assert "not applicable to a singleton" in renderer.lower()
    assert ".sequence-view" in prefix
    assert "sequence_msa=sequence_msa" in builder
    assert "structural_msa=msa" in builder


def test_domain_family_mode_and_structure_search_are_exposed():
    prefix = (ROOT / "workflow" / "builders" / "template" / "prefix.html").read_text()
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()
    builder = (ROOT / "workflow" / "builders" / "html_builder.py").read_text()
    portal = (ROOT / "portal" / "suss_portal.py").read_text()
    snakefile = (ROOT / "workflow" / "Snakefile").read_text()

    assert 'id="modedomains"' in prefix
    assert 'id="domains"' in prefix
    assert "function renderDomainTable" in renderer
    assert "function dlDomainDiagnostics" in renderer
    assert "Borderline local hit" in renderer
    assert "domain_match_diagnostics.csv" in snakefile
    assert "function showDomain" in renderer
    assert renderer.count("function showDomain(id)") == 1
    assert "function showDomainSegment" in renderer
    assert "function renderDomainStructure" in renderer
    assert "function toggleDomainSelection" in renderer
    assert "function selectDomainMembers" in renderer
    assert "function dlDomainSelectedSuperposition" in renderer
    assert 'Superpose selected (' in renderer
    assert "Selected full proteins aligned" in renderer
    assert "All domain structures ZIP" in renderer
    assert "All parent structures ZIP" in renderer
    assert "Complete D-family ZIP" in renderer
    assert "Detector-native pocket records stay in the workbook/package" in builder
    assert "function buildDomainTreesPane" in renderer
    assert "function buildDomainConservationPane" in renderer
    assert "function setDomainBackground" in renderer
    assert "function setDomainRep" in renderer
    assert "function setDomainPocket" in renderer
    assert "Domain-segment BLASTp identity" in renderer
    assert "Sequence conservation" in renderer
    assert "Rate4Site uses only independently searched D-segment" in renderer
    assert "function setViewerBackground" in renderer
    assert "Foldseek TM matrix was not available" in renderer
    assert "Independent US-align matrix was not available" in renderer
    assert "FoldTree structural relationship" in renderer
    assert "Independent D-segment BLASTp subgroup" in renderer
    assert "FoldMason structural MSA" in renderer
    assert "domain-architecture" in renderer
    assert "DOMAIN_EDGES" in renderer
    assert "DNET.edges" in renderer
    assert "DOMAIN_FAMILIES" in renderer
    assert "action=/search-structure" in portal
    assert "structure_search_index.csv" in portal
    assert "protein={urllib.parse.quote(hit['acc'])}" in portal
    assert "&open={urllib.parse.quote(family)}" in portal
    assert 'segment={urllib.parse.quote(hit["acc"])}' in portal
    assert "function openAtlasTargetFromUrl" in renderer
    assert 'params.get("protein")' in renderer
    assert 'params.get("open")' in renderer
    assert 'params.get("segment")' in renderer
    assert "except subprocess.TimeoutExpired" in portal
    assert "shutil.rmtree(search_dir, ignore_errors=True)" in portal
    assert "Domain-aware Foldseek" in portal
    assert 'name=analysis_scope value=both checked' in portal
    assert 'name=analysis_scope value=domain' in portal
    assert 'analysis_scope in ("full", "both")' in portal
    assert "Sequence MSA / tree (F proteins + independent D segments)" in portal
    assert 'kind == "domain_package"' in portal
    assert 'kind == "domain_parents"' in portal


def test_tree_svg_titles_report_the_actual_evidence_source():
    sequence_svg = html_builder._newick_to_svg(
        "(A:0.1,B:0.2);", title="Sequence tree · MAFFT + FastTree"
    )
    structural_svg = html_builder._newick_to_svg(
        "(A:0.1,B:0.2);", title="FoldTree structural tree · lddt"
    )

    assert "Sequence tree · MAFFT + FastTree" in sequence_svg
    assert "FoldTree structural tree · lddt" in structural_svg


def test_domain_downloads_include_segments_parents_and_manifest():
    members = [
        {"segment_id": "A:1-2", "acc": "A", "start": 1, "end": 2},
        {"segment_id": "B:2-3", "acc": "B", "start": 2, "end": 3},
    ]
    structures = {
        "A": _pdb([[0, 0, 0], [1, 0, 0], [2, 0, 0]]),
        "B": _pdb([[0, 1, 0], [1, 1, 0], [2, 1, 0]]),
    }

    segment_zip = base64.b64decode(
        html_builder._domain_structures_zip_b64(
            "D0", members, structures, segments_only=True
        )
    )
    parent_zip = base64.b64decode(
        html_builder._domain_structures_zip_b64(
            "D0", members, structures, segments_only=False
        )
    )

    with zipfile.ZipFile(io.BytesIO(segment_zip)) as archive:
        assert "D0_domain_structures/A_1-2.pdb" in archive.namelist()
        assert "D0_domain_structures/manifest.csv" in archive.namelist()
    with zipfile.ZipFile(io.BytesIO(parent_zip)) as archive:
        assert "D0_parent_structures/A.pdb" in archive.namelist()
        assert "D0_parent_structures/B.pdb" in archive.namelist()

    package = base64.b64decode(
        html_builder._domain_package_b64(
            "D0",
            members,
            [],
            {
                "hub": "A:1-2",
                "transforms": {"A:1-2": {"rotation": [], "translation": []}},
                "fit_stats": {},
                "structural_msa": {"A:1-2": "AA"},
                "three_di_msa": {"A:1-2": "xx"},
                "sequence_msa": {},
                "sequence_subgroups": [
                    {
                        "id": "S0",
                        "msa": {"A:1-2": "AA"},
                        "newick": "(A__1-2:0.1);",
                    }
                ],
                "foldmason_guide_newick": "(A__1-2:0.1);",
                "foldtree_trees": {},
                "foldtree_status": {},
                "structural_conservation": [0.8, None],
                "sequence_identity_labels": ["A:1-2"],
                "sequence_identity_matrix": [[1.0]],
                "structural_alignment_identity_labels": ["A:1-2"],
                "structural_alignment_identity_matrix": [[1.0]],
                "domain_blast_edges": [
                    {
                        "q": "A:1-2",
                        "t": "B:2-3",
                        "pident": 25.0,
                        "min_coverage": 0.8,
                    }
                ],
                "sequence_conservation": {
                    "A:1-2": {"1": 0.8, "2": -0.2}
                },
                "usalign_labels": ["A:1-2"],
                "usalign_matrix": [[1.0]],
                "status": {"foldmason": "complete"},
            },
            structures,
            {"A": "AAA", "B": "BBB"},
            base64.b64encode(b"workbook").decode(),
        )
    )
    with zipfile.ZipFile(io.BytesIO(package)) as archive:
        names = archive.namelist()
        assert "D0_domain_family/trees/FoldMason_guide.nwk" in names
        assert (
            "D0_domain_family/alignments/sequence_subgroups/S0_MAFFT.fasta"
            in names
        )
        assert "D0_domain_family/superposition/transforms.json" in names
        assert (
            "D0_domain_family/tables/foldmason_structural_conservation.csv"
            in names
        )
        assert "D0_domain_family/tables/domain_sequence_identity.csv" in names
        assert "D0_domain_family/tables/foldmason_AA_identity.csv" in names
        assert "D0_domain_family/tables/domain_blastp_hits.csv" in names
        assert (
            "D0_domain_family/tables/rate4site_sequence_conservation.csv"
            in names
        )


def test_domain_workbook_keeps_explicit_sequence_msa_status_when_not_applicable():
    encoded = html_builder._domain_xlsx_b64(
        "D0",
        [
            {
                "segment_id": "A:1-2",
                "acc": "A",
                "start": 1,
                "end": 2,
                "expression": {},
                "p2rank": {},
                "fpocket": {},
                "esm_values": {},
                "overlap_annotations": [],
            }
        ],
        [],
        {
            "sequence_msa": {},
            "structural_msa": {"A:1-2": "AA"},
            "three_di_msa": {"A:1-2": "xx"},
            "sequence_subgroups": [],
            "usalign_labels": ["A:1-2"],
            "usalign_matrix": [[1.0]],
            "fit_stats": {},
            "structural_conservation": [0.8, None],
            "sequence_identity_labels": ["A:1-2"],
            "sequence_identity_matrix": [[1.0]],
            "structural_alignment_identity_labels": ["A:1-2"],
            "structural_alignment_identity_matrix": [[1.0]],
            "domain_blast_edges": [
                {
                    "q": "A:1-2",
                    "t": "A:1-2",
                    "pident": 100.0,
                }
            ],
            "sequence_conservation": {"A:1-2": {"1": 0.8}},
        },
    )

    workbook = pd.ExcelFile(io.BytesIO(base64.b64decode(encoded)))
    assert "sequence_MSA" in workbook.sheet_names
    assert "sequence_identity" in workbook.sheet_names
    assert "foldmason_AA_identity" in workbook.sheet_names
    assert "domain_blastp_hits" in workbook.sheet_names
    assert "sequence_conservation" in workbook.sheet_names
    sequence_msa = workbook.parse("sequence_MSA")
    assert sequence_msa.loc[0, "status"] == "not_applicable"
    assert "No reciprocal-coverage" in sequence_msa.loc[0, "reason"]


def test_checkpoint_member_lists_do_not_own_family_analysis_directory():
    snakefile = (ROOT / "workflow" / "Snakefile").read_text()
    builder = (
        ROOT / "workflow" / "builders" / "html_builder.py"
    ).read_text()

    assert 'MEMBER_DIR = f"{RESULTS}/family_members"' in snakefile
    assert "famdir=directory(MEMBER_DIR)" in snakefile
    assert snakefile.count(
        'famfile=f"{MEMBER_DIR}/{{fam}}.members.txt"'
    ) == 7
    assert 'famfile=f"{FAM_DIR}/{{fam}}.members.txt"' not in snakefile
    signature = (
        ROOT / "workflow" / "scripts" / "signature.py"
    ).read_text()
    assert "famfile = snakemake.input.famfile" in signature
    assert (
        'os.path.dirname(os.path.dirname(r4s_path)), f"{fam}.members.txt"'
        not in signature
    )
    pocket_script = (
        ROOT / "workflow" / "scripts" / "sasa_pocket.py"
    ).read_text()
    esm_script = (
        ROOT / "workflow" / "scripts" / "esm_scan.py"
    ).read_text()
    assert "famdir = snakemake.input.famdir" in pocket_script
    assert "family_member_dir = snakemake.input.famdir" in esm_script
    assert (
        'os.path.join(os.path.dirname(out_sasa), "families")'
        not in pocket_script
    )
    assert (
        'os.path.join(os.path.dirname(out_csv), "families")'
        not in esm_script
    )
    assert '"family_members", f"{fam}.members.txt"' in builder


def test_small_disconnected_domain_overview_uses_fixed_layout():
    renderer = (
        ROOT / "workflow" / "builders" / "template" / "renderer.js"
    ).read_text()

    assert "fixedSmallNetwork=edges.length===0" in renderer
    assert "node.fixed={x:true,y:true}" in renderer
    assert "physics:fixedSmallNetwork?false:" in renderer


def test_p2rank_alphafold_profile_is_configurable_and_reference_safe():
    snakefile = (ROOT / "workflow" / "Snakefile").read_text()
    pocket_script = (ROOT / "workflow" / "scripts" / "sasa_pocket.py").read_text()
    esm_script = (ROOT / "workflow" / "scripts" / "esm_scan.py").read_text()
    portal = (ROOT / "portal" / "suss_portal.py").read_text()

    assert 'get("p2rank_profile", "alphafold")' in snakefile
    assert "threads: 4" in snakefile
    assert 'p2cmd.extend(["-c", p2rank_profile])' in pocket_script
    assert "shutil.rmtree(out_dir, ignore_errors=True)" in pocket_script
    assert "expected one P2Rank predictions CSV" in pocket_script
    assert '"n_pockets": 0' in pocket_script
    assert '"pockets": []' in pocket_script
    assert "targets = sorted(keep)" in pocket_script
    assert "ThreadPoolExecutor(max_workers=workers)" in pocket_script
    assert "profile_marker" in pocket_script
    assert "cached_profiles" in pocket_script
    assert 'scope == "all_proteins"' in esm_script
    assert 'scope == "domain_members"' in esm_script
    assert 'scope == "family_representatives"' in esm_script
    assert 'scope == "representatives"' in esm_script
    assert "set(family_refs.values()) | singletons" in esm_script
    assert 'f"*_{accession}-res-in-matrix.csv"' in esm_script
    domain_script = (
        ROOT / "workflow" / "scripts" / "domain_workbench.py"
    ).read_text()
    assert "json.dumps(json_safe(payload)" in domain_script
    assert "not math.isfinite(value)" in domain_script
    assert "domain workbench: disabled" in domain_script
    assert 'payload = {"schema_version": 4' in domain_script
    assert "snakemake.input.domain_blastp" in domain_script
    assert "snakemake.input.blastp" not in domain_script
    assert "AlphaFold / predicted structures" in portal
    assert 'cfg.setdefault("pocket", {}).update' in portal
    assert 'domain_foldtree=ck("domain_foldtree")' in portal


def test_full_length_rnaseq_renderer_has_explicit_missing_state():
    renderer = (
        ROOT / "workflow" / "builders" / "template" / "renderer.js"
    ).read_text(encoding="utf-8")

    assert "RNA-seq data were not available for this family." in renderer
    assert "innerHTML=a.rna_svg?" in renderer


def test_structural_color_direction_and_missing_coverage_are_explicit():
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()

    assert 'min:100,max:0' in renderer
    assert 'min:qmax,max:qmin' in renderer
    assert 'min:d.esm_max,max:d.esm_min' in renderer
    assert "structural_scored_resi" in renderer
    assert "insufficient pair coverage" in renderer


def test_viewer_context_conservation_mapping_and_analysis_axes_are_explicit():
    renderer = (
        ROOT / "workflow" / "builders" / "template" / "renderer.js"
    ).read_text()
    prefix = (
        ROOT / "workflow" / "builders" / "template" / "prefix.html"
    ).read_text()

    assert "function domainStructuralValues" in renderer
    assert "wb.structural_msa" in renderer
    assert "wb.structural_conservation" in renderer
    assert "domainStructuralValues(seg,wb)" in renderer
    assert 'seg.segment_id===wb.hub?' not in renderer
    assert 'color:"#b9c1c5",opacity:domainRepMode==="surface"?.52:.9' in renderer
    assert 'preferred="quality"' in renderer
    assert 'id="bquality" class="on"' in renderer
    assert "function setAnalysisAxis" in renderer
    assert "parent-context ESM-1b" in renderer
    assert "D-segment MAFFT + Rate4Site" in renderer
    assert "Domain-segment BLASTp identity" in renderer
    assert 'id="scopefull"' in prefix
    assert 'id="scopedomain"' in prefix
    assert 'id="fullModeTabs"' in prefix
    assert 'id="domainModeTabs"' in prefix
    assert "Full-length analysis" in prefix
    assert "Domain analysis" in prefix


def test_portal_and_atlas_share_the_scientific_console_visual_contract():
    portal = (ROOT / "portal" / "suss_portal.py").read_text()
    prefix = (
        ROOT / "workflow" / "builders" / "template" / "prefix.html"
    ).read_text()
    renderer = (
        ROOT / "workflow" / "builders" / "template" / "renderer.js"
    ).read_text()

    assert "--cyan:#82ddf6" in portal
    assert "STRUCTURE-FIRST COMPARATIVE PROTEOMICS" in portal
    assert "class=brand-mark" in portal
    assert "class=recent-table" in portal
    assert "@media(max-width:620px)" in portal

    assert "--s-cyan:#82ddf6" in prefix
    assert "STRUCTURE-FIRST ATLAS" in prefix
    assert ".mode-tabs.domain-axis button.on{background:var(--s-amber)" in prefix
    assert "#primary{border-right:0;background:var(--s-paper)}" in prefix
    assert ".dlbar{margin-top:8px" in prefix
    assert "background:#f7f9fb" not in renderer
    assert "@media(max-width:700px)" in prefix


def test_domain_sequence_search_is_a_declared_independent_workflow():
    snakefile = (ROOT / "workflow" / "Snakefile").read_text()

    assert "rule domain_sequences:" in snakefile
    assert "rule domain_blastp:" in snakefile
    assert 'domain_blastp=f"{RESULTS}/domain_blastp_allvsall.tsv"' in snakefile
    assert 'domain_manifest=f"{RESULTS}/domain_sequence_manifest.csv"' in snakefile


def test_domain_annotations_require_coordinate_overlap():
    row = pd.Series({
        "pfam_domains": "Thioredoxin(12-99) | Remote domain(150-220)",
        "interpro_entries": "IPR000001(20-80)",
    })

    overlaps = html_builder._overlapping_annotations(row, 2, 127)

    assert {item["label"] for item in overlaps} == {
        "Thioredoxin", "IPR000001"
    }
    assert all(item["overlap"] > 0 for item in overlaps)
