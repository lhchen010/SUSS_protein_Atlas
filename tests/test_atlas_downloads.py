import base64
import io
import sys
import tempfile
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "builders"))

import html_builder


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
            "p2rank": {"top_score": 0.8, "n_pockets": 1, "lining_residues": [2, 3],
                       "pockets": [{"pocket_id": 1, "score": 0.8, "lining_residues": [2, 3]}]},
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
        "pfam_domains": "PF12345", "interpro_entries": "IPR012345",
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
        '{"PROT1":{"ref":"PROT1","fpocket":{"top_score":2.4,'
        '"n_pockets":1,"lining_residues":[2]}}}'
    )
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
        ]).to_csv(p2dir / "test_predictions.csv", index=False)
        (fpdir.parent / "A1_info.txt").write_text(
            "Pocket 1 :\n Score : 3.2\nPocket 2 :\n Score : 1.1\n"
        )
        (fpdir / "pocket1_atm.pdb").write_text(_pdb([[0, 0, 0], [1, 0, 0]]))

        enriched = html_builder._enrich_pocket_entry(
            tmp, "F0", {"ref": "A1", "fpocket": {"top_score": 3.2},
                        "p2rank": {"top_score": 0.9}}
        )

        assert len(enriched["p2rank"]["pockets"]) == 2
        assert len(enriched["fpocket"]["pockets"]) == 2
        assert enriched["fpocket"]["pockets"][0]["lining_residues"] == [1, 2]
        raw = html_builder._pocket_raw_tables(tmp, "F0", enriched)
        assert list(raw["p2rank_pockets"].columns) == ["rank", "score", "residue_ids"]
        assert set(raw["fpocket_pockets"]["pocket_id"]) == {1, 2}


def test_renderer_uses_aligned_payload_and_zip_download():
    renderer = (ROOT / "workflow" / "builders" / "template" / "renderer.js").read_text()
    assert "PAY[fam].transforms" in renderer
    assert "alignedPdb(curFam,m)" in renderer
    assert "structures_zip_b64" in renderer
    assert "All structures (ZIP)" in renderer
    assert "All structures (multi-PDB)" not in renderer


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
    assert "FoldMason structure-guided alignment" in renderer
    assert "singletonTab(3)" in renderer
    assert "not applicable to a singleton" in renderer
    assert ".sequence-view" in prefix
    assert "seq=seq, msa=msa" in builder
