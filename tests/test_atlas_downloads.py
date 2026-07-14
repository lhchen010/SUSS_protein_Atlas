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
