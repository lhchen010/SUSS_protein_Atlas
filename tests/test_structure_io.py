"""Tests for unified PDB/mmCIF reading, sequence extraction, and QC staging."""

import gzip
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

import runtime_utils
import structure_io
import validate_inputs


MINIMAL_MMCIF = """data_test
_entry.id test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
ATOM 1 C CA . ALA A 1 1 ? 0.000 0.000 0.000 1.00 85.50 1 ALA A CA 1
ATOM 2 C CA . CYS A 1 2 ? 3.800 0.000 0.000 1.00 92.50 2 CYS A CA 1
ATOM 3 C CA . GLY A 1 3 ? 7.600 0.000 0.000 1.00 78.00 3 GLY A CA 1
"""


def test_mmcif_metrics_sequence_and_staging_conversion(tmp_path):
    in_dir = tmp_path / "input_structs"
    out_dir = tmp_path / "qc_pdb"
    in_dir.mkdir()

    cif_path = in_dir / "bc_P12345.cif"
    cif_path.write_text(MINIMAL_MMCIF, encoding="utf-8")

    gz_cif_path = in_dir / "bc_Q9Y261.mmcif.gz"
    with gzip.open(gz_cif_path, "wt", encoding="utf-8") as fh:
        fh.write(MINIMAL_MMCIF)

    ca_count, bfactors = structure_io.read_ca_metrics(cif_path)
    assert ca_count == 3
    assert bfactors == [85.5, 92.5, 78.0]
    assert structure_io.read_ca_sequence(cif_path) == "ACG"
    assert structure_io.read_ca_sequence(gz_cif_path) == "ACG"

    # Preflight accepts mmCIF & gzipped mmCIF
    cfg = {
        "strain": {"code": "bc"},
        "input": {"pdb_dir": str(in_dir)},
        "qc": {"min_length": 3, "max_length": 1000, "min_plddt": 50},
    }
    errors, _, info = validate_inputs.validate(cfg)
    assert errors == []
    assert info["n_structures"] == 2
    assert info["n_mmcif_structures"] == 2

    # Staging converts mmCIF -> canonical .pdb preserving CA coords & pLDDT
    qc_csv = tmp_path / "qc.csv"
    qc_csv.write_text(
        "acc,fn,length,mean_plddt,pass\n"
        "P12345,bc_P12345.cif,3,85.33,True\n"
        "Q9Y261,bc_Q9Y261.mmcif.gz,3,85.33,True\n",
        encoding="utf-8",
    )
    staged = runtime_utils.stage_qc_structures(qc_csv, in_dir, out_dir)
    assert staged == ["bc_P12345.pdb", "bc_Q9Y261.pdb"]
    staged_ca, staged_bf = structure_io.read_ca_metrics(out_dir / "bc_P12345.pdb")
    assert staged_ca == 3
    assert staged_bf == [85.5, 92.5, 78.0]
    assert structure_io.read_ca_sequence(out_dir / "bc_P12345.pdb") == "ACG"
