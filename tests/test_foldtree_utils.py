from pathlib import Path
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from foldtree_utils import TreeValidationError, read_valid_tree, recover_metric


LEAVES = ["TDZ15773.1", "TDZ21235.1", "TDZ14607.1"]
PRE_ROOT = "(TDZ15773.1:0,TDZ21235.1:0.6794,TDZ14607.1:0.025);"


def _write(root, name, text):
    path = root / name
    path.write_text(text, encoding="utf-8")
    return path


def test_zero_length_branch_is_valid(tmp_path):
    tree_path = _write(tmp_path, "tree.nwk", PRE_ROOT)

    tree = read_valid_tree(tree_path, LEAVES)

    assert {leaf.name for leaf in tree.get_terminals()} == set(LEAVES)


def test_final_mad_tree_is_selected_before_raw_rooted(tmp_path):
    _write(tmp_path, "foldtree_struct_tree.PP.nwk.rooted.final", PRE_ROOT)
    _write(tmp_path, "foldtree_struct_tree.PP.nwk.rooted", "")
    destination = tmp_path / "F53_foldtree.nwk"

    status = recover_metric(tmp_path, "foldtree", destination, LEAVES)

    assert status["status"] == "complete"
    assert status["source_stage"] == "mad_final"
    assert status["rooting_method"] == "mad"
    read_valid_tree(destination, LEAVES)


def test_empty_mad_output_uses_midpoint_fallback(tmp_path):
    _write(tmp_path, "alntmscore_struct_tree.PP.nwk.rooted", "")
    _write(tmp_path, "alntmscore_struct_tree.PP.nwk", PRE_ROOT)
    destination = tmp_path / "F53_alntmscore.nwk"

    status = recover_metric(tmp_path, "alntmscore", destination, LEAVES)

    assert status["status"] == "complete_with_fallback"
    assert status["rooting_method"] == "midpoint"
    assert status["reason"] == "mad_unavailable_or_invalid"
    read_valid_tree(destination, LEAVES)


def test_small_family_policy_uses_midpoint_even_if_old_mad_tree_exists(tmp_path):
    _write(tmp_path, "lddt_struct_tree.PP.nwk.rooted.final", PRE_ROOT)
    _write(tmp_path, "lddt_struct_tree.PP.nwk", PRE_ROOT)
    destination = tmp_path / "F53_lddt.nwk"

    status = recover_metric(
        tmp_path, "lddt", destination, LEAVES, small_family=True
    )

    assert status["status"] == "complete"
    assert status["rooting_method"] == "midpoint"
    assert status["reason"] == "small_family_policy"


def test_leaf_set_mismatch_is_not_published(tmp_path):
    _write(tmp_path, "foldtree_struct_tree.PP.nwk", "(TDZ15773.1:1,WRONG:1);")
    destination = tmp_path / "F53_foldtree.nwk"

    with pytest.raises(TreeValidationError, match="terminal set mismatch"):
        recover_metric(tmp_path, "foldtree", destination, LEAVES)

    assert not destination.exists()


def test_negative_branch_length_is_rejected(tmp_path):
    tree_path = _write(
        tmp_path,
        "tree.nwk",
        "(TDZ15773.1:-0.1,TDZ21235.1:0.2,TDZ14607.1:0.3);",
    )

    with pytest.raises(TreeValidationError, match="invalid branch length"):
        read_valid_tree(tree_path, LEAVES)


def test_missing_newick_terminator_is_rejected(tmp_path):
    tree_path = _write(tmp_path, "tree.nwk", PRE_ROOT.rstrip(";"))

    with pytest.raises(TreeValidationError, match="terminator"):
        read_valid_tree(tree_path, LEAVES)
