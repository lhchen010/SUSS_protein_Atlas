"""Validation and recovery helpers for FoldTree Newick outputs."""

from __future__ import annotations

import io
import math
import os
import tempfile
from pathlib import Path

from Bio import Phylo


class TreeValidationError(ValueError):
    """Raised when a FoldTree output is not a usable family tree."""


def read_valid_tree(path, expected_leaves):
    """Read the first Newick record and validate its leaves and branch lengths."""
    path = Path(path)
    if not path.is_file() or path.stat().st_size == 0:
        raise TreeValidationError(f"{path}: missing or empty")

    text = path.read_text(encoding="utf-8", errors="replace").strip()
    terminator = text.find(";")
    if terminator < 0:
        raise TreeValidationError(f"{path}: Newick terminator is missing")

    try:
        tree = Phylo.read(io.StringIO(text[: terminator + 1]), "newick")
    except Exception as exc:
        raise TreeValidationError(f"{path}: invalid Newick ({exc})") from exc

    leaves = [str(leaf.name or "") for leaf in tree.get_terminals()]
    expected = [str(leaf) for leaf in expected_leaves]
    if len(leaves) != len(set(leaves)):
        raise TreeValidationError(f"{path}: duplicate terminal names")
    if set(leaves) != set(expected) or len(leaves) != len(expected):
        missing = sorted(set(expected) - set(leaves))
        extra = sorted(set(leaves) - set(expected))
        raise TreeValidationError(
            f"{path}: terminal set mismatch; missing={missing}, extra={extra}"
        )

    for clade in tree.find_clades():
        length = clade.branch_length
        if length is None:
            continue
        if not math.isfinite(float(length)) or float(length) < 0:
            raise TreeValidationError(
                f"{path}: invalid branch length {length!r} on {clade.name or 'internal node'}"
            )
    return tree


def write_tree_atomic(tree, destination):
    """Write one normalized Newick tree without exposing a partial output."""
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=".tmp", dir=destination.parent
    )
    os.close(fd)
    try:
        Phylo.write(tree, temporary, "newick")
        os.replace(temporary, destination)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def _candidate_paths(family_root, metric):
    prefix = Path(family_root) / f"{metric}_struct_tree.PP.nwk"
    return {
        "mad_final": Path(f"{prefix}.rooted.final"),
        "mad_rooted": Path(f"{prefix}.rooted"),
        "pre_root": prefix,
    }


def recover_metric(
    family_root,
    metric,
    destination,
    expected_leaves,
    *,
    small_family=False,
    fallback="midpoint",
):
    """Select, validate, and normalize one FoldTree metric output.

    Small families use the pre-root tree and midpoint rooting by policy. Larger
    families prefer MAD's finalized tree, then its raw rooted tree, and finally
    midpoint-root the valid pre-root tree when configured.
    """
    candidates = _candidate_paths(family_root, metric)
    errors = {}

    if not small_family:
        for stage in ("mad_final", "mad_rooted"):
            source = candidates[stage]
            try:
                tree = read_valid_tree(source, expected_leaves)
            except TreeValidationError as exc:
                errors[stage] = str(exc)
                continue
            write_tree_atomic(tree, destination)
            return {
                "status": "complete" if stage == "mad_final" else "complete_with_recovery",
                "rooting_method": "mad",
                "source_stage": stage,
                "source": str(source),
                "reason": None if stage == "mad_final" else "mad_final_missing_or_invalid",
            }

    if fallback == "midpoint" or small_family:
        source = candidates["pre_root"]
        try:
            tree = read_valid_tree(source, expected_leaves)
            tree.root_at_midpoint()
            write_tree_atomic(tree, destination)
            # Validate the serialized result too, including its exact terminal set.
            read_valid_tree(destination, expected_leaves)
        except Exception as exc:
            errors["pre_root"] = str(exc)
        else:
            return {
                "status": "complete" if small_family else "complete_with_fallback",
                "rooting_method": "midpoint",
                "source_stage": "pre_root",
                "source": str(source),
                "reason": "small_family_policy" if small_family else "mad_unavailable_or_invalid",
            }

    details = "; ".join(f"{stage}: {message}" for stage, message in errors.items())
    raise TreeValidationError(
        f"{metric}: no valid rooted tree or recoverable pre-root tree ({details})"
    )
