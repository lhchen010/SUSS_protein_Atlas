#!/usr/bin/env python3
"""Compare a small F23 atlas run with direct tools and the full reference run."""

from __future__ import annotations

import argparse
import io
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from Bio import Phylo


FOLDSEEK_COLUMNS = [
    "query",
    "target",
    "alntmscore",
    "qtmscore",
    "ttmscore",
    "lddt",
    "fident",
    "alnlen",
    "qlen",
    "tlen",
    "evalue",
    "bits",
]
DOMAIN_COLUMNS = [
    "query",
    "target",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "alnlen",
    "qlen",
    "tlen",
    "evalue",
    "prob",
    "bits",
    "alntmscore",
    "lddt",
    "fident",
    "qcov",
    "tcov",
]
ACCESSION_RE = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")


def accession(value: object) -> str:
    match = ACCESSION_RE.search(str(value))
    return match.group(0) if match else str(value)


def canonical_pair(left: object, right: object) -> str:
    return "|".join(sorted((accession(left), accession(right))))


class Audit:
    def __init__(self) -> None:
        self.checks: list[dict[str, object]] = []

    def add(self, name: str, passed: bool, **details: object) -> None:
        self.checks.append({"name": name, "passed": bool(passed), **details})

    def write(self, path: Path) -> None:
        passed = sum(bool(item["passed"]) for item in self.checks)
        payload = {
            "status": "pass" if passed == len(self.checks) else "fail",
            "passed": passed,
            "failed": len(self.checks) - passed,
            "checks": self.checks,
        }
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def sorted_foldseek(path: Path, columns: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", names=columns)
    frame["query"] = frame["query"].map(accession)
    frame["target"] = frame["target"].map(accession)
    return frame.sort_values(["query", "target"]).reset_index(drop=True)


def frames_equal(left: pd.DataFrame, right: pd.DataFrame) -> bool:
    try:
        pd.testing.assert_frame_equal(
            left,
            right,
            check_dtype=False,
            check_like=True,
            rtol=1e-6,
            atol=1e-6,
        )
    except AssertionError:
        return False
    return True


def tree_distances(path: Path) -> tuple[list[str], dict[str, float]]:
    tree = Phylo.read(io.StringIO(path.read_text().strip()), "newick")
    leaves = sorted(terminal.name for terminal in tree.get_terminals())
    distances = {
        canonical_pair(left, right): float(tree.distance(left, right))
        for index, left in enumerate(leaves)
        for right in leaves[index + 1 :]
    }
    return leaves, distances


def compare_trees(audit: Audit, name: str, left: Path, right: Path) -> None:
    left_leaves, left_distances = tree_distances(left)
    right_leaves, right_distances = tree_distances(right)
    common = sorted(set(left_distances) & set(right_distances))
    max_delta = max(
        (abs(left_distances[key] - right_distances[key]) for key in common),
        default=float("inf"),
    )
    audit.add(
        name,
        left_leaves == right_leaves
        and set(left_distances) == set(right_distances)
        and max_delta <= 1e-8,
        n_leaves=len(left_leaves),
        n_pairwise_distances=len(common),
        max_patristic_delta=max_delta,
    )


def parse_usalign(path: Path) -> dict[str, float]:
    values: dict[str, float] = {}
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("DIRECT"):
            continue
        fields = line.split("\t")
        if len(fields) < 4:
            continue
        values[canonical_pair(fields[0], fields[1])] = min(
            float(fields[2]), float(fields[3])
        )
    return values


def matrix_pairs(path: Path) -> dict[str, float]:
    matrix = pd.read_csv(path, index_col=0)
    return {
        canonical_pair(left, right): float(matrix.loc[left, right])
        for index, left in enumerate(matrix.index)
        for right in matrix.index[index + 1 :]
    }


def parse_fpocket(info_path: Path) -> dict[str, object]:
    scores = {
        int(number): float(score)
        for number, score in re.findall(
            r"Pocket\s+(\d+)\s*:\s*\n\s*Score\s*:\s*([\-\d.]+)",
            info_path.read_text(),
        )
    }
    top_id = max(scores, key=scores.get)
    pocket_dir = info_path.parent / "pockets"
    residues: set[int] = set()
    for line in (pocket_dir / f"pocket{top_id}_atm.pdb").read_text().splitlines():
        if line.startswith(("ATOM", "HETATM")):
            try:
                residues.add(int(line[22:26]))
            except ValueError:
                pass
    return {
        "top_score": round(scores[top_id], 3),
        "n_pockets": len(scores),
        "lining_residues": sorted(residues),
    }


def parse_p2rank(path: Path) -> dict[str, object]:
    frame = pd.read_csv(path)
    frame.columns = [column.strip() for column in frame.columns]
    pockets = []
    for _, row in frame.iterrows():
        residues = sorted(
            {
                int(item.split("_")[-1])
                for item in str(row.get("residue_ids", "")).split()
                if item.split("_")[-1].isdigit()
            }
        )
        pockets.append(
            {
                "score": float(row.get("score", 0)),
                "probability": (
                    float(row.get("probability"))
                    if pd.notna(row.get("probability"))
                    else None
                ),
                "lining_residues": residues,
            }
        )
    if not pockets:
        return {
            "top_score": None,
            "top_probability": None,
            "n_pockets": 0,
            "lining_residues": [],
        }
    top = max(pockets, key=lambda pocket: pocket["score"])
    return {
        "top_score": top["score"],
        "top_probability": top["probability"],
        "n_pockets": len(pockets),
        "lining_residues": top["lining_residues"],
    }


def pocket_equal(left: dict[str, object], right: dict[str, object]) -> bool:
    scalar_keys = ("top_score", "top_probability")
    for key in scalar_keys:
        left_value = left.get(key)
        right_value = right.get(key)
        if left_value is None or right_value is None:
            if left_value != right_value:
                return False
        elif not np.isclose(float(left_value), float(right_value), atol=1e-9):
            return False
    return (
        int(left.get("n_pockets", -1)) == int(right.get("n_pockets", -2))
        and list(left.get("lining_residues", []))
        == list(right.get("lining_residues", []))
    )


def compare_direct_pockets(audit: Audit, root: Path, members: list[str]) -> None:
    atlas = json.loads((root / "results/pockets.json").read_text())
    fpocket_matches = 0
    p2rank_matches = 0
    for member in members:
        direct_root = root / "controls/pockets_direct" / member
        fpocket_info = next(
            (direct_root / "fpocket").glob("*_out/*_info.txt")
        )
        p2rank_csv = next(
            (direct_root / "p2rank/out").glob("*_predictions.csv")
        )
        if pocket_equal(atlas[member]["fpocket"], parse_fpocket(fpocket_info)):
            fpocket_matches += 1
        if pocket_equal(atlas[member]["p2rank"], parse_p2rank(p2rank_csv)):
            p2rank_matches += 1
    audit.add(
        "fpocket_direct_match",
        fpocket_matches == len(members),
        matched=fpocket_matches,
        expected=len(members),
    )
    audit.add(
        "p2rank_direct_match",
        p2rank_matches == len(members),
        matched=p2rank_matches,
        expected=len(members),
        profile="alphafold",
    )
    audit.add(
        "pocket_family_alias",
        "F0" in atlas and atlas["F0"].get("ref") == members[0],
        alias_ref=atlas.get("F0", {}).get("ref"),
        expected_ref=members[0],
    )


def raw_rnaseq_means(root: Path, members: list[str]) -> pd.DataFrame:
    mapping = pd.read_excel(root / "input/rnaseq.xlsx", sheet_name="id_mapping")
    expression = pd.read_excel(root / "input/rnaseq.xlsx", sheet_name="expression")
    accession_to_gene = dict(
        zip(
            mapping["protein_accession"].astype(str),
            mapping["gene_id"].astype(str),
        )
    )
    expression = expression.set_index(expression["gene_id"].astype(str))
    output = pd.read_csv(root / "results/rnaseq_expression.csv")
    rows = []
    for member in members:
        gene = accession_to_gene[member]
        row: dict[str, object] = {"acc": member}
        for condition in output.columns[1:]:
            columns = [
                column
                for column in expression.columns
                if column == condition or column.startswith(condition + ".")
            ]
            row[condition] = float(
                pd.to_numeric(expression.loc[gene, columns]).mean()
            )
        rows.append(row)
    return pd.DataFrame(rows)


def normalized_pair_frame(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output.insert(
        0,
        "pair",
        [
            canonical_pair(left, right)
            for left, right in zip(output["q"], output["t"])
        ],
    )
    return output.drop(columns=["q", "t"]).sort_values("pair").reset_index(drop=True)


def compare_reference(
    audit: Audit, root: Path, reference: Path, members: list[str]
) -> None:
    member_set = set(members)
    reference_members = pd.read_csv(reference / "results/members.csv")
    reference_f23 = set(
        reference_members.loc[
            reference_members["family"] == "F23", "acc"
        ].astype(str)
    )
    audit.add(
        "reference_f23_membership",
        member_set == reference_f23,
        audit_members=len(member_set),
        reference_members=len(reference_f23),
    )

    for filename, name in (
        ("edges.csv", "reference_f23_edges"),
        ("classification.csv", "reference_f23_classification"),
    ):
        small = pd.read_csv(root / "results" / filename)
        full = pd.read_csv(reference / "results" / filename)
        full = full[
            full["q"].map(accession).isin(member_set)
            & full["t"].map(accession).isin(member_set)
        ]
        left = normalized_pair_frame(small)
        right = normalized_pair_frame(full)
        if filename == "edges.csv":
            semantic_columns = [
                "pair",
                "tm",
                "fident",
                "alnlen",
                "evalue",
            ]
        else:
            semantic_columns = [
                "pair",
                "tm",
                "fident",
                "alnlen",
                "evalue",
                "pident",
                "blast_alignment_length",
                "blast_evalue",
                "blast_min_coverage",
                "blast_detected",
                "class",
            ]
        left = left[semantic_columns]
        right = right[semantic_columns]
        audit.add(
            name,
            frames_equal(left, right),
            audit_rows=len(left),
            reference_rows=len(right),
        )

    for small_name, reference_name, check_name in (
        ("F0_TM.csv", "F23_TM.csv", "reference_foldseek_tm_matrix"),
        ("F0_TM_usalign.csv", "F23_TM_usalign.csv", "reference_usalign_matrix"),
        (
            "F0_structural_conservation.csv",
            "F23_structural_conservation.csv",
            "reference_structural_conservation",
        ),
    ):
        left = pd.read_csv(root / "results/families/F0" / small_name)
        right = pd.read_csv(reference / "results/families/F23" / reference_name)
        if check_name == "reference_structural_conservation":
            left = left.drop(columns=["family"])
            right = right.drop(columns=["family"])
        audit.add(
            check_name,
            frames_equal(left, right),
            audit_rows=len(left),
            reference_rows=len(right),
        )

    for suffix, name in (
        ("fm_aa.fa", "reference_foldmason_aa"),
        ("fm_3di.fa", "reference_foldmason_3di"),
    ):
        left = (root / "results/families/F0" / f"F0_{suffix}").read_bytes()
        right = (
            reference / "results/families/F23" / f"F23_{suffix}"
        ).read_bytes()
        audit.add(name, left == right, bytes=len(left))

    for metric in ("foldtree", "alntmscore", "lddt"):
        compare_trees(
            audit,
            f"reference_{metric}_tree",
            root / "results/families/F0" / f"F0_{metric}.nwk",
            reference / "results/families/F23" / f"F23_{metric}.nwk",
        )

    left_expression = pd.read_csv(root / "results/rnaseq_expression.csv")
    right_expression = pd.read_csv(reference / "results/rnaseq_expression.csv")
    right_expression = right_expression[
        right_expression["acc"].astype(str).isin(member_set)
    ]
    left_expression = left_expression.sort_values("acc").reset_index(drop=True)
    right_expression = right_expression.sort_values("acc").reset_index(drop=True)
    audit.add(
        "reference_rnaseq",
        frames_equal(left_expression, right_expression),
        rows=len(left_expression),
        conditions=len(left_expression.columns) - 1,
    )

    left_annotation = pd.read_csv(root / "results/member_annotation.csv")
    right_annotation = pd.read_csv(reference / "results/member_annotation.csv")
    right_annotation = right_annotation[
        right_annotation["acc"].astype(str).isin(member_set)
    ]
    common = sorted(
        (set(left_annotation.columns) & set(right_annotation.columns))
        - {"family"}
    )
    left_annotation = (
        left_annotation[common]
        .sort_values("acc")
        .reset_index(drop=True)
        .fillna("")
    )
    right_annotation = (
        right_annotation[common]
        .sort_values("acc")
        .reset_index(drop=True)
        .fillna("")
    )
    audit.add(
        "reference_annotation",
        frames_equal(left_annotation, right_annotation),
        rows=len(left_annotation),
        compared_columns=len(common),
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("audit_root", type=Path)
    parser.add_argument("reference_root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.audit_root.resolve()
    reference = args.reference_root.resolve()
    output = args.output.resolve()
    audit = Audit()

    members = [
        accession(line)
        for line in (root / "results/family_members/F0.members.txt")
        .read_text()
        .splitlines()
        if line.strip()
    ]
    audit.add(
        "cohort_shape",
        len(members) == 10 and len(set(members)) == 10,
        members=len(members),
    )

    for atlas_path, direct_path, columns, name in (
        (
            root / "results/foldseek_allvsall.tsv",
            root / "controls/foldseek_full_direct.tsv",
            FOLDSEEK_COLUMNS,
            "foldseek_full_direct_match",
        ),
        (
            root / "results/foldseek_domain_allvsall.tsv",
            root / "controls/foldseek_domain_direct.tsv",
            DOMAIN_COLUMNS,
            "foldseek_domain_direct_match",
        ),
    ):
        left = sorted_foldseek(atlas_path, columns)
        right = sorted_foldseek(direct_path, columns)
        audit.add(
            name,
            frames_equal(left, right),
            atlas_rows=len(left),
            direct_rows=len(right),
        )

    for atlas_path, direct_path, name in (
        (
            root / "results/families/F0/F0_fm_aa.fa",
            root / "controls/F0_direct_aa.fa",
            "foldmason_aa_direct_match",
        ),
        (
            root / "results/families/F0/F0_fm_3di.fa",
            root / "controls/F0_direct_3di.fa",
            "foldmason_3di_direct_match",
        ),
        (
            root / "results/families/F0/F0_fm.nw",
            root / "controls/F0_direct.nw",
            "foldmason_guide_tree_direct_match",
        ),
    ):
        audit.add(
            name,
            atlas_path.read_bytes() == direct_path.read_bytes(),
            bytes=atlas_path.stat().st_size,
        )

    direct_usalign = parse_usalign(root / "controls/usalign_direct.tsv")
    atlas_usalign = matrix_pairs(
        root / "results/families/F0/F0_TM_usalign.csv"
    )
    common_pairs = sorted(set(direct_usalign) & set(atlas_usalign))
    max_usalign_delta = max(
        (
            abs(direct_usalign[pair] - atlas_usalign[pair])
            for pair in common_pairs
        ),
        default=float("inf"),
    )
    audit.add(
        "usalign_direct_match",
        len(common_pairs) == 45
        and set(direct_usalign) == set(atlas_usalign)
        and max_usalign_delta <= 0.00051,
        pairs=len(common_pairs),
        max_absolute_delta_before_display_rounding=max_usalign_delta,
    )

    for metric in ("foldtree", "alntmscore", "lddt"):
        compare_trees(
            audit,
            f"{metric}_tree_direct_match",
            root / "results/families/F0" / f"F0_{metric}.nwk",
            root
            / "controls/foldtree_direct/CTRL_F0"
            / f"CTRL_F0_{metric}.nwk",
        )

    compare_direct_pockets(audit, root, members)

    raw_means = raw_rnaseq_means(root, members).sort_values("acc").reset_index(
        drop=True
    )
    atlas_means = (
        pd.read_csv(root / "results/rnaseq_expression.csv")
        .sort_values("acc")
        .reset_index(drop=True)
    )
    numeric_delta = (
        raw_means.iloc[:, 1:].to_numpy(dtype=float)
        - atlas_means.iloc[:, 1:].to_numpy(dtype=float)
    )
    audit.add(
        "rnaseq_raw_replicate_means",
        raw_means["acc"].tolist() == atlas_means["acc"].tolist()
        and float(np.max(np.abs(numeric_delta))) <= 1e-10,
        proteins=len(raw_means),
        conditions=len(raw_means.columns) - 1,
        values=numeric_delta.size,
        max_absolute_delta=float(np.max(np.abs(numeric_delta))),
    )

    annotation = pd.read_csv(root / "results/member_annotation.csv")
    status_columns = [
        "annotation_status",
        "interpro_status",
        "foldseek_pdb_status",
        "foldseek_afdb_status",
        "effectorp_status",
        "deeptmhmm_status",
    ]
    complete = int(
        annotation[status_columns]
        .eq("complete")
        .all(axis=1)
        .sum()
    )
    audit.add(
        "annotation_complete",
        complete == len(members),
        complete=complete,
        expected=len(members),
    )
    audit.add(
        "atlas_html_generated",
        (root / "results/f23_audit_atlas.html").stat().st_size > 100_000,
        bytes=(root / "results/f23_audit_atlas.html").stat().st_size,
    )

    compare_reference(audit, root, reference, members)
    audit.write(output)
    failed = sum(not bool(check["passed"]) for check in audit.checks)
    print(
        f"audit checks: {len(audit.checks) - failed}/{len(audit.checks)} "
        f"passed; report={output}"
    )
    for check in audit.checks:
        print(("PASS" if check["passed"] else "FAIL"), check["name"])
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
