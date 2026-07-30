"""Pure helpers for the v3 sequence/structure analysis layers."""

from __future__ import annotations

import math
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from runtime_utils import symmetric_tm


ACCESSION_RE = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")


def protein_id(value: object) -> str:
    """Return a stable protein identifier from a Foldseek/FASTA/PDB label."""
    text = Path(str(value)).name
    match = ACCESSION_RE.search(text)
    if match:
        return match.group(0)
    text = re.sub(r"\.(pdb|cif|mmcif)(\.gz)?$", "", text, flags=re.IGNORECASE)
    return text.split()[0]


def coverage(alnlen: object, length: object) -> float:
    """Alignment coverage, guarded against missing and zero lengths."""
    try:
        aligned = float(alnlen)
        total = float(length)
    except (TypeError, ValueError):
        return math.nan
    if total <= 0:
        return math.nan
    return min(1.0, max(0.0, aligned / total))


def _canonicalize_pair_direction(
    data: pd.DataFrame,
    swaps: tuple[tuple[str, str], ...] = (),
) -> pd.DataFrame:
    """Orient pair rows lexically and keep directional fields attached to IDs."""
    output = data.copy()
    reverse = output["q"].astype(str) > output["t"].astype(str)
    if not reverse.any():
        return output
    identifiers = output.loc[reverse, ["q", "t"]].copy()
    output.loc[reverse, "q"] = identifiers["t"].to_numpy()
    output.loc[reverse, "t"] = identifiers["q"].to_numpy()
    for left, right in swaps:
        if left not in output.columns or right not in output.columns:
            continue
        values = output.loc[reverse, [left, right]].copy()
        output.loc[reverse, left] = values[right].to_numpy()
        output.loc[reverse, right] = values[left].to_numpy()
    return output


def select_blast_relationships(
    table: pd.DataFrame,
    *,
    evalue_threshold: float,
    coverage_threshold: float,
) -> pd.DataFrame:
    """Select the best HSP that passes both relationship thresholds."""
    data = table.copy()
    if data.empty:
        data["pair"] = pd.Series(dtype=object)
        return data
    data["evalue"] = pd.to_numeric(data["evalue"], errors="coerce")
    data["bitscore"] = pd.to_numeric(data["bitscore"], errors="coerce")
    data["min_coverage"] = pd.to_numeric(
        data["min_coverage"], errors="coerce"
    )
    data = data[
        (data.evalue <= float(evalue_threshold))
        & (data.min_coverage >= float(coverage_threshold))
    ].copy()
    data = _canonicalize_pair_direction(
        data,
        swaps=(("qcov", "scov"), ("qlen", "slen")),
    )
    data["pair"] = list(zip(data.q.astype(str), data.t.astype(str)))
    tie_columns = [
        column
        for column in ("pident", "min_coverage", "q", "t")
        if column in data.columns
    ]
    ascending = [True, False] + [
        False if column in {"pident", "min_coverage"} else True
        for column in tie_columns
    ]
    return (
        data.sort_values(
            ["evalue", "bitscore", *tie_columns],
            ascending=ascending,
        )
        .drop_duplicates("pair", keep="first")
        .reset_index(drop=True)
    )


def domain_sequence_records(
    members: pd.DataFrame,
    parent_sequences: dict[str, str],
) -> tuple[dict[str, str], pd.DataFrame]:
    """Extract coordinate-defined domain sequences with an auditable manifest."""
    records: dict[str, str] = {}
    rows = []
    required = {"domain_family", "segment_id", "acc", "start", "end"}
    missing = required - set(members.columns)
    if missing:
        raise ValueError(
            "domain member table is missing columns: " + ", ".join(sorted(missing))
        )
    ordered = members.sort_values(
        ["domain_family", "segment_id"], kind="stable"
    )
    for row in ordered.itertuples():
        family = str(row.domain_family)
        segment_id = str(row.segment_id)
        accession = str(row.acc)
        start, end = int(row.start), int(row.end)
        sequence = parent_sequences.get(accession, "")
        if not sequence:
            raise ValueError(f"{segment_id}: parent sequence {accession} is missing")
        if start < 1 or end < start or end > len(sequence):
            raise ValueError(
                f"{segment_id}: invalid coordinates {start}-{end} "
                f"for {accession} length {len(sequence)}"
            )
        segment = sequence[start - 1:end]
        if segment_id in records and records[segment_id] != segment:
            raise ValueError(f"{segment_id}: duplicate identifier has different sequence")
        records[segment_id] = segment
        rows.append(
            {
                "domain_family": family,
                "segment_id": segment_id,
                "acc": accession,
                "start": start,
                "end": end,
                "segment_length": len(segment),
                "parent_length": len(sequence),
            }
        )
    return records, pd.DataFrame(
        rows,
        columns=[
            "domain_family",
            "segment_id",
            "acc",
            "start",
            "end",
            "segment_length",
            "parent_length",
        ],
    )


def relationship_components(
    labels: list[str],
    relationships: pd.DataFrame,
) -> list[list[str]]:
    """Connected components from exact segment-to-segment relationships."""
    ordered = list(dict.fromkeys(map(str, labels)))
    parent = {label: label for label in ordered}

    def find(item: str) -> str:
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left: str, right: str) -> None:
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    allowed = set(ordered)
    if relationships is not None and len(relationships):
        for row in relationships.itertuples():
            left, right = str(row.q), str(row.t)
            if left in allowed and right in allowed and left != right:
                union(left, right)
    groups: dict[str, list[str]] = defaultdict(list)
    for label in ordered:
        groups[find(label)].append(label)
    return sorted(
        (sorted(group) for group in groups.values()),
        key=lambda group: (-len(group), group[0]),
    )


def blast_identity_matrix(
    table: pd.DataFrame,
    labels: list[str],
) -> list[list[float]]:
    """Best-HSP identity matrix for exact domain-segment identifiers."""
    ordered = list(dict.fromkeys(map(str, labels)))
    index = {label: position for position, label in enumerate(ordered)}
    matrix = np.zeros((len(ordered), len(ordered)), dtype=float)
    np.fill_diagonal(matrix, 1.0)
    if table is None or not len(table):
        return matrix.tolist()
    for row in table.itertuples():
        left, right = str(row.q), str(row.t)
        if left not in index or right not in index:
            continue
        try:
            identity = float(row.pident) / 100.0
        except (TypeError, ValueError):
            continue
        if not math.isfinite(identity):
            continue
        left_index, right_index = index[left], index[right]
        value = max(matrix[left_index, right_index], identity)
        matrix[left_index, right_index] = value
        matrix[right_index, left_index] = value
    return matrix.round(6).tolist()


def whole_fold_edges(
    table: pd.DataFrame,
    *,
    tm_threshold: float,
    coverage_threshold: float,
    symmetry: str,
) -> pd.DataFrame:
    """Collapse directed Foldseek hits into auditable undirected whole-fold edges."""
    data = table.copy()
    data["q"] = data["query"].map(protein_id)
    data["t"] = data["target"].map(protein_id)
    data = data[data.q != data.t].copy()
    data["qcov"] = [coverage(a, n) for a, n in zip(data.alnlen, data.qlen)]
    data["tcov"] = [coverage(a, n) for a, n in zip(data.alnlen, data.tlen)]
    data["min_coverage"] = data[["qcov", "tcov"]].min(axis=1)
    data["tm"] = [
        symmetric_tm(float(q), float(t), symmetry)
        for q, t in zip(data.qtmscore, data.ttmscore)
    ]
    data = data[
        (data.tm >= float(tm_threshold))
        & (data.min_coverage >= float(coverage_threshold))
    ].copy()
    if data.empty:
        return pd.DataFrame(
            columns=[
                "q",
                "t",
                "tm",
                "qtmscore",
                "ttmscore",
                "qcov",
                "tcov",
                "min_coverage",
                "lddt",
                "fident",
                "alnlen",
                "evalue",
            ]
        )
    data = _canonicalize_pair_direction(
        data,
        swaps=(("qtmscore", "ttmscore"), ("qcov", "tcov")),
    )
    data["pair"] = list(zip(data.q, data.t))
    numeric = [
        "tm",
        "qtmscore",
        "ttmscore",
        "qcov",
        "tcov",
        "min_coverage",
        "lddt",
        "fident",
        "alnlen",
    ]
    for column in numeric:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data["evalue"] = pd.to_numeric(data["evalue"], errors="coerce")
    best = data.sort_values(
        [
            "tm",
            "min_coverage",
            "evalue",
            "lddt",
            "fident",
            "alnlen",
            "q",
            "t",
        ],
        ascending=[False, False, True, False, False, False, True, True],
    ).drop_duplicates("pair", keep="first").reset_index(drop=True)
    return best[
        [
            "q",
            "t",
            "tm",
            "qtmscore",
            "ttmscore",
            "qcov",
            "tcov",
            "min_coverage",
            "lddt",
            "fident",
            "alnlen",
            "evalue",
        ]
    ]


def merge_intervals(
    intervals: list[tuple[int, int]], reciprocal_overlap: float = 0.5
) -> list[tuple[int, int]]:
    """Merge strongly overlapping local-hit intervals for one protein."""
    clean = sorted((min(a, b), max(a, b)) for a, b in intervals if a and b)
    groups: list[list[tuple[int, int]]] = []
    for interval in clean:
        assigned = False
        for group in groups:
            start = min(x[0] for x in group)
            end = max(x[1] for x in group)
            overlap = max(0, min(end, interval[1]) - max(start, interval[0]) + 1)
            shorter = min(end - start + 1, interval[1] - interval[0] + 1)
            if shorter and overlap / shorter >= reciprocal_overlap:
                group.append(interval)
                assigned = True
                break
        if not assigned:
            groups.append([interval])
    return [
        (
            int(round(np.median([item[0] for item in group]))),
            int(round(np.median([item[1] for item in group]))),
        )
        for group in groups
    ]


def domain_segments(
    table: pd.DataFrame,
    *,
    evalue_threshold: float,
    probability_threshold: float,
    min_aligned_residues: int,
    min_shorter_coverage: float,
    min_lddt: float = 0.0,
    min_alntm: float = 0.0,
    interval_overlap: float = 0.5,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert local Foldseek hits into segment nodes and segment-to-segment edges."""
    data = table.copy()
    data["q"] = data["query"].map(protein_id)
    data["t"] = data["target"].map(protein_id)
    data = data[data.q != data.t].copy()
    for column in (
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
    ):
        if column not in data:
            data[column] = 0.0
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data["shorter_coverage"] = (
        data.alnlen / data[["qlen", "tlen"]].min(axis=1)
    ).clip(upper=1.0)
    data = data[
        (data.evalue <= float(evalue_threshold))
        & (data.prob >= float(probability_threshold))
        & (data.alnlen >= int(min_aligned_residues))
        & (data.shorter_coverage >= float(min_shorter_coverage))
        & (data.lddt >= float(min_lddt))
        & (data.alntmscore >= float(min_alntm))
    ].copy()
    if data.empty:
        return (
            pd.DataFrame(columns=["segment_id", "acc", "start", "end", "length"]),
            pd.DataFrame(
                columns=[
                    "source",
                    "target",
                    "evalue",
                    "prob",
                    "bits",
                    "alntmscore",
                    "lddt",
                    "fident",
                    "qcov",
                    "tcov",
                    "alnlen",
                    "shorter_coverage",
                ]
            ),
        )

    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    for row in data.itertuples():
        intervals[row.q].append((int(row.qstart), int(row.qend)))
        intervals[row.t].append((int(row.tstart), int(row.tend)))
    consolidated: dict[str, list[tuple[int, int]]] = {
        acc: merge_intervals(items, interval_overlap)
        for acc, items in intervals.items()
    }
    segment_rows = []
    for acc, items in consolidated.items():
        for start, end in items:
            segment_rows.append(
                {
                    "segment_id": f"{acc}:{start}-{end}",
                    "acc": acc,
                    "start": start,
                    "end": end,
                    "length": end - start + 1,
                }
            )
    segments = pd.DataFrame(segment_rows)

    def nearest(acc: str, start: int, end: int) -> str:
        candidates = consolidated[acc]
        chosen = max(
            candidates,
            key=lambda item: max(
                0, min(item[1], end) - max(item[0], start) + 1
            ),
        )
        return f"{acc}:{chosen[0]}-{chosen[1]}"

    data["source"] = [
        nearest(acc, int(start), int(end))
        for acc, start, end in zip(data.q, data.qstart, data.qend)
    ]
    data["target"] = [
        nearest(acc, int(start), int(end))
        for acc, start, end in zip(data.t, data.tstart, data.tend)
    ]
    data = data[data.source != data.target].copy()
    data["pair"] = [
        tuple(sorted((source, target)))
        for source, target in zip(data.source, data.target)
    ]
    best = data.sort_values(
        ["evalue", "bits", "prob"], ascending=[True, False, False]
    ).groupby("pair", as_index=False).first()
    best[["source", "target"]] = pd.DataFrame(
        best.pair.tolist(), index=best.index
    )
    return segments, best[
        [
            "source",
            "target",
            "evalue",
            "prob",
            "bits",
            "alntmscore",
            "lddt",
            "fident",
            "qcov",
            "tcov",
            "alnlen",
            "shorter_coverage",
        ]
    ]


def domain_match_diagnostics(
    table: pd.DataFrame,
    accessions: list[str],
    assigned_accessions: set[str],
    *,
    evalue_threshold: float,
    probability_threshold: float,
    min_aligned_residues: int,
    min_shorter_coverage: float,
    min_lddt: float,
    min_alntm: float,
    borderline_lddt_margin: float = 0.05,
) -> pd.DataFrame:
    """Report the best local hit and exact filter outcome for every protein.

    This diagnostic does not alter D-family membership. A borderline hit must
    pass every configured filter except local lDDT and fall within the lDDT
    margin.
    """
    output_columns = [
        "acc",
        "status",
        "best_match",
        "protein_start",
        "protein_end",
        "match_start",
        "match_end",
        "protein_length",
        "match_length",
        "evalue",
        "prob",
        "bits",
        "alntmscore",
        "lddt",
        "fident",
        "protein_coverage",
        "match_coverage",
        "shorter_coverage",
        "alnlen",
        "raw_hit_count",
        "passed_filter_count",
        "passes_all",
        "failed_filters",
        "summary",
    ]
    ordered_accessions = list(dict.fromkeys(map(str, accessions)))
    data = table.copy()
    for column in (
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
    ):
        if column not in data:
            data[column] = np.nan
    if data.empty:
        return pd.DataFrame(
            [
                {
                    "acc": accession,
                    "status": (
                        "retained_family"
                        if accession in assigned_accessions
                        else "no_raw_hit"
                    ),
                    "raw_hit_count": 0,
                    "passed_filter_count": 0,
                    "passes_all": False,
                    "failed_filters": "",
                    "summary": (
                        "Protein has a retained D-family segment."
                        if accession in assigned_accessions
                        else "Foldseek reported no non-self local hit."
                    ),
                }
                for accession in ordered_accessions
            ],
            columns=output_columns,
        )

    data["q"] = data["query"].map(protein_id)
    data["t"] = data["target"].map(protein_id)
    data = data[data.q != data.t].copy()
    numeric_columns = [
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
    for column in numeric_columns:
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data["shorter_coverage"] = (
        data.alnlen / data[["qlen", "tlen"]].min(axis=1)
    ).clip(lower=0.0, upper=1.0)

    oriented_rows = []
    for row in data.itertuples():
        common = {
            "evalue": row.evalue,
            "prob": row.prob,
            "bits": row.bits,
            "alntmscore": row.alntmscore,
            "lddt": row.lddt,
            "fident": row.fident,
            "shorter_coverage": row.shorter_coverage,
            "alnlen": row.alnlen,
        }
        oriented_rows.extend(
            [
                {
                    **common,
                    "acc": row.q,
                    "best_match": row.t,
                    "protein_start": row.qstart,
                    "protein_end": row.qend,
                    "match_start": row.tstart,
                    "match_end": row.tend,
                    "protein_length": row.qlen,
                    "match_length": row.tlen,
                    "protein_coverage": row.qcov,
                    "match_coverage": row.tcov,
                },
                {
                    **common,
                    "acc": row.t,
                    "best_match": row.q,
                    "protein_start": row.tstart,
                    "protein_end": row.tend,
                    "match_start": row.qstart,
                    "match_end": row.qend,
                    "protein_length": row.tlen,
                    "match_length": row.qlen,
                    "protein_coverage": row.tcov,
                    "match_coverage": row.qcov,
                },
            ]
        )
    oriented = pd.DataFrame(oriented_rows)
    if oriented.empty:
        return domain_match_diagnostics(
            pd.DataFrame(),
            ordered_accessions,
            assigned_accessions,
            evalue_threshold=evalue_threshold,
            probability_threshold=probability_threshold,
            min_aligned_residues=min_aligned_residues,
            min_shorter_coverage=min_shorter_coverage,
            min_lddt=min_lddt,
            min_alntm=min_alntm,
            borderline_lddt_margin=borderline_lddt_margin,
        )

    checks = [
        ("evalue", oriented.evalue <= float(evalue_threshold)),
        ("probability", oriented.prob >= float(probability_threshold)),
        ("aligned_residues", oriented.alnlen >= int(min_aligned_residues)),
        (
            "shorter_coverage",
            oriented.shorter_coverage >= float(min_shorter_coverage),
        ),
        ("lddt", oriented.lddt >= float(min_lddt)),
        ("alntm", oriented.alntmscore >= float(min_alntm)),
    ]
    for name, values in checks:
        oriented[f"pass_{name}"] = values.fillna(False)
    pass_columns = [f"pass_{name}" for name, _ in checks]
    oriented["passed_filter_count"] = oriented[pass_columns].sum(axis=1)
    oriented["passes_all"] = oriented[pass_columns].all(axis=1)
    non_lddt_columns = [column for column in pass_columns if column != "pass_lddt"]
    oriented["borderline_lddt"] = (
        oriented[non_lddt_columns].all(axis=1)
        & ~oriented.pass_lddt
        & (
            oriented.lddt
            >= float(min_lddt) - max(0.0, float(borderline_lddt_margin))
        )
    )

    def number(value: object, digits: int = 3) -> str:
        try:
            value = float(value)
        except (TypeError, ValueError):
            return "NA"
        return f"{value:.{digits}g}" if math.isfinite(value) else "NA"

    def failed_filters(row: pd.Series) -> str:
        failures = []
        if not row.pass_evalue:
            failures.append(
                f"e-value {number(row.evalue)} > {number(evalue_threshold)}"
            )
        if not row.pass_probability:
            failures.append(
                "probability "
                f"{number(row.prob)} < {number(probability_threshold)}"
            )
        if not row.pass_aligned_residues:
            failures.append(
                f"aligned residues {number(row.alnlen, 6)} "
                f"< {int(min_aligned_residues)}"
            )
        if not row.pass_shorter_coverage:
            failures.append(
                "shorter coverage "
                f"{number(row.shorter_coverage)} "
                f"< {number(min_shorter_coverage)}"
            )
        if not row.pass_lddt:
            failures.append(
                f"local lDDT {number(row.lddt)} < {number(min_lddt)}"
            )
        if not row.pass_alntm:
            failures.append(
                f"alignment TM {number(row.alntmscore)} < {number(min_alntm)}"
            )
        return "; ".join(failures)

    oriented["failed_filters"] = oriented.apply(failed_filters, axis=1)
    oriented = oriented.sort_values(
        [
            "acc",
            "passes_all",
            "borderline_lddt",
            "passed_filter_count",
            "lddt",
            "alntmscore",
            "prob",
            "alnlen",
            "evalue",
        ],
        ascending=[True, False, False, False, False, False, False, False, True],
        kind="stable",
    )
    best_by_accession = {
        str(row.acc): row
        for row in oriented.drop_duplicates("acc", keep="first").itertuples()
    }
    hit_counts = oriented.groupby("acc").size().to_dict()

    records = []
    for accession in ordered_accessions:
        best = best_by_accession.get(accession)
        if best is None:
            status = (
                "retained_family"
                if accession in assigned_accessions
                else "no_raw_hit"
            )
            records.append(
                {
                    "acc": accession,
                    "status": status,
                    "raw_hit_count": 0,
                    "passed_filter_count": 0,
                    "passes_all": False,
                    "failed_filters": "",
                    "summary": (
                        "Protein has a retained D-family segment."
                        if status == "retained_family"
                        else "Foldseek reported no non-self local hit."
                    ),
                }
            )
            continue
        if accession in assigned_accessions:
            status = "retained_family"
            summary = "Protein has a retained D-family segment."
        elif bool(best.passes_all):
            status = "unclustered_retained_hit"
            summary = (
                "A local hit passed the edge filters but did not produce a retained "
                "D family under the community and minimum-size settings."
            )
        elif bool(best.borderline_lddt):
            status = "borderline"
            summary = (
                "Best local hit passed every filter except local lDDT and is within "
                f"{float(borderline_lddt_margin):.3f} of that cutoff."
            )
        else:
            status = "filtered_hit"
            summary = "Foldseek reported a local hit that failed one or more filters."
        record = {
            column: getattr(best, column)
            for column in output_columns
            if hasattr(best, column)
        }
        record.update(
            {
                "acc": accession,
                "status": status,
                "raw_hit_count": int(hit_counts.get(accession, 0)),
                "passed_filter_count": int(best.passed_filter_count),
                "passes_all": bool(best.passes_all),
                "failed_filters": str(best.failed_filters),
                "summary": summary,
            }
        )
        records.append(record)
    return pd.DataFrame(records, columns=output_columns)


def aggregate_domain_bridges(
    edges: pd.DataFrame, segment_family: dict[str, str]
) -> pd.DataFrame:
    """Aggregate local Foldseek edges that connect two different D families."""
    columns = [
        "source_family",
        "target_family",
        "n_edges",
        "mean_probability",
        "max_probability",
        "mean_lddt",
        "max_lddt",
        "mean_aligned_residues",
    ]
    if edges.empty:
        return pd.DataFrame(columns=columns)
    data = edges.copy()
    data["source_family"] = data.source.map(segment_family)
    data["target_family"] = data.target.map(segment_family)
    data = data[
        data.source_family.notna()
        & data.target_family.notna()
        & (data.source_family != data.target_family)
    ].copy()
    if data.empty:
        return pd.DataFrame(columns=columns)
    ordered = [
        tuple(sorted((source, target)))
        for source, target in zip(data.source_family, data.target_family)
    ]
    data[["source_family", "target_family"]] = pd.DataFrame(
        ordered, index=data.index
    )
    grouped = data.groupby(["source_family", "target_family"], as_index=False)
    result = grouped.agg(
        n_edges=("source", "size"),
        mean_probability=("prob", "mean"),
        max_probability=("prob", "max"),
        mean_lddt=("lddt", "mean"),
        max_lddt=("lddt", "max"),
        mean_aligned_residues=("alnlen", "mean"),
    )
    return result[columns].sort_values(
        ["n_edges", "mean_lddt"], ascending=[False, False]
    )


def foldmason_column_scores(payload: dict[str, Any], length: int) -> list[float]:
    """Validate and normalize FoldMason msa2lddtjson column scores."""
    raw = payload.get("scores")
    if not isinstance(raw, list):
        raise ValueError("FoldMason LDDT JSON does not contain a scores list")
    if len(raw) != int(length):
        raise ValueError(
            f"FoldMason returned {len(raw)} scores for an alignment of {length} columns"
        )
    scores = []
    for value in raw:
        try:
            score = float(value)
        except (TypeError, ValueError):
            score = math.nan
        scores.append(score if 0.0 <= score <= 1.0 else math.nan)
    return scores


def fasta_records(path: str | Path) -> dict[str, str]:
    """Read unaligned or aligned FASTA records keyed by normalized protein id."""
    records: dict[str, str] = {}
    key: str | None = None
    chunks: list[str] = []
    for raw in Path(path).read_text(encoding="utf-8", errors="replace").splitlines():
        line = raw.strip()
        if line.startswith(">"):
            if key is not None:
                records[key] = "".join(chunks)
            key = protein_id(line[1:])
            chunks = []
        elif key is not None:
            chunks.append(line)
    if key is not None:
        records[key] = "".join(chunks)
    return records


def write_fasta(records: dict[str, str], path: str | Path) -> None:
    """Write deterministic wrapped FASTA."""
    with Path(path).open("w", encoding="utf-8") as handle:
        for key, sequence in records.items():
            handle.write(f">{key}\n")
            for offset in range(0, len(sequence), 80):
                handle.write(sequence[offset : offset + 80] + "\n")


def shannon(values: list[str]) -> float:
    """Shannon entropy in bits for one alignment column, excluding gaps."""
    clean = [value for value in values if value not in {"-", ".", "X", "x"}]
    if not clean:
        return math.nan
    counts = Counter(clean)
    total = len(clean)
    return -sum((count / total) * math.log2(count / total) for count in counts.values())
