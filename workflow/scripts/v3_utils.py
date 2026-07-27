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
    return aligned / total if total > 0 else math.nan


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
    data["pair"] = [
        tuple(sorted((query, target))) for query, target in zip(data.q, data.t)
    ]
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
        ["tm", "min_coverage", "evalue"], ascending=[False, False, True]
    ).groupby("pair", as_index=False).first()
    best[["q", "t"]] = pd.DataFrame(best.pair.tolist(), index=best.index)
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
        "lddt",
        "fident",
    ):
        data[column] = pd.to_numeric(data[column], errors="coerce")
    data["shorter_coverage"] = (
        data.alnlen / data[["qlen", "tlen"]].min(axis=1)
    ).clip(upper=1.0)
    data = data[
        (data.evalue <= float(evalue_threshold))
        & (data.prob >= float(probability_threshold))
        & (data.alnlen >= int(min_aligned_residues))
        & (data.shorter_coverage >= float(min_shorter_coverage))
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
                    "lddt",
                    "fident",
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
            "lddt",
            "fident",
            "alnlen",
            "shorter_coverage",
        ]
    ]


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
