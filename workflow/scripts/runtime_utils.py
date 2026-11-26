"""Shared runtime helpers for external tools and symmetric TM scores."""

from __future__ import annotations

import csv
import os
import shutil
from pathlib import Path
from shutil import which


VALID_TM_SYMMETRY = {"min", "max", "mean"}


def resolve_executable(value: str | os.PathLike | None, name: str, *, required: bool = True) -> str | None:
    """Resolve an executable supplied as a path or a command on PATH."""
    raw = str(value or "").strip()
    if raw:
        candidate = Path(raw).expanduser()
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate.resolve())
        resolved = which(raw)
        if resolved:
            return str(Path(resolved).resolve())
    if required:
        raise FileNotFoundError(f"Required executable not found for {name}: {raw or '<empty>'}")
    return None


def resolve_file(value: str | os.PathLike | None, name: str, *, required: bool = True) -> str | None:
    """Resolve a configured non-executable file, such as a Python entry script."""
    raw = str(value or "").strip()
    if raw:
        candidate = Path(raw).expanduser()
        if candidate.is_file():
            return str(candidate.resolve())
    if required:
        raise FileNotFoundError(f"Required file not found for {name}: {raw or '<empty>'}")
    return None


def symmetric_tm(qtm: float, ttm: float, mode: str) -> float:
    """Combine directional TM scores using the configured convention."""
    if mode not in VALID_TM_SYMMETRY:
        raise ValueError(f"tm_symmetric must be one of {sorted(VALID_TM_SYMMETRY)}; got {mode!r}")
    if mode == "min":
        return min(qtm, ttm)
    if mode == "max":
        return max(qtm, ttm)
    return (qtm + ttm) / 2.0


def analysis_status(enabled: bool, statuses: dict[str, str], required: tuple[str, ...]) -> str:
    """Return not_run, complete, or partial from component states."""
    if not enabled:
        return "not_run"
    return "complete" if all(statuses.get(key) == "complete" for key in required) else "partial"


def novel_call(has_domain: bool, has_known_fold: bool, evidence_complete: bool) -> bool | None:
    """Only make a novelty call when both domain and fold evidence completed."""
    if not evidence_complete:
        return None
    return not has_domain and not has_known_fold


def stage_qc_structures(
    qc_csv: str | os.PathLike,
    pdb_dir: str | os.PathLike,
    output_dir: str | os.PathLike,
) -> list[str]:
    """Stage QC-passing PDB/mmCIF files into a canonical PDB search directory.

    Plain `.pdb` files are hard-linked (or copied); `.pdb.gz`, `.cif`, `.mmcif`,
    and gzipped mmCIF files are converted into standardized `.pdb` files so all
    downstream tools receive uniform PDB coordinates.
    """
    from structure_io import convert_structure_to_pdb, staged_pdb_name

    source_dir = Path(pdb_dir)
    destination = Path(output_dir)
    destination.mkdir(parents=True, exist_ok=True)
    staged: list[str] = []
    with Path(qc_csv).open(newline="", encoding="utf-8") as handle:
        rows = csv.DictReader(handle)
        if not rows.fieldnames or not {"fn", "pass"}.issubset(rows.fieldnames):
            raise ValueError("QC table must contain fn and pass columns")
        for row in rows:
            if str(row.get("pass", "")).strip().lower() not in {"true", "1", "yes"}:
                continue
            filename = Path(str(row["fn"])).name
            source = source_dir / filename
            if not source.is_file():
                raise FileNotFoundError(f"QC-passing structure is missing: {source}")
            out_name = staged_pdb_name(filename)
            target = destination / out_name
            if target.exists():
                raise ValueError(f"duplicate QC structure filename: {out_name}")
            convert_structure_to_pdb(source, target)
            staged.append(out_name)
    if not staged:
        raise ValueError(
            "No structures passed QC; adjust qc.min_length, qc.max_length, or "
            "qc.min_plddt and rerun"
        )
    return sorted(staged)
