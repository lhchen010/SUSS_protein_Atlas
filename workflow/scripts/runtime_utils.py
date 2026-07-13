"""Shared runtime helpers for external tools and symmetric TM scores."""

from __future__ import annotations

import os
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
