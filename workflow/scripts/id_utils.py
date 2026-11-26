"""Centralized protein accession parsing and structure filename helpers.

Supports GenBank (e.g. TDZ13209.1), RefSeq (XP_001234.1), AlphaFold DB
(AF-P12345-F1-model_v4), UniProt (sp|P12345|... or bare P12345 / A0A0K2H345),
user-supplied regex patterns, and arbitrary custom stem identifiers.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

GENBANK_ACCESSION_RE = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
REFSEQ_ACCESSION_RE = re.compile(
    r"(?<![A-Za-z0-9])(?:XP|NP|WP|YP|AP|ZP)_\d+(?:\.\d+)?(?![A-Za-z0-9])"
)
AFDB_ACCESSION_RE = re.compile(
    r"(?<![A-Za-z0-9])AF-([A-Z0-9]{6,10})-F\d+(?:-model_v\d+)?(?![A-Za-z0-9])"
)
UNIPROT_PIPE_RE = re.compile(
    r"(?:sp|tr)\|([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?:-\d+)?\|"
)
UNIPROT_BARE_RE = re.compile(
    r"(?<![A-Za-z0-9])([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2})(?![A-Za-z0-9])"
)
STRUCTURE_EXT_RE = re.compile(r"\.(pdb|cif|mmcif)(\.gz)?$", flags=re.IGNORECASE)
SAFE_ID_RE = re.compile(r"^[A-Za-z0-9._:@|+=-]+$")

# Backward-compatible alias for legacy callers.
ACCESSION_RE = GENBANK_ACCESSION_RE


def strip_structure_ext(filename: str) -> str:
    """Remove supported structure file extensions (.pdb/.cif/.mmcif + optional .gz)."""
    return STRUCTURE_EXT_RE.sub("", str(filename))


def structure_ext(filename: str) -> str:
    """Return the normalized lowercase structure extension, or empty string."""
    match = STRUCTURE_EXT_RE.search(str(filename))
    return match.group(0).lower() if match else ""


def extract_accession(
    value: object,
    *,
    strain: str = "",
    custom_pattern: str | None = None,
    strip_strain_prefix: bool = False,
) -> str:
    """Extract a canonical protein accession from a filename, path, or FASTA header.

    Resolution order:
    1. Optional user-supplied regex (`custom_pattern`)
    2. GenBank versioned accession (e.g., TDZ13209.1)
    3. RefSeq accession (e.g., XP_001234.1)
    4. AlphaFold DB filename (AF-<UniProt>-F1 -> <UniProt>)
    5. UniProt header/token (sp|P12345|..., P12345, A0A0K2H345)
    6. Fallback: first whitespace token / stem without structure extension
       (optionally removing leading `<strain>_` when requested or when `strain` matches)
    """
    raw = str(value or "").strip()
    if not raw:
        return ""
    if raw.startswith(">"):
        raw = raw[1:].strip()

    # Keep basename for filesystem paths, while preserving space-separated FASTA headers.
    first_token = raw.split()[0]
    looks_like_path = ("/" in first_token or "\\" in first_token) and " " not in raw
    search_text = os.path.basename(raw) if looks_like_path else raw

    if custom_pattern:
        compiled = re.compile(custom_pattern)
        match = compiled.search(search_text)
        if match:
            return match.group(1) if match.lastindex else match.group(0)

    for regex, group_idx in (
        (GENBANK_ACCESSION_RE, 0),
        (REFSEQ_ACCESSION_RE, 0),
        (AFDB_ACCESSION_RE, 1),
        (UNIPROT_PIPE_RE, 1),
        (UNIPROT_BARE_RE, 1),
    ):
        match = regex.search(search_text)
        if match:
            return match.group(group_idx)

    stem = strip_structure_ext(Path(first_token).name)
    if strain and stem.startswith(f"{strain}_") and len(stem) > len(strain) + 1:
        return stem[len(strain) + 1 :]
    if strip_strain_prefix and "_" in stem:
        prefix, rest = stem.split("_", 1)
        if prefix and rest:
            return rest
    return stem


def normalize_structure_filename(
    filename: str,
    strain_code: str,
    *,
    custom_pattern: str | None = None,
) -> str:
    """Normalize an uploaded structure filename to `<strain>_<accession><ext>`."""
    base = os.path.basename(filename)
    ext = structure_ext(base) or ".pdb"
    acc = extract_accession(
        base,
        strain=strain_code,
        custom_pattern=custom_pattern,
        strip_strain_prefix=True,
    )
    if not acc:
        acc = strip_structure_ext(base)
    return f"{strain_code}_{acc}{ext}"
