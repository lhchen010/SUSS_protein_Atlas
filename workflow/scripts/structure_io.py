"""Unified structure discovery, CA metric extraction, and mmCIF-to-PDB conversion."""

from __future__ import annotations

import gzip
import os
import shutil
from pathlib import Path
from typing import TextIO

from id_utils import STRUCTURE_EXT_RE, strip_structure_ext, structure_ext

SUPPORTED_STRUCTURE_SUFFIXES = (
    ".pdb",
    ".cif",
    ".mmcif",
    ".pdb.gz",
    ".cif.gz",
    ".mmcif.gz",
)

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "SEC": "U", "PYL": "O",
}


def is_structure_file(path: str | os.PathLike) -> bool:
    """Return True if the filename has a supported PDB or mmCIF extension."""
    return bool(STRUCTURE_EXT_RE.search(str(path)))


def is_mmcif_file(path: str | os.PathLike) -> bool:
    """Return True if the filename is an mmCIF structure (plain or gzipped)."""
    ext = structure_ext(str(path))
    return ext in {".cif", ".mmcif", ".cif.gz", ".mmcif.gz"}


def is_gzipped(path: str | os.PathLike) -> bool:
    return str(path).lower().endswith(".gz")


def staged_pdb_name(filename: str) -> str:
    """Map any supported structure filename to its canonical staged `.pdb` name."""
    base = Path(str(filename)).name
    return f"{strip_structure_ext(base)}.pdb"


def list_structure_files(directory: str | os.PathLike) -> list[Path]:
    """Return sorted supported structure files inside `directory`."""
    root = Path(directory)
    if not root.is_dir():
        return []
    files = [
        entry
        for entry in root.iterdir()
        if entry.is_file()
        and not entry.name.startswith("._")
        and is_structure_file(entry.name)
    ]
    return sorted(files, key=lambda p: p.name)


def _open_text(path: Path) -> TextIO:
    if is_gzipped(path):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _parse_mmcif_structure(path: Path):
    from Bio.PDB import MMCIFParser

    parser = MMCIFParser(QUIET=True)
    with _open_text(path) as handle:
        return parser.get_structure(strip_structure_ext(path.name), handle)


def read_ca_metrics(path: str | os.PathLike) -> tuple[int, list[float]]:
    """Return (n_ca_atoms, ca_bfactors) for a PDB or mmCIF structure."""
    file_path = Path(path)
    if not is_mmcif_file(file_path):
        ca_count = 0
        bfactors: list[float] = []
        with _open_text(file_path) as handle:
            for line in handle:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    ca_count += 1
                    try:
                        bfactors.append(float(line[60:66]))
                    except ValueError:
                        pass
        return ca_count, bfactors

    structure = _parse_mmcif_structure(file_path)
    ca_count = 0
    bfactors = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":
                    continue
                if "CA" in residue:
                    ca_count += 1
                    bfactors.append(float(residue["CA"].get_bfactor()))
        break  # First model only (matches standard AF/Boltz single-model convention)
    return ca_count, bfactors


def read_ca_sequence(path: str | os.PathLike) -> str:
    """Extract one-letter amino-acid sequence from CA atoms in PDB or mmCIF."""
    file_path = Path(path)
    if not is_mmcif_file(file_path):
        residues: dict[int, str] = {}
        with _open_text(file_path) as handle:
            for line in handle:
                if line.startswith("ATOM") and line[12:16].strip() == "CA":
                    resname = line[17:20].strip()
                    try:
                        resi = int(line[22:26])
                    except ValueError:
                        continue
                    residues[resi] = THREE_TO_ONE.get(resname, "X")
        return "".join(residues[k] for k in sorted(residues))

    structure = _parse_mmcif_structure(file_path)
    seq: list[str] = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " " or "CA" not in residue:
                    continue
                seq.append(THREE_TO_ONE.get(residue.get_resname().strip(), "X"))
        break
    return "".join(seq)


def _normalize_chain_ids_for_pdb(structure) -> None:
    """Ensure every chain ID is a single character as required by PDB format."""
    alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789abcdefghijklmnopqrstuvwxyz"
    for model in structure:
        used = {ch.id for ch in model if len(ch.id) == 1}
        pool = [c for c in alphabet if c not in used]
        for chain in list(model):
            if len(chain.id) <= 1:
                continue
            new_id = pool.pop(0) if pool else "A"
            chain.id = new_id


def convert_structure_to_pdb(source: str | os.PathLike, dest_pdb: str | os.PathLike) -> Path:
    """Stage a PDB or mmCIF structure into `dest_pdb` as a plain `.pdb` file."""
    src = Path(source)
    dst = Path(dest_pdb)
    dst.parent.mkdir(parents=True, exist_ok=True)
    ext = structure_ext(src.name)

    if ext == ".pdb":
        try:
            os.link(src, dst)
        except OSError:
            shutil.copy2(src, dst)
        return dst

    if ext == ".pdb.gz":
        with gzip.open(src, "rb") as src_fh, dst.open("wb") as dst_fh:
            shutil.copyfileobj(src_fh, dst_fh)
        return dst

    if is_mmcif_file(src):
        from Bio.PDB import PDBIO

        structure = _parse_mmcif_structure(src)
        _normalize_chain_ids_for_pdb(structure)
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(dst))
        return dst

    raise ValueError(f"Unsupported structure format: {src.name}")
