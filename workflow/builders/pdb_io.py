"""PDB coordinate parsing, B-factor mapping, Kabsch superposition, and ZIP bundling."""

import base64
import glob
import io
import math
import os
import re
import sys
import zipfile

import numpy as np
import pandas as pd

_SCRIPTS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "scripts"))
if _SCRIPTS_DIR not in sys.path:
    sys.path.insert(0, _SCRIPTS_DIR)
from id_utils import extract_accession

_AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLU": "E", "GLN": "Q", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def _seq_from_pdb(pdbtext):
    """One-letter sequence from CA atoms of a PDB string (fallback when no fasta)."""
    seq = {}
    for ln in pdbtext.split("\n"):
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            try: resi = int(ln[22:26])
            except ValueError: continue
            seq[resi] = _AA3TO1.get(ln[17:20].strip(), "X")
    return "".join(seq[k] for k in sorted(seq))


def _bfactor_range(pdbtext):
    """Return the finite B-factor range, used as the singleton pLDDT scale."""
    values = []
    for line in pdbtext.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                value = float(line[60:66])
            except ValueError:
                continue
            if math.isfinite(value):
                values.append(value)
    return (min(values), max(values)) if values else (0.0, 100.0)


def _pdb_with_bfactors(pdbtext, values, default=0.0):
    """Return PDB text with residue-indexed values written to the B-factor column."""
    output = []
    for line in pdbtext.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                residue = int(line[22:26])
            except ValueError:
                output.append(line)
                continue
            output.append(f"{line[:60]}{float(values.get(residue, default)):6.2f}{line[66:]}")
        else:
            output.append(line)
    return "\n".join(output)


def _ca_coordinates(pdbtext):
    coords = []
    for line in pdbtext.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except (ValueError, IndexError):
                continue
    return np.asarray(coords, dtype=float)


def _aligned_ca_pairs(ref_pdb, mobile_pdb, ref_aln=None, mobile_aln=None):
    """Return corresponding CA coordinates, preferably from the FoldMason MSA."""
    ref_ca = _ca_coordinates(ref_pdb)
    mob_ca = _ca_coordinates(mobile_pdb)
    if ref_aln and mobile_aln and len(ref_aln) == len(mobile_aln):
        pairs = []
        ri = mi = 0
        for ra, ma in zip(ref_aln, mobile_aln):
            rpos = ri if ra not in "-." else None
            mpos = mi if ma not in "-." else None
            if rpos is not None:
                ri += 1
            if mpos is not None:
                mi += 1
            if rpos is not None and mpos is not None and rpos < len(ref_ca) and mpos < len(mob_ca):
                pairs.append((rpos, mpos))
        if len(pairs) >= 3:
            return (ref_ca[[p[0] for p in pairs]], mob_ca[[p[1] for p in pairs]], "foldmason")
    n = min(len(ref_ca), len(mob_ca))
    if n < 3:
        raise ValueError("at least three paired CA atoms are required for superposition")
    return ref_ca[:n], mob_ca[:n], "ca_order"


def _superpose_pdb(mobile_pdb, ref_pdb, mobile_aln=None, ref_aln=None):
    """Rigidly align a PDB to a reference and return transformed text plus fit metadata."""
    ref_xyz, mob_xyz, method = _aligned_ca_pairs(ref_pdb, mobile_pdb, ref_aln, mobile_aln)
    def fit(mask):
        ref_fit = ref_xyz[mask]
        mob_fit = mob_xyz[mask]
        ref_center = ref_fit.mean(axis=0)
        mob_center = mob_fit.mean(axis=0)
        u, _, vt = np.linalg.svd((mob_fit - mob_center).T @ (ref_fit - ref_center))
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        return rotation, ref_center - mob_center @ rotation

    # FoldMason columns provide correspondence, while iterative rejection prevents long
    # flexible loops from pulling the conserved core away from the hub. Four angstroms
    # is deliberately permissive; if an initial fit has too few inliers, retain the
    # closest half and iterate rather than failing a divergent but valid family.
    mask = np.ones(len(ref_xyz), dtype=bool)
    for _ in range(8):
        rot, tran = fit(mask)
        distances = np.linalg.norm(mob_xyz @ rot + tran - ref_xyz, axis=1)
        new_mask = distances <= 4.0
        if new_mask.sum() < 3:
            keep = max(3, int(math.ceil(len(ref_xyz) * 0.5)))
            new_mask = np.zeros(len(ref_xyz), dtype=bool)
            new_mask[np.argsort(distances)[:keep]] = True
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask
    rot, tran = fit(mask)
    fitted = mob_xyz @ rot + tran
    squared = np.sum((fitted - ref_xyz) ** 2, axis=1)
    rmsd = float(np.sqrt(np.mean(squared[mask])))
    rmsd_all = float(np.sqrt(np.mean(squared)))
    output = []
    for line in mobile_pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                xyz = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                output.append(line)
                continue
            x, y, z = xyz @ rot + tran
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        output.append(line)
    return "\n".join(output) + "\n", {
        "method": method, "n_ca": int(mask.sum()), "n_ca_total": len(ref_xyz),
        "rmsd": round(rmsd, 4), "rmsd_all": round(rmsd_all, 4),
        "rotation": rot.round(10).tolist(), "translation": tran.round(10).tolist(),
    }


def _structures_zip_b64(fam, structures):
    """Create a family ZIP containing one independently usable PDB per member."""
    if not structures:
        return ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for member, pdbtext in structures.items():
            archive.writestr(f"{fam}_structures/{member}.pdb", pdbtext.rstrip() + "\n")
        manifest = "family\tmember\tfile\n" + "".join(
            f"{fam}\t{member}\t{member}.pdb\n" for member in structures)
        archive.writestr(f"{fam}_structures/manifest.tsv", manifest)
    return base64.b64encode(buf.getvalue()).decode()


def _structure_path(accession, input_dir="", strain="", family_dir=None, fallback_dir=None):
    """Resolve one structure by accession without assuming a filename prefix."""
    candidates = []
    if family_dir:
        candidates.append(os.path.join(family_dir, f"{accession}.pdb"))
    for root in (input_dir, fallback_dir):
        if not root:
            continue
        if strain:
            candidates.append(os.path.join(root, f"{strain}_{accession}.pdb"))
        candidates.append(os.path.join(root, f"{accession}.pdb"))

    for candidate in dict.fromkeys(candidates):
        if os.path.isfile(candidate):
            return candidate

    prefixed_matches = []
    for root in dict.fromkeys(root for root in (input_dir, fallback_dir) if root):
        prefixed_matches.extend(
            glob.glob(os.path.join(root, f"*_{glob.escape(accession)}.pdb"))
        )
    prefixed_matches = sorted(set(prefixed_matches))
    if len(prefixed_matches) == 1:
        return prefixed_matches[0]
    if len(prefixed_matches) > 1:
        raise RuntimeError(
            f"Ambiguous structures for {accession}: "
            + ", ".join(prefixed_matches)
        )
    # Fallback: match by extracted accession (e.g. AF-<UniProt>-F1-model_v4.pdb)
    extracted_matches = []
    for root in dict.fromkeys(root for root in (input_dir, fallback_dir) if root and os.path.isdir(root)):
        for entry in sorted(os.listdir(root)):
            if not entry.lower().endswith(".pdb"):
                continue
            if extract_accession(entry, strain=strain, strip_strain_prefix=True) == accession:
                extracted_matches.append(os.path.join(root, entry))
    if len(extracted_matches) == 1:
        return extracted_matches[0]
    if len(extracted_matches) > 1:
        raise RuntimeError(
            f"Ambiguous structures for {accession}: "
            + ", ".join(extracted_matches)
        )
    return ""


def _domain_structures_zip_b64(fam, members, parent_structures, segments_only):
    """Create a ZIP of cropped domain members or their complete parent proteins."""
    buf = io.BytesIO()
    written = set()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        manifest_rows = []
        for member in members:
            accession = str(member["acc"])
            parent = parent_structures.get(accession, "")
            if not parent:
                continue
            if segments_only:
                name = re.sub(r"[^A-Za-z0-9_.-]+", "_", member["segment_id"])
                text = "\n".join(
                    line for line in parent.splitlines()
                    if line.startswith(("ATOM", "HETATM"))
                    and line[22:26].strip().lstrip("-").isdigit()
                    and int(line[22:26]) >= int(member["start"])
                    and int(line[22:26]) <= int(member["end"])
                )
                filename = f"{name}.pdb"
                archive.writestr(
                    f"{fam}_domain_structures/{filename}",
                    text.rstrip() + "\nTER\nEND\n",
                )
                manifest_rows.append(
                    {
                        "segment_id": member["segment_id"],
                        "accession": accession,
                        "start": member["start"],
                        "end": member["end"],
                        "file": filename,
                    }
                )
            elif accession not in written:
                written.add(accession)
                filename = f"{accession}.pdb"
                archive.writestr(
                    f"{fam}_parent_structures/{filename}",
                    parent.rstrip() + "\n",
                )
                manifest_rows.append(
                    {
                        "segment_id": "",
                        "accession": accession,
                        "start": "",
                        "end": "",
                        "file": filename,
                    }
                )
        manifest = pd.DataFrame(
            manifest_rows,
            columns=["segment_id", "accession", "start", "end", "file"],
        ).to_csv(index=False)
        root = (
            f"{fam}_domain_structures"
            if segments_only
            else f"{fam}_parent_structures"
        )
        archive.writestr(f"{root}/manifest.csv", manifest)
    return base64.b64encode(buf.getvalue()).decode()

