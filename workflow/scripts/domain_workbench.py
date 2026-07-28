"""Build representative structures, superpositions, MSA, and trees for D families."""

import json
import math
import re
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from runtime_utils import resolve_executable
from v3_utils import fasta_records


def segment_pdb(path: Path, start: int, end: int) -> str:
    lines = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            residue = int(line[22:26])
        except (ValueError, IndexError):
            continue
        if start <= residue <= end:
            lines.append(line)
    return "\n".join(lines) + ("\nTER\nEND\n" if lines else "")


def safe_name(segment_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", segment_id.replace(":", "_"))


def parse_usalign(stdout: str) -> tuple[dict, dict]:
    """Parse US-align outfmt 2 summary and the structure1 -> structure2 matrix."""
    rows = []
    in_matrix = False
    summary = {}
    for line in stdout.splitlines():
        if line.startswith("#PDBchain1"):
            continue
        if "rotation matrix" in line.lower():
            in_matrix = True
            continue
        fields = line.split()
        if in_matrix and len(fields) >= 5 and fields[0] in {"0", "1", "2"}:
            rows.append([float(value) for value in fields[1:5]])
        elif (
            not line.startswith("#")
            and len(fields) >= 5
            and not summary
            and all(re.fullmatch(r"[-+.\deE]+", value) for value in fields[2:5])
        ):
            summary = {
                "tm_mobile": float(fields[2]),
                "tm_reference": float(fields[3]),
                "rmsd": float(fields[4]),
            }
    if len(rows) != 3:
        raise ValueError("US-align output did not contain a 3x4 transform matrix")
    translation = [row[0] for row in rows]
    # US-align prints x' = t + Ux. The atlas viewer applies row vectors xR + t.
    rotation = np.asarray([row[1:] for row in rows], dtype=float).T
    transform = {
        "rotation": rotation.round(10).tolist(),
        "translation": np.asarray(translation).round(10).tolist(),
    }
    return transform, summary


members = pd.read_csv(snakemake.input.members)
edges = pd.read_csv(snakemake.input.edges)
qc = pd.read_csv(snakemake.input.qc)
sequences = fasta_records(snakemake.input.seqs)
pdb_dir = Path(snakemake.input.pdb_dir)
output = Path(snakemake.output.json)
output.parent.mkdir(parents=True, exist_ok=True)

def configured_tool(value, name):
    raw = str(value or "").strip()
    return resolve_executable(raw, name, required=bool(raw))


mafft = configured_tool(snakemake.params.mafft, "MAFFT")
fasttree = configured_tool(snakemake.params.fasttree, "FastTree")
usalign = configured_tool(snakemake.params.usalign, "US-align")
threads = max(1, int(snakemake.threads))

pdb_paths = {}
for row in qc.itertuples():
    accession = str(row.acc)
    filename = getattr(row, "fn", "")
    candidates = [
        pdb_dir / str(filename),
        pdb_dir / f"{accession}.pdb",
    ]
    candidates.extend(sorted(pdb_dir.glob(f"*{accession}*.pdb")))
    pdb_paths[accession] = next(
        (candidate for candidate in candidates if candidate.is_file()), None
    )

payload = {"schema_version": 1, "families": {}}
with tempfile.TemporaryDirectory(prefix="suss-domain-workbench-") as tmp:
    tmpdir = Path(tmp)
    for family, group in members.groupby("domain_family", sort=False):
        family = str(family)
        family_edges = edges[edges.domain_family.astype(str) == family].copy()
        segment_ids = sorted(group.segment_id.astype(str))
        scores = {segment_id: 0.0 for segment_id in segment_ids}
        for edge in family_edges.itertuples():
            weight = float(edge.lddt) if pd.notna(edge.lddt) else 0.0
            scores[str(edge.source)] = scores.get(str(edge.source), 0.0) + weight
            scores[str(edge.target)] = scores.get(str(edge.target), 0.0) + weight
        hub = min(segment_ids, key=lambda segment_id: (-scores[segment_id], segment_id))

        segment_files = {}
        segment_sequences = {}
        id_map = {}
        for row in group.itertuples():
            segment_id = str(row.segment_id)
            path = pdb_paths.get(str(row.acc))
            if path:
                text = segment_pdb(path, int(row.start), int(row.end))
                if text:
                    segment_path = tmpdir / f"{family}_{safe_name(segment_id)}.pdb"
                    segment_path.write_text(text, encoding="utf-8")
                    segment_files[segment_id] = segment_path
            sequence = sequences.get(str(row.acc), "")
            if sequence:
                segment_sequences[segment_id] = sequence[int(row.start) - 1:int(row.end)]
            id_map[segment_id] = safe_name(segment_id)

        status = {
            "msa": "not_run_tool_missing" if not mafft else "not_applicable",
            "tree": "not_run_tool_missing" if not fasttree else "not_applicable",
            "superposition": "not_run_tool_missing" if not usalign else "not_applicable",
        }
        msa_records = {}
        tree = ""
        if mafft and len(segment_sequences) >= 2:
            fasta_path = tmpdir / f"{family}.fasta"
            fasta_path.write_text(
                "".join(
                    f">{id_map[segment_id]}\n{sequence}\n"
                    for segment_id, sequence in segment_sequences.items()
                ),
                encoding="utf-8",
            )
            result = subprocess.run(
                [mafft, "--auto", "--thread", str(threads), str(fasta_path)],
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            if result.returncode != 0 or not result.stdout.lstrip().startswith(">"):
                raise RuntimeError(
                    f"{family}: domain MAFFT failed ({result.returncode}): "
                    + "\n".join(result.stderr.splitlines()[-20:])
                )
            reverse_ids = {value: key for key, value in id_map.items()}
            current = None
            chunks = []
            for line in result.stdout.splitlines():
                if line.startswith(">"):
                    if current is not None:
                        msa_records[reverse_ids.get(current, current)] = "".join(chunks)
                    current = line[1:].split()[0]
                    chunks = []
                elif current is not None:
                    chunks.append(line.strip())
            if current is not None:
                msa_records[reverse_ids.get(current, current)] = "".join(chunks)
            status["msa"] = "complete"
            if fasttree:
                msa_path = tmpdir / f"{family}.aln"
                msa_path.write_text(result.stdout, encoding="utf-8")
                tree_result = subprocess.run(
                    [fasttree, "-wag", str(msa_path)],
                    capture_output=True,
                    text=True,
                    timeout=900,
                    check=False,
                )
                if tree_result.returncode != 0 or ";" not in tree_result.stdout:
                    raise RuntimeError(
                        f"{family}: domain FastTree failed ({tree_result.returncode})"
                    )
                tree = tree_result.stdout.strip()
                status["tree"] = "complete"

        transforms = {}
        fit_stats = {}
        identity = {
            "rotation": np.eye(3).tolist(),
            "translation": [0.0, 0.0, 0.0],
        }
        if hub in segment_files:
            transforms[hub] = identity
            fit_stats[hub] = {
                "reference": hub,
                "method": "reference",
                "tm_mobile": 1.0,
                "tm_reference": 1.0,
                "rmsd": 0.0,
            }
            if usalign:
                status["superposition"] = "complete"
                for segment_id, mobile_path in segment_files.items():
                    if segment_id == hub:
                        continue
                    result = subprocess.run(
                        [
                            usalign,
                            str(mobile_path),
                            str(segment_files[hub]),
                            "-outfmt",
                            "2",
                            "-m",
                            "-",
                        ],
                        capture_output=True,
                        text=True,
                        timeout=180,
                        check=False,
                    )
                    if result.returncode != 0:
                        raise RuntimeError(
                            f"{family}/{segment_id}: US-align failed "
                            f"({result.returncode}): {result.stderr[-500:]}"
                        )
                    transform, stats = parse_usalign(result.stdout)
                    transforms[segment_id] = transform
                    fit_stats[segment_id] = {
                        "reference": hub,
                        "method": "US-align domain segment",
                        **stats,
                    }

        lddt_values = pd.to_numeric(
            family_edges.get("lddt", pd.Series(dtype=float)), errors="coerce"
        ).dropna()
        alntm_values = pd.to_numeric(
            family_edges.get("alntmscore", pd.Series(dtype=float)),
            errors="coerce",
        ).dropna()
        payload["families"][family] = {
            "hub": hub,
            "members": segment_ids,
            "transforms": transforms,
            "fit_stats": fit_stats,
            "sequence_msa": msa_records,
            "sequence_newick": tree,
            "status": status,
            "mean_lddt": (
                round(float(lddt_values.mean()), 4) if len(lddt_values) else None
            ),
            "mean_alntm": (
                round(float(alntm_values.mean()), 4) if len(alntm_values) else None
            ),
        }

output.write_text(
    json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n",
    encoding="utf-8",
)
print(
    f"domain workbench: {len(payload['families'])} families -> {output}"
)
