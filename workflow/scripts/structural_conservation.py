"""Map official FoldMason per-column LDDT scores to the hub structure."""

import json
import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from runtime_utils import resolve_executable
from v3_utils import (
    fasta_records,
    foldmason_column_scores,
    protein_id,
    shannon,
    write_fasta,
)


def ca_atoms(path: str | Path) -> tuple[list[int], np.ndarray]:
    residue_numbers = []
    coordinates = []
    for line in Path(path).read_text(
        encoding="utf-8", errors="replace"
    ).splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                residue_numbers.append(int(line[22:26]))
                coordinates.append(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                )
            except ValueError:
                continue
    return residue_numbers, np.asarray(coordinates, dtype=float)


family = snakemake.wildcards.fam
paths = [
    line.strip()
    for line in Path(snakemake.input.famfile).read_text(encoding="utf-8").splitlines()
    if line.strip()
]
hub = protein_id(paths[0])
aa_alignment = fasta_records(snakemake.input.aa)
di_alignment = fasta_records(snakemake.input.di)
path_by_acc = {protein_id(path): path for path in paths}
members = [protein_id(path) for path in paths if protein_id(path) in aa_alignment]

if hub not in aa_alignment or hub not in path_by_acc:
    raise RuntimeError(f"{family}: FoldMason alignment does not contain hub {hub}")
alignment_length = len(aa_alignment[hub])
if any(len(aa_alignment[member]) != alignment_length for member in members):
    raise RuntimeError(f"{family}: FoldMason AA alignment has inconsistent lengths")

foldmason = resolve_executable(snakemake.params.foldmason, "FoldMason")
pair_threshold = float(snakemake.params.pair_threshold)
with tempfile.TemporaryDirectory(prefix=f"{family}_lddt_") as tmp:
    tmp_path = Path(tmp)
    structure_dir = tmp_path / "structures"
    structure_dir.mkdir()
    for member in members:
        shutil.copy2(path_by_acc[member], structure_dir / f"{member}.pdb")
    normalized_msa = tmp_path / "alignment.fasta"
    write_fasta({member: aa_alignment[member] for member in members}, normalized_msa)
    database = tmp_path / "structures_db"
    output_json = tmp_path / "column_lddt.json"
    subprocess.run(
        [foldmason, "createdb", str(structure_dir), str(database), "--threads", "1"],
        check=True,
        capture_output=True,
        text=True,
    )
    subprocess.run(
        [
            foldmason,
            "msa2lddtjson",
            str(database),
            str(normalized_msa),
            str(output_json),
            "--pair-threshold",
            str(pair_threshold),
            "--threads",
            str(max(1, int(snakemake.threads))),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    column_scores = foldmason_column_scores(
        json.loads(output_json.read_text(encoding="utf-8")), alignment_length
    )

hub_numbers, _ = ca_atoms(path_by_acc[hub])
rows = []
hub_residue_index = 0
for column, hub_symbol in enumerate(aa_alignment[hub]):
    if hub_symbol in {"-", "."}:
        continue
    if hub_residue_index >= len(hub_numbers):
        break
    aa_values = [aa_alignment[member][column] for member in members]
    di_values = [
        di_alignment.get(member, "")[column]
        for member in members
        if len(di_alignment.get(member, "")) > column
    ]
    present = sum(value not in {"-", "."} for value in aa_values)
    score = column_scores[column]
    rows.append(
        {
            "family": family,
            "reference": hub,
            "resi": hub_numbers[hub_residue_index],
            "aa": hub_symbol,
            "alignment_column": column + 1,
            "occupancy": present / len(members) if members else np.nan,
            "aa_entropy": shannon(aa_values),
            "three_di_entropy": shannon(di_values),
            "structural_lddt": score,
            "n_structural_members": present,
            "pair_threshold": pair_threshold,
            "scoring_method": "FoldMason msa2lddtjson",
        }
    )
    hub_residue_index += 1

result = pd.DataFrame(rows)
os.makedirs(os.path.dirname(snakemake.output.csv), exist_ok=True)
result.to_csv(snakemake.output.csv, index=False)

score_by_residue = dict(zip(result.resi, result.structural_lddt))
output_lines = []
for line in Path(path_by_acc[hub]).read_text(
    encoding="utf-8", errors="replace"
).splitlines(keepends=True):
    if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
        try:
            residue = int(line[22:26])
        except ValueError:
            output_lines.append(line)
            continue
        score = score_by_residue.get(residue)
        bfactor = 0.0 if pd.isna(score) else 100.0 * float(score)
        output_lines.append(f"{line[:60]}{bfactor:6.2f}{line[66:]}")
    else:
        output_lines.append(line)
Path(snakemake.output.pdb).write_text("".join(output_lines), encoding="utf-8")
print(
    f"{family} FoldMason structural conservation: ref={hub}, columns={len(result)}, "
    f"scored={(result.structural_lddt.notna()).sum()}, "
    f"pair_threshold={pair_threshold}"
)
