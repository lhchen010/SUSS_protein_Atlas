"""Map FoldMason AA/3Di columns to hub residues and quantify structural variation."""

import os
from pathlib import Path

import numpy as np
import pandas as pd

from v3_utils import fasta_records, protein_id, shannon


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

residue_numbers = {}
coordinates = {}
column_to_index = {}
for member in members:
    numbers, coords = ca_atoms(path_by_acc[member])
    residue_numbers[member] = numbers
    coordinates[member] = coords
    mapping = {}
    sequence_index = 0
    for column, symbol in enumerate(aa_alignment[member]):
        if symbol not in {"-", "."}:
            if sequence_index < len(coords):
                mapping[column] = sequence_index
            sequence_index += 1
    column_to_index[member] = mapping

hub_numbers = residue_numbers[hub]
hub_coords = coordinates[hub]
hub_mapping = column_to_index[hub]
per_column_scores: dict[int, list[float]] = {
    column: [] for column in hub_mapping
}
thresholds = np.asarray([0.5, 1.0, 2.0, 4.0])

for member in members:
    if member == hub:
        continue
    shared_columns = sorted(set(hub_mapping) & set(column_to_index[member]))
    if len(shared_columns) < 3:
        continue
    hub_indices = np.asarray([hub_mapping[column] for column in shared_columns])
    member_indices = np.asarray(
        [column_to_index[member][column] for column in shared_columns]
    )
    hc = hub_coords[hub_indices]
    mc = coordinates[member][member_indices]
    hub_dist = np.linalg.norm(hc[:, None, :] - hc[None, :, :], axis=2)
    member_dist = np.linalg.norm(mc[:, None, :] - mc[None, :, :], axis=2)
    delta = np.abs(hub_dist - member_dist)
    for local_index, column in enumerate(shared_columns):
        neighbors = (hub_dist[local_index] > 0) & (hub_dist[local_index] < 15.0)
        if not neighbors.any():
            continue
        differences = delta[local_index, neighbors]
        score = float((differences[:, None] < thresholds[None, :]).mean())
        per_column_scores[column].append(score)

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
    scores = per_column_scores.get(column, [])
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
            "structural_lddt": float(np.mean(scores)) if scores else np.nan,
            "n_structural_comparisons": len(scores),
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
    f"{family} structural conservation: ref={hub}, columns={len(result)}, "
    f"scored={(result.structural_lddt.notna()).sum()}"
)
