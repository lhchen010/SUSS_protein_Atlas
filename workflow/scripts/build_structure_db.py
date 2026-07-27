"""Create the searchable atlas Foldseek database and protein-to-family index."""

import json
import os
import shutil
import subprocess
from pathlib import Path

import pandas as pd

from runtime_utils import resolve_executable


foldseek = resolve_executable(snakemake.params.foldseek, "Foldseek structure database")
database = str(Path(snakemake.params.database))
database_path = Path(database)
database_path.parent.mkdir(parents=True, exist_ok=True)

for candidate in database_path.parent.glob(database_path.name + "*"):
    if candidate.is_file() or candidate.is_symlink():
        candidate.unlink()
    elif candidate.is_dir():
        shutil.rmtree(candidate)

result = subprocess.run(
    [foldseek, "createdb", str(snakemake.input.pdb_dir), database],
    capture_output=True,
    text=True,
    timeout=3600,
    check=False,
)
if result.returncode != 0 or not Path(database + ".dbtype").exists():
    raise RuntimeError(
        "Foldseek createdb failed: "
        + "\n".join((result.stderr + result.stdout).splitlines()[-30:])
    )

members = pd.read_csv(snakemake.input.members)
subgroups = (
    pd.read_csv(snakemake.input.subgroups)
    if Path(snakemake.input.subgroups).exists()
    else pd.DataFrame()
)
domains = (
    pd.read_csv(snakemake.input.domains)
    if Path(snakemake.input.domains).exists()
    else pd.DataFrame()
)
index = members[["acc", "family"]].copy()
if len(subgroups):
    index = index.merge(
        subgroups[["acc", "sequence_subgroup", "homology_status"]],
        on="acc",
        how="left",
    )
else:
    index["sequence_subgroup"] = pd.NA
    index["homology_status"] = pd.NA
if len(domains):
    grouped = domains.groupby("acc").agg(
        domain_families=("domain_family", lambda values: "|".join(sorted(set(values)))),
        domain_segments=(
            "segment_id",
            lambda values: "|".join(sorted(set(values))),
        ),
    )
    index = index.merge(grouped, left_on="acc", right_index=True, how="left")
else:
    index["domain_families"] = pd.NA
    index["domain_segments"] = pd.NA

Path(snakemake.output.index).parent.mkdir(parents=True, exist_ok=True)
index.to_csv(snakemake.output.index, index=False)
Path(snakemake.output.marker).write_text(
    json.dumps(
        {
            "status": "complete",
            "database": database,
            "n_proteins": len(index),
        },
        indent=2,
    )
    + "\n",
    encoding="utf-8",
)
print(f"structure database: {database} ({len(index)} indexed proteins)")
