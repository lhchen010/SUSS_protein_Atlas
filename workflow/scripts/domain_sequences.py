"""Extract domain-segment amino-acid sequences for independent D-family analysis."""

from pathlib import Path

import pandas as pd

from v3_utils import domain_sequence_records, fasta_records, write_fasta


members = pd.read_csv(snakemake.input.members)
parents = fasta_records(snakemake.input.seqs)
records, manifest = domain_sequence_records(members, parents)

fasta_path = Path(snakemake.output.fasta)
manifest_path = Path(snakemake.output.manifest)
fasta_path.parent.mkdir(parents=True, exist_ok=True)
write_fasta(records, fasta_path)
manifest.to_csv(manifest_path, index=False)
print(
    f"domain sequences: {len(records)} segments from "
    f"{manifest.acc.nunique() if len(manifest) else 0} parent proteins"
)
