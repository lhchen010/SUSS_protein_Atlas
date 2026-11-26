"""Rule 1 — QC. Parse mature AF2/AF3/Boltz PDB or mmCIF structures, compute
per-structure pLDDT/length, and filter.
"""
import os
import numpy as np
import pandas as pd

from id_utils import extract_accession
from structure_io import list_structure_files, read_ca_metrics, staged_pdb_name

pdb_dir = snakemake.input.pdb_dir
out_csv = snakemake.output[0]
min_plddt = snakemake.params.min_plddt
max_len = snakemake.params.max_length
min_len = snakemake.params.min_length
strain = snakemake.config.get("strain", {}).get("code", "")
custom_pattern = (
    str(snakemake.config.get("input", {}).get("accession_pattern", "") or "").strip() or None
)

rows = []
for path in list_structure_files(pdb_dir):
    acc = extract_accession(
        path.name,
        strain=strain,
        custom_pattern=custom_pattern,
        strip_strain_prefix=True,
    )
    _, ca_bf = read_ca_metrics(path)
    if ca_bf:
        bf = np.array(ca_bf, dtype=float)
        rows.append(
            dict(
                acc=acc,
                fn=path.name,
                staged_fn=staged_pdb_name(path.name),
                length=len(bf),
                mean_plddt=round(float(bf.mean()), 2),
                frac_conf=round(float((bf >= 70).mean()), 3),
            )
        )

qc = pd.DataFrame(rows)
qc["pass"] = (qc.mean_plddt >= min_plddt) & (qc.length >= min_len) & (qc.length <= max_len)
qc = qc.sort_values("mean_plddt").reset_index(drop=True)
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
qc.to_csv(out_csv, index=False)
print(
    f"parsed {len(qc)} structures; pass QC: {int(qc['pass'].sum())} "
    f"(min_plddt={min_plddt}, len {min_len}-{max_len})"
)
