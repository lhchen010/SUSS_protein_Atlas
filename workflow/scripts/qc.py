"""Rule 1 — QC. Parse mature AF2 PDBs, compute per-structure pLDDT/length, filter.
Input dir is prepared by the upload portal (no Drive download here — see io_contract R2).
Filename convention: <strain>_<accession>.pdb ; accession = first token before '_' after strain.
"""
import os, glob, re
import numpy as np
import pandas as pd

pdb_dir   = snakemake.input.pdb_dir
out_csv   = snakemake.output[0]
min_plddt = snakemake.params.min_plddt
max_len   = snakemake.params.max_length
min_len   = snakemake.params.min_length

def parse_acc(fn):
    match = re.search(r"[A-Z]{2,3}\d{4,}\.\d+", os.path.basename(fn))
    if match:
        return match.group(0)
    # Fallback for custom IDs: <strain>_<accession>.pdb or bare accession.pdb.
    base = os.path.basename(fn)[:-4] if fn.endswith(".pdb") else os.path.basename(fn)
    parts = base.split("_", 1)
    return parts[1] if len(parts) == 2 else parts[0]

rows = []
for fn in sorted(glob.glob(os.path.join(pdb_dir, "*.pdb"))):
    acc = parse_acc(fn)
    ca_bf = []
    with open(fn, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try: ca_bf.append(float(line[60:66]))
                except: pass
    if ca_bf:
        bf = np.array(ca_bf)
        rows.append(dict(acc=acc, fn=os.path.basename(fn), length=len(bf),
                         mean_plddt=round(float(bf.mean()), 2),
                         frac_conf=round(float((bf >= 70).mean()), 3)))

qc = pd.DataFrame(rows)
qc["pass"] = (qc.mean_plddt >= min_plddt) & (qc.length >= min_len) & (qc.length <= max_len)
qc = qc.sort_values("mean_plddt").reset_index(drop=True)
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
qc.to_csv(out_csv, index=False)
print(f"parsed {len(qc)} structures; pass QC: {int(qc['pass'].sum())} "
      f"(min_plddt={min_plddt}, len {min_len}-{max_len})")
