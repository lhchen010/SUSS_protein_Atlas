"""Build per-family Foldseek TM and BLAST identity matrices.

The TM matrix uses NA for a pair that Foldseek did not report. A missing
measurement is not evidence for TM=0. BLAST keeps zero for an undetected HSP
because that matrix represents detected sequence relationships rather than a
global pairwise identity measurement.
"""
import os
import numpy as np, pandas as pd
from runtime_utils import symmetric_tm
from v3_utils import protein_id

fs_tsv   = snakemake.input.foldseek
bl_tsv   = snakemake.input.blastp
famfile  = snakemake.input.famfile
out_tm   = snakemake.output.tm
out_id   = snakemake.output.idm
fam      = snakemake.wildcards.fam
sym_mode = snakemake.params.sym

members = [protein_id(l.strip()) for l in open(famfile) if l.strip()]
members = list(dict.fromkeys(members))  # unique, keep order (member[0] = ref)
mset = set(members)

# --- TM from foldseek ---
fc = ["q","t","alntm","qtm","ttm","lddt","fident","aln","ql","tl","e","b"]
fs = pd.read_csv(fs_tsv, sep="\t", names=fc)
fs["qa"] = fs.q.map(protein_id); fs["ta"] = fs.t.map(protein_id)
fs = fs[fs.qa.isin(mset) & fs.ta.isin(mset)].copy()
fs["tm"] = [symmetric_tm(q, t, sym_mode) for q, t in zip(fs.qtm, fs.ttm)]
TM = pd.DataFrame(np.nan, index=members, columns=members)
for _, r in fs.iterrows():
    current = TM.loc[r.qa, r.ta]
    value = r.tm if pd.isna(current) else max(float(current), r.tm)
    TM.loc[r.qa, r.ta] = value
    TM.loc[r.ta, r.qa] = value
for m in members:
    TM.loc[m, m] = 1.0

# --- identity from blastp (pident, 0-100 -> 0-1) ---
ID = pd.DataFrame(0.0, index=members, columns=members)
if os.path.exists(bl_tsv) and os.path.getsize(bl_tsv) > 0:
    bc = ["q","t","pident","length","evalue","bitscore","qlen","slen"]
    bl = pd.read_csv(bl_tsv, sep="\t", names=bc)
    bl["qa"] = bl.q.map(protein_id); bl["ta"] = bl.t.map(protein_id)
    bl = bl[bl.qa.isin(mset) & bl.ta.isin(mset)]
    for _, r in bl.iterrows():
        v = r.pident / 100.0
        ID.loc[r.qa, r.ta] = max(ID.loc[r.qa, r.ta], v)
        ID.loc[r.ta, r.qa] = ID.loc[r.qa, r.ta]
for m in members:
    ID.loc[m, m] = 1.0

os.makedirs(os.path.dirname(out_tm), exist_ok=True)
TM.round(3).to_csv(out_tm)
ID.round(3).to_csv(out_id)
print(f"matrices {fam}: {len(members)} members, TM+ID written")
