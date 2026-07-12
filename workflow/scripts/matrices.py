"""matrices rule — per-family TM and sequence-identity square matrices, from the
whole-set foldseek TSV (symmetric TM = min(qtm,ttm)) and blastp TSV (pident). Written
as labeled CSVs {fam}_TM.csv / {fam}_ID.csv that cards.py and html_builder consume.
Diagonal = 1.0; missing pairs = 0 (TM) / 0 (identity)."""
import os, re
import numpy as np, pandas as pd

fs_tsv   = snakemake.input.foldseek
bl_tsv   = snakemake.input.blastp
famfile  = snakemake.input.famfile
out_tm   = snakemake.output.tm
out_id   = snakemake.output.idm
fam      = snakemake.wildcards.fam

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(s):
    m = accre.search(str(s)); return m.group(0) if m else str(s)

members = [acc_of(l) for l in open(famfile) if l.strip()]
members = list(dict.fromkeys(members))  # unique, keep order (member[0] = ref)
mset = set(members)

# --- TM from foldseek ---
fc = ["q","t","alntm","qtm","ttm","lddt","fident","aln","ql","tl","e","b"]
fs = pd.read_csv(fs_tsv, sep="\t", names=fc)
fs["qa"] = fs.q.map(acc_of); fs["ta"] = fs.t.map(acc_of)
fs = fs[fs.qa.isin(mset) & fs.ta.isin(mset)].copy()
fs["tm"] = fs[["qtm","ttm"]].min(axis=1)
TM = pd.DataFrame(0.0, index=members, columns=members)
for _, r in fs.iterrows():
    TM.loc[r.qa, r.ta] = max(TM.loc[r.qa, r.ta], r.tm)
    TM.loc[r.ta, r.qa] = TM.loc[r.qa, r.ta]
for m in members:
    TM.loc[m, m] = 1.0

# --- identity from blastp (pident, 0-100 -> 0-1) ---
ID = pd.DataFrame(0.0, index=members, columns=members)
if os.path.exists(bl_tsv) and os.path.getsize(bl_tsv) > 0:
    bc = ["q","t","pident","length","evalue","bitscore","qcovs"]
    bl = pd.read_csv(bl_tsv, sep="\t", names=bc)
    bl["qa"] = bl.q.map(acc_of); bl["ta"] = bl.t.map(acc_of)
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
