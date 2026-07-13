"""cross_edges rule — inter-family structural edges for the network view. For each pair
of families, aggregates member-vs-member symmetric TM from the whole-set
foldseek TSV: edge weight = mean TM over member pairs with TM>=foldseek_tm, tm_max = max,
n = count of such pairs. Only family pairs with >=1 qualifying pair get an edge. Emits
results/cross_family_edges.csv (from,to,tm,tm_max,n)."""
import os, re
import numpy as np, pandas as pd
from runtime_utils import symmetric_tm

fs_tsv   = snakemake.input.foldseek
members  = snakemake.input.members
out_csv  = snakemake.output[0]
TM_THR   = float(snakemake.params.tm_thr)
sym_mode = snakemake.params.sym

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(s):
    m = accre.search(str(s)); return m.group(0) if m else str(s)

mem = pd.read_csv(members)
mem = mem[mem.family != "singleton"].copy()
fam_of = dict(zip(mem.acc.astype(str), mem.family.astype(str)))

fc = ["q","t","alntm","qtm","ttm","lddt","fident","aln","ql","tl","e","b"]
fs = pd.read_csv(fs_tsv, sep="\t", names=fc)
fs["qa"] = fs.q.map(acc_of); fs["ta"] = fs.t.map(acc_of)
fs["fq"] = fs.qa.map(fam_of); fs["ft"] = fs.ta.map(fam_of)
fs = fs.dropna(subset=["fq","ft"])
fs = fs[fs.fq != fs.ft].copy()                       # inter-family only
fs["tm"] = [symmetric_tm(q, t, sym_mode) for q, t in zip(fs.qtm, fs.ttm)]
fs = fs[fs.tm >= TM_THR]
# undirected family pair key
fs["pair"] = fs.apply(lambda r: tuple(sorted([r.fq, r.ft])), axis=1)

rows = []
for (fa, fb), g in fs.groupby("pair"):
    rows.append(dict(**{"from": fa, "to": fb},
                     tm=round(float(g.tm.mean()), 3),
                     tm_max=round(float(g.tm.max()), 3),
                     n=int(len(g) // 2)))          # //2: directed pairs -> undirected
edges = pd.DataFrame(rows).sort_values("tm_max", ascending=False) if rows else \
        pd.DataFrame(columns=["from","to","tm","tm_max","n"])
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
edges.to_csv(out_csv, index=False)
print(f"cross_family_edges: {len(edges)} inter-family edges (TM>={TM_THR})")
