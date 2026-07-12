"""master rule — assemble the family-level master table: per family n, TM stats,
identity, SUSS%, conservation-SASA r, effector/novel counts, consensus annotation.
Joins families.csv + classification + signature (per family) + annotation. No dN/dS
(removed per project ruling). (recovered from master.py, de-hardcoded.)
"""
import os, glob, re
import numpy as np, pandas as pd

fam_csv  = snakemake.input.families
anno_csv = snakemake.input.annotation
cls_csv  = getattr(snakemake.input, "classification", None)   # gated off → absent
out_csv  = snakemake.output[0]
resdir   = os.path.dirname(out_csv)
sigdir   = os.path.join(resdir, "families")

fam = pd.read_csv(fam_csv)
anno = pd.read_csv(anno_csv)
cls = pd.read_csv(cls_csv) if cls_csv and os.path.exists(cls_csv) else pd.DataFrame(columns=["q","t","class"])

# SUSS% per family = fraction of within-family structural edges that are core_SUSS
def suss_pct(f):
    accs = set(anno[anno.family == f].acc)
    sub = cls[(cls.q.isin(accs)) & (cls.t.isin(accs))] if {"q","t"}.issubset(cls.columns) else cls.iloc[0:0]
    return round(100 * (sub["class"] == "core_SUSS").mean(), 1) if len(sub) else np.nan

# conservation-SASA r per family from signature csv
def cons_r(f):
    p = os.path.join(sigdir, f, f"{f}_signature.csv")
    if not os.path.exists(p): return np.nan
    m = pd.read_csv(p); v = m.dropna(subset=["rel_sasa","conservation"])
    return round(float(np.corrcoef(v.conservation, v.rel_sasa)[0,1]), 3) if len(v) > 2 else np.nan

rows = []
for _, r in fam.iterrows():
    f = r.family; g = anno[anno.family == f]
    rows.append(dict(
        family=f, n_members=int(r.n_members),
        mean_TM=r.get("mean_TM"), mean_identity=r.get("mean_identity"), max_identity=r.get("max_identity"),
        mean_pLDDT=r.get("mean_pLDDT"), mean_len=r.get("mean_len"),
        suss_pct=suss_pct(f), cons_sasa_r=cons_r(f),
        n_effector=int(g.is_effector.sum()) if len(g) else 0,
        pct_effector=round(100*g.is_effector.mean(),1) if len(g) else np.nan,
        n_novel=int(g.novel.sum()) if len(g) else 0,
        pct_novel=round(100*g.novel.mean(),1) if len(g) else np.nan,
        n_small=("n small" if int(r.n_members) <= 3 else ""),
    ))
master = pd.DataFrame(rows).sort_values("n_members", ascending=False)
master.to_csv(out_csv, index=False)
print(f"master table: {len(master)} families -> {out_csv}")

