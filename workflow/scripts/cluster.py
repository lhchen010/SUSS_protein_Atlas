"""Rule 3 — cluster (CHECKPOINT). Foldseek all-vs-all -> symmetric TM -> Leiden.
Writes families.csv, members.csv, edges.csv, and one <fam>.members.txt per family
into famdir (this is what the checkpoint-aware DAG globs to expand per-family rules).
Config: foldseek_tm, tm_symmetric, leiden_resolution, leiden_seed, min_family_size.
"""
import os
import numpy as np
import pandas as pd
import igraph as ig
import leidenalg

tsv      = snakemake.input.tsv
qc_csv   = snakemake.input.qc
out_fam  = snakemake.output.families
out_mem  = snakemake.output.members
out_edge = snakemake.output.edges
famdir   = snakemake.output.famdir
TM_THR   = float(snakemake.params.tm)
SYM      = snakemake.params.sym
RES      = float(snakemake.params.res)
SEED     = int(snakemake.params.seed)
MIN_SIZE = int(snakemake.params.min_size)

cols = ["query","target","alntmscore","qtmscore","ttmscore","lddt","fident",
        "alnlen","qlen","tlen","evalue","bits"]
df = pd.read_csv(tsv, sep="\t", names=cols)
# Foldseek names may carry either a strain prefix (cor_TDZ...) or an AF2 suffix
# (TDZ..._unrelaxed_rank_...). Accession = the TDZ\d+\.\d+ token; fall back to basename.
import re as _re
def norm(s):
    s = str(s).replace(".pdb", "")
    m = _re.search(r"[A-Z]{2,3}\d{4,}\.\d+", s)   # GenBank-style accession
    if m: return m.group(0)
    parts = s.split("_")
    return parts[1] if len(parts) == 2 else parts[0]
df["q"] = df["query"].map(norm); df["t"] = df["target"].map(norm)
if SYM == "min":   df["tm"] = df[["qtmscore","ttmscore"]].min(axis=1)
elif SYM == "max": df["tm"] = df[["qtmscore","ttmscore"]].max(axis=1)
else:              df["tm"] = df[["qtmscore","ttmscore"]].mean(axis=1)
off = df[df.q != df.t]

edges = off[off.tm >= TM_THR][["q","t","tm","fident"]].copy()
edges["pair"] = edges.apply(lambda r: tuple(sorted((r.q, r.t))), axis=1)
edges = edges.groupby("pair", as_index=False).agg(tm=("tm","max"), fident=("fident","max"))
edges[["q","t"]] = pd.DataFrame(edges.pair.tolist(), index=edges.index)

qc = pd.read_csv(qc_csv)
qc = qc[qc["pass"]] if "pass" in qc.columns else qc
plddt = dict(zip(qc.acc, qc.mean_plddt)); length = dict(zip(qc.acc, qc.length))
all_nodes = sorted(qc.acc.unique())

g = ig.Graph(); g.add_vertices(all_nodes)
keep = edges[(edges.q.isin(all_nodes)) & (edges.t.isin(all_nodes))]
g.add_edges(list(zip(keep.q, keep.t))); g.es["weight"] = keep.tm.tolist()
part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition,
                                weights="weight", resolution_parameter=RES, seed=SEED)
memb = np.array(part.membership)
comm = pd.DataFrame({"acc": g.vs["name"], "community": memb, "deg": g.degree()})
sizes = comm.community.value_counts()
fam_ids = sizes[sizes >= MIN_SIZE].index.tolist()
# rank families by size -> F0, F1, ...
fam_sizes = comm[comm.community.isin(fam_ids)].community.value_counts()
fam_rank = {fid: f"F{i}" for i, fid in enumerate(fam_sizes.index)}
comm["family"] = comm.community.map(fam_rank).fillna("singleton")
comm["plddt"] = comm.acc.map(plddt); comm["length"] = comm.acc.map(length)

# families.csv (one row per family)
rows = []
for fid in fam_sizes.index:
    fam = fam_rank[fid]; accs = comm[comm.community == fid].acc.tolist()
    e_in = keep[(keep.q.isin(accs)) & (keep.t.isin(accs))]
    rows.append(dict(family=fam, community=int(fid), n_members=len(accs),
                     n_edges=len(e_in),
                     mean_TM=round(e_in.tm.mean(),3) if len(e_in) else np.nan,
                     mean_identity=round(e_in.fident.mean(),3) if len(e_in) else np.nan,
                     max_identity=round(e_in.fident.max(),3) if len(e_in) else np.nan,
                     mean_pLDDT=round(comm[comm.community==fid].plddt.mean(),1),
                     mean_len=int(comm[comm.community==fid].length.mean())))
fam = pd.DataFrame(rows).sort_values("n_members", ascending=False).reset_index(drop=True)

os.makedirs(os.path.dirname(out_fam), exist_ok=True)
os.makedirs(famdir, exist_ok=True)
fam.to_csv(out_fam, index=False)
comm[["acc","family","community","deg","plddt","length"]].to_csv(out_mem, index=False)
keep[["q","t","tm","fident"]].to_csv(out_edge, index=False)

# per-family member-file (paths to PDBs) for downstream per-family rules.
# HUB = member with highest mean within-family TM ("most like everyone"); written FIRST
# so FoldMason MSA, rate4site -a, pocket ref, and hub marking all use the SAME protein.
pdb_dir = snakemake.config["input"]["pdb_dir"]
acc2fn = dict(zip(qc.acc, qc.fn)) if "fn" in qc.columns else {}
# symmetric TM per undirected pair for hub scoring
tm_pair = keep.copy()
for fid in fam_sizes.index:
    famname = fam_rank[fid]; accs = comm[comm.community == fid].acc.tolist()
    aset = set(accs)
    sub = tm_pair[tm_pair.q.isin(aset) & tm_pair.t.isin(aset)]
    if len(sub) and len(accs) > 2:
        meantm = {}
        for a in accs:
            v = sub[(sub.q == a) | (sub.t == a)].tm
            meantm[a] = float(v.mean()) if len(v) else 0.0
        hub = max(accs, key=lambda a: meantm.get(a, 0.0))
    else:
        hub = accs[0]
    ordered = [hub] + [a for a in accs if a != hub]
    with open(os.path.join(famdir, f"{famname}.members.txt"), "w") as fh:
        for a in ordered:
            fn = acc2fn.get(a, f"*_{a}.pdb")
            fh.write(os.path.join(pdb_dir, fn) + "\n")
print(f"clustered: {len(fam)} families (min_size={MIN_SIZE}), "
      f"{int(fam.n_members.sum())} proteins in families, "
      f"{(comm.family=='singleton').sum()} singletons; edges={len(keep)} (TM>={TM_THR})")
