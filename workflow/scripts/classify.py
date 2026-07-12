"""Rule 4 — classify. Label each structural edge on the BLAST-divergence spectrum.
NOTE (io_contract R1): blast_evalue is config-driven. Lab default 1.0 (design plateau).
core_SUSS = structurally similar (TM>=thr) but NOT BLAST-detected at this evalue.
"""
import pandas as pd

blastp   = snakemake.input.blastp
edges_in = snakemake.input.edges
out_csv  = snakemake.output[0]
EVAL     = float(snakemake.params.evalue)

bcols = ["q","t","pident","length","evalue","bitscore","qcovs"]
bl = pd.read_csv(blastp, sep="\t", names=bcols)
import re as _re
def norm(s):
    s = str(s).replace(".pdb", "")
    m = _re.search(r"[A-Z]{2,3}\d{4,}\.\d+", s)
    if m: return m.group(0)
    parts = s.split("_")
    return parts[1] if len(parts) == 2 else parts[0]
bl["q"] = bl["q"].map(norm); bl["t"] = bl["t"].map(norm)
bl = bl[bl.q != bl.t].copy()
bl["pair"] = bl.apply(lambda r: tuple(sorted((r.q, r.t))), axis=1)
blp = bl.sort_values("evalue").groupby("pair", as_index=False).first()

edges = pd.read_csv(edges_in)
edges["pair"] = edges.apply(lambda r: tuple(sorted((r.q, r.t))), axis=1)
m = edges.merge(blp[["pair","pident","evalue","bitscore","qcovs"]], on="pair", how="left")
m["blast_detected"] = m.evalue.notna() & (m.evalue <= EVAL)
m["blast_any"] = m.evalue.notna()

def cls(r):
    if not r.blast_detected: return "core_SUSS"
    if r.pident < 30: return "diverged_paralog"
    if r.pident < 50: return "moderate_paralog"
    return "recent_duplicate"
m["class"] = m.apply(cls, axis=1)
m.drop(columns=["pair"]).to_csv(out_csv, index=False)
vc = m["class"].value_counts().to_dict()
print(f"classified {len(m)} edges (blast_evalue={EVAL}): {vc}")
