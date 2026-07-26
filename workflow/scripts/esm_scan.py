"""ESM-1b per-residue variant effects for family references and singletons.

Singleton records use their accession in the ``family`` column so the atlas can
attach the result without pretending that all singletons form one family.
"""
import os, glob, re, subprocess
import pandas as pd
from runtime_utils import resolve_executable, resolve_file

seqs_fa = snakemake.input.seqs
members_csv = snakemake.input.members
out_csv = snakemake.output[0]
esmpy   = snakemake.params.script
model   = snakemake.params.model
strategy= snakemake.params.strategy
resdir  = os.path.join(os.path.dirname(out_csv), "families")
if not snakemake.params.enabled:
    raise RuntimeError("ESM rule was scheduled while steps.esm is disabled")
esmpy = resolve_file(esmpy, "ESM-Scan")
python = resolve_executable(snakemake.params.get("python", "python"), "ESM Python")

# read seqs
seqs = {}
name = None
for line in open(seqs_fa, encoding="utf-8", errors="replace"):
    line = line.strip()
    if line.startswith(">"): name = line[1:].split()[0]; seqs[name] = ""
    elif name: seqs[name] += line
accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(s):
    m = accre.search(s); return m.group(0) if m else s

# reference per family
refs = {}
for ff in sorted(glob.glob(os.path.join(resdir, "*.members.txt"))):
    fam = os.path.basename(ff).split(".")[0]
    m = accre.search(open(ff).readline())
    if m: refs[fam] = m.group(0)
members = pd.read_csv(members_csv)
for acc in sorted(members.loc[members.family == "singleton", "acc"].astype(str)):
    refs[acc] = acc

seq_by_acc = {acc_of(k): v for k, v in seqs.items()}
outdir = os.path.join(os.path.dirname(out_csv), "esm_out"); os.makedirs(outdir, exist_ok=True)
frames = []
for fam, ref in sorted(refs.items()):
    seq = seq_by_acc.get(ref)
    if not seq: continue
    pref = os.path.join(outdir, f"{fam}_{ref}")
    mat = pref + "-res-in-matrix.csv"
    if not os.path.exists(mat):
        subprocess.run([python, esmpy,
                        "--model-location", model, "--sequence", seq,
                        "--scoring-strategy", strategy, "--output-prefix", pref],
                       capture_output=True, text=True, timeout=1800, check=True)
        if not os.path.exists(mat):
            raise RuntimeError(f"ESM-Scan completed without expected output for {fam}")
    if os.path.exists(mat):
        df = pd.read_csv(mat).rename(columns={df_c: df_c for df_c in []})
        df.insert(0, "family", fam); df.insert(1, "ref", ref)
        frames.append(df)
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
if frames:
    pd.concat(frames, ignore_index=True).to_csv(out_csv, index=False)
else:
    open(out_csv, "w").write("family,ref\n")
n_singletons = int((members.family == "singleton").sum())
print(f"ESM: {len(frames)} references ({n_singletons} singleton targets) -> {out_csv}")
