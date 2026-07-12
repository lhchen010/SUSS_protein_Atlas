"""esm rule — ESM-1b per-residue variant-effect (masked-marginals) for each family
reference sequence. Runs esmscan.py on 4070 (boltzgen env, GPU); weights auto-download
first run. Emits esm_all.csv = long-form (fam, resi, wt, per-AA LLR). (recovered from
esm_driver.py.)
"""
import os, glob, re, subprocess
import pandas as pd

seqs_fa = snakemake.input.seqs
out_csv = snakemake.output[0]
esmpy   = snakemake.params.script
model   = snakemake.params.model
strategy= snakemake.params.strategy
resdir  = os.path.join(os.path.dirname(out_csv), "families")

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

seq_by_acc = {acc_of(k): v for k, v in seqs.items()}
outdir = os.path.join(os.path.dirname(out_csv), "esm_out"); os.makedirs(outdir, exist_ok=True)
frames = []
for fam, ref in sorted(refs.items(), key=lambda x: int(x[0][1:]) if x[0][1:].isdigit() else 0):
    seq = seq_by_acc.get(ref)
    if not seq: continue
    pref = os.path.join(outdir, f"{fam}_{ref}")
    mat = pref + "-res-in-matrix.csv"
    if not os.path.exists(mat):
        subprocess.run([snakemake.params.get("python", "python"), esmpy,
                        "--model-location", model, "--sequence", seq,
                        "--scoring-strategy", strategy, "--output-prefix", pref],
                       capture_output=True, text=True, timeout=1800)
    if os.path.exists(mat):
        df = pd.read_csv(mat).rename(columns={df_c: df_c for df_c in []})
        df.insert(0, "family", fam); df.insert(1, "ref", ref)
        frames.append(df)
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
if frames:
    pd.concat(frames, ignore_index=True).to_csv(out_csv, index=False)
else:
    open(out_csv, "w").write("family,ref\n")
print(f"ESM: {len(frames)} family references scored -> {out_csv}")

