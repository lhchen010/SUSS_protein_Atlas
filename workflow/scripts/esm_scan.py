"""ESM-1b per-residue variant effects for configured protein targets.

The default ``family_representatives`` scope scans F-family references, protein
singletons, and each D family's independently selected structural hub.
``representatives`` retains the Full-only evidence scope. ``domain_members``
scans every parent protein
represented in a D family.
``all_proteins`` is available for users who explicitly accept its much higher
masked-marginals cost.
"""
import os, glob, re, shutil, subprocess
import json
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
from runtime_utils import resolve_executable, resolve_file

seqs_fa = snakemake.input.seqs
members_csv = snakemake.input.members
out_csv = snakemake.output[0]
esmpy   = snakemake.params.script
model   = snakemake.params.model
strategy= snakemake.params.strategy
scope = str(
    snakemake.params.scope or "family_representatives"
).strip().lower()
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

# Full-family aliases preserve the existing family-reference contract.
family_refs = {}
for ff in sorted(glob.glob(os.path.join(resdir, "*.members.txt"))):
    fam = os.path.basename(ff).split(".")[0]
    m = accre.search(open(ff).readline())
    if m: family_refs[fam] = m.group(0)
members = pd.read_csv(members_csv)
singletons = set(
    members.loc[members.family == "singleton", "acc"].astype(str)
)
try:
    domain_members = pd.read_csv(snakemake.input.domain_members)
except (pd.errors.EmptyDataError, FileNotFoundError):
    domain_members = pd.DataFrame(columns=["acc"])
domain_proteins = (
    set(domain_members["acc"].dropna().astype(str))
    if "acc" in domain_members.columns
    else set()
)
segment_to_accession = (
    dict(
        zip(
            domain_members["segment_id"].astype(str),
            domain_members["acc"].astype(str),
        )
    )
    if {"segment_id", "acc"}.issubset(domain_members.columns)
    else {}
)
try:
    domain_workbench = json.load(open(snakemake.input.domain_workbench))
except (FileNotFoundError, json.JSONDecodeError):
    domain_workbench = {"families": {}}
domain_hubs = {
    segment_to_accession.get(str(family.get("hub", "")))
    for family in domain_workbench.get("families", {}).values()
}
domain_hubs.discard(None)
if scope == "all_proteins":
    targets = sorted(set(members.acc.astype(str)))
elif scope == "domain_members":
    targets = sorted(
        set(family_refs.values()) | singletons | domain_proteins
    )
elif scope == "family_representatives":
    targets = sorted(
        set(family_refs.values()) | singletons | domain_hubs
    )
elif scope == "representatives":
    targets = sorted(set(family_refs.values()) | singletons)
else:
    raise ValueError(
        "signals.esm_scope must be representatives, family_representatives, "
        f"domain_members, or all_proteins, got {scope}"
    )

seq_by_acc = {acc_of(k): v for k, v in seqs.items()}
outdir = os.path.join(os.path.dirname(out_csv), "esm_out"); os.makedirs(outdir, exist_ok=True)
workers = max(1, int(snakemake.params.get("workers", 1)))


def scan_accession(accession):
    seq = seq_by_acc.get(accession)
    if not seq:
        return accession, None
    pref = os.path.join(outdir, f"{accession}_{accession}")
    mat = pref + "-res-in-matrix.csv"
    if not os.path.exists(mat):
        legacy = sorted(
            path for path in glob.glob(
                os.path.join(outdir, f"*_{accession}-res-in-matrix.csv")
            )
            if os.path.abspath(path) != os.path.abspath(mat)
        )
        if len(legacy) == 1:
            shutil.copy2(legacy[0], mat)
    if not os.path.exists(mat):
        subprocess.run([python, esmpy,
                        "--model-location", model, "--sequence", seq,
                        "--scoring-strategy", strategy, "--output-prefix", pref],
                       capture_output=True, text=True, timeout=1800, check=True)
        if not os.path.exists(mat):
            raise RuntimeError(
                f"ESM-Scan completed without expected output for {accession}"
            )
    if os.path.exists(mat):
        df = pd.read_csv(mat).rename(columns={df_c: df_c for df_c in []})
        df.insert(0, "family", accession)
        df.insert(1, "ref", accession)
        return accession, df
    return accession, None


frames_by_acc = {}
with ThreadPoolExecutor(max_workers=workers) as executor:
    for accession, frame in executor.map(scan_accession, targets):
        if frame is not None:
            frames_by_acc[accession] = frame
frames = list(frames_by_acc.values())
for family, reference in sorted(family_refs.items()):
    if reference in frames_by_acc:
        alias = frames_by_acc[reference].copy()
        alias["family"] = family
        frames.append(alias)
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
if frames:
    pd.concat(frames, ignore_index=True).to_csv(out_csv, index=False)
else:
    open(out_csv, "w").write("family,ref\n")
n_singletons = len(singletons)
print(
    f"ESM: {len(frames_by_acc)}/{len(targets)} proteins "
    f"({len(family_refs)} full-family aliases, {n_singletons} singletons, "
    f"{len(domain_hubs)} D-family hubs, {len(domain_proteins)} D-family parent "
    f"proteins, scope={scope}, workers={workers}) -> {out_csv}"
)
