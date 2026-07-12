"""annotate rule — 4-tool annotation for ALL proteins (incl. singletons):
Foldseek vs pdb100 + afdb_swissprot (known fold + real protein NAME), InterProScan
(Pfam domains + InterPro entries, multi-domain flag), EffectorP (effector call),
DeepTMHMM (TM region count). No SignalP (SP already removed pre-AF2 per project rule).
Runs the four tools on 4070 (this driver assumes their raw outputs are staged under
results/anno/), then merges to member_annotation.csv + cluster_annotation.csv.
novel = no known fold AND no domain. (recovered from annotation.py, de-hardcoded.)
"""
import os, re
from collections import defaultdict
import numpy as np, pandas as pd

pdb_dir  = snakemake.input.pdb_dir
seqs_fa  = snakemake.input.seqs
mem_csv  = snakemake.input.members
out_mem  = snakemake.output.member
out_clu  = snakemake.output.cluster
anno_dir = os.path.join(os.path.dirname(out_mem), "anno")   # where tool outputs are staged
os.makedirs(anno_dir, exist_ok=True)

# ---- run the four annotation tools if their outputs are not already staged ----
# (previously this driver only READ results/anno/*.tsv and nothing produced them, so a
#  fresh run had empty annotation. Now it runs them, mirroring the commands that produced
#  the 517 annotation. EffectorP/DeepTMHMM run under an explicit interpreter from config
#  (DeepTMHMM needs torch -> use the GPU env python, not the suss env).)
import subprocess, shutil
P = snakemake.params
def _sh(cmd, log, timeout=5400):
    with open(log, "a") as lf:
        lf.write("\n$ " + cmd + "\n")
        return subprocess.run(cmd, shell=True, stdout=lf, stderr=subprocess.STDOUT, timeout=timeout).returncode

alog = os.path.join(anno_dir, "annotate_tools.log")
# sanitize sequences once (strip '*' stop chars that break InterProScan/EffectorP)
clean_fa = os.path.join(anno_dir, "clean.fasta")
if not os.path.exists(clean_fa):
    _sh(f"sed '/^>/!s/[*]//g' {seqs_fa} > {clean_fa}", alog, timeout=300)

# InterProScan -> interpro.tsv (Pfam + CDD + Gene3D)
ipr_tsv = os.path.join(anno_dir, "interpro.tsv")
if not os.path.exists(ipr_tsv) and getattr(P, "ips", "") and os.path.exists(str(P.ips)):
    try: _sh(f"{P.ips} -i {clean_fa} -f TSV -o {ipr_tsv} -appl Pfam,CDD,Gene3D --cpu 24 --goterms", alog)
    except Exception as e: open(alog,"a").write(f"\nInterProScan failed: {e}\n")

# Foldseek easy-search vs pdb100 and afdb_swissprot -> foldseek_pdb100.tsv / foldseek_afdbsp.tsv
fs = getattr(P, "fs", ""); fsdb = getattr(P, "fsdb", "")
fmt = "query,target,evalue,alntmscore,fident,alnlen,prob,qstart,qend,tstart,tend"
if fs and fsdb and os.path.exists(str(fs)):
    fsq = os.path.join(anno_dir, "fsq_db")
    if not os.path.exists(fsq):
        _sh(f"{fs} createdb {pdb_dir} {fsq} 2>&1", alog, timeout=1800)
    for db, out in [("pdb100", "foldseek_pdb100.tsv"), ("afdb_swissprot", "foldseek_afdbsp.tsv")]:
        op = os.path.join(anno_dir, out); dbp = os.path.join(str(fsdb), db)
        db_present = os.path.exists(dbp) or os.path.exists(dbp + ".dbtype")
        if not os.path.exists(op) and db_present:
            try: _sh(f"{fs} easy-search {fsq} {dbp} {op} {anno_dir}/tmp_{db} "
                     f'--format-output "{fmt}" -e 0.01 --max-seqs 10 --threads 24 2>&1', alog)
            except Exception as e: open(alog,"a").write(f"\nFoldseek {db} failed: {e}\n")

# tool interpreters (config-driven): DeepTMHMM needs torch+h5py -> base miniconda python;
# EffectorP runs fine under the suss env python3. Fall back to PATH python3 if unset.
effp_py  = str(getattr(P, "effp_python", "") or "") or shutil.which("python3") or "python3"
tmhmm_py = str(getattr(P, "tmhmm_python", "") or "") or shutil.which("python3") or "python3"

# EffectorP -> effectorp.txt
effp_out = os.path.join(anno_dir, "effectorp.txt")
effp = getattr(P, "effp", "")
if not os.path.exists(effp_out) and effp and os.path.exists(str(effp)):
    try: _sh(f"{effp_py} {effp} -i {clean_fa} > {effp_out} 2>&1", alog, timeout=2400)
    except Exception as e: open(alog,"a").write(f"\nEffectorP failed: {e}\n")

# DeepTMHMM -> anno/tmhmm/TMRs.gff3 (run from its install dir; tolerate final matplotlib crash)
tmhmm = getattr(P, "tmhmm", "")
tm_gff = os.path.join(anno_dir, "tmhmm_TMRs.gff3")
if not os.path.exists(tm_gff) and tmhmm and os.path.exists(str(tmhmm)):
    tmdir = os.path.dirname(str(tmhmm)); pybin = tmhmm_py
    tmout = os.path.join(anno_dir, "tmhmm")
    try:
        _sh(f'cd {tmdir} && {pybin} predict.py --fasta "{os.path.abspath(clean_fa)}" '
            f'--output-dir "{os.path.abspath(tmout)}" 2>&1 | tail -4 || true', alog, timeout=3600)
        for cand in ("TMRs.gff3", "predicted_topologies.3line"):
            src = os.path.join(tmout, cand)
            if os.path.exists(src) and cand == "TMRs.gff3": shutil.copy(src, tm_gff)
    except Exception as e: open(alog,"a").write(f"\nDeepTMHMM failed: {e}\n")

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(s):
    m = accre.search(s); return m.group(0) if m else s.split()[0]

# accession universe from seqs fasta
accs = []
for line in open(seqs_fa, encoding="utf-8", errors="replace"):
    if line.startswith(">"): accs.append(acc_of(line[1:]))
accs = sorted(set(accs))

# --- InterProScan (results/anno/interpro.tsv) ---
member_domains = defaultdict(list); ipr_entries = defaultdict(set)
ipr_path = os.path.join(anno_dir, "interpro.tsv")
if os.path.exists(ipr_path):
    ipr = pd.read_csv(ipr_path, sep="\t", header=None,
                      names=["acc","md5","len","db","sig","desc","start","end","evalue","status","date","ipr","ipr_desc","go"],
                      usecols=range(14), on_bad_lines="skip")
    ipr["acc"] = ipr["acc"].map(acc_of)
    for _, r in ipr[ipr.db == "Pfam"].iterrows():
        try: member_domains[r.acc].append({"desc": r.desc, "start": int(r.start), "end": int(r.end)})
        except: pass
    for _, r in ipr.iterrows():
        if isinstance(r.ipr_desc, str) and r.ipr_desc not in ("-", ""): ipr_entries[r.acc].add(r.ipr_desc)

# --- Foldseek best hit vs pdb100 / afdb_swissprot ---
fscols = ["query","target","evalue","tm","fident","alnlen","prob","qs","qe","ts","te"]
def best_fs(name):
    p = os.path.join(anno_dir, name)
    if not os.path.exists(p): return {}
    df = pd.read_csv(p, sep="\t", header=None, names=fscols)
    df["query"] = df["query"].map(acc_of)
    return df.loc[df.groupby("query").evalue.idxmin()].set_index("query").to_dict("index")
fs_pdb = best_fs("foldseek_pdb100.tsv"); fs_afdb = best_fs("foldseek_afdbsp.tsv")

# --- EffectorP ---
eff = {}
effp = os.path.join(anno_dir, "effectorp.txt")
if os.path.exists(effp):
    for line in open(effp):
        p = [x.strip() for x in line.split("\t")]
        if p and accre.search(p[0]): eff[acc_of(p[0])] = p[-1]

# --- DeepTMHMM ---
tm_regions = defaultdict(int)
tmp = os.path.join(anno_dir, "tmhmm_TMRs.gff3")
if os.path.exists(tmp):
    for line in open(tmp):
        m = re.match(r"#\s*(\S+)\s+Number of predicted TMRs:\s*(\d+)", line)
        if m: tm_regions[acc_of(m.group(1))] = int(m.group(2))

fam_of = dict(zip(*[pd.read_csv(mem_csv)[c] for c in ["acc","family"]]))
def dom_str(a):
    ds = sorted(member_domains.get(a, []), key=lambda d: d["start"])
    return " | ".join(f"{d['desc']}({d['start']}-{d['end']})" for d in ds)

rows = []
for a in accs:
    fp = fs_pdb.get(a); fa = fs_afdb.get(a); nd = len(member_domains.get(a, []))
    rows.append(dict(acc=a, family=fam_of.get(a, "unassigned"), effectorp=eff.get(a, "NA"),
        is_effector=("effector" in eff.get(a, "").lower() and "non" not in eff.get(a, "").lower()),
        n_TMR=tm_regions.get(a, 0), pfam_domains=dom_str(a), n_pfam_dom=nd, multi_domain=nd >= 2,
        interpro_entries="; ".join(sorted(ipr_entries.get(a, []))),
        pdb_hit=fp["target"] if fp else "", pdb_tm=round(float(fp["tm"]),3) if fp else np.nan,
        pdb_fident=round(float(fp["fident"]),3) if fp else np.nan,
        afdbsp_hit=fa["target"] if fa else "", afdbsp_tm=round(float(fa["tm"]),3) if fa else np.nan,
        has_known_fold=bool(fp) or bool(fa),
        has_any_domain=nd > 0 or len(ipr_entries.get(a, [])) > 0))
member = pd.DataFrame(rows)
member["novel"] = ~member.has_any_domain & ~member.has_known_fold
os.makedirs(os.path.dirname(out_mem), exist_ok=True)
member.to_csv(out_mem, index=False)

# cluster-level consensus annotation
clu = []
for fam, g in member[member.family != "unassigned"].groupby("family"):
    doms = [d for s in g.pfam_domains if s for d in s.split(" | ")]
    top = pd.Series([re.sub(r"\(.*", "", d) for d in doms]).value_counts()
    clu.append(dict(family=fam, n=len(g), n_effector=int(g.is_effector.sum()),
                    n_novel=int(g.novel.sum()), pct_novel=round(100*g.novel.mean(),1),
                    consensus_domain=top.index[0] if len(top) else "novel/unknown",
                    n_with_pdb=int((g.pdb_hit != "").sum())))
pd.DataFrame(clu).to_csv(out_clu, index=False)
print(f"annotate: {len(member)} proteins ({int(member.novel.sum())} novel), {len(clu)} clusters")

