"""Run configured annotation tools and emit explicit completion states."""

import os
import re
import shutil
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from runtime_utils import analysis_status, novel_call, resolve_executable, resolve_file


pdb_dir = snakemake.input.pdb_dir
seqs_fa = snakemake.input.seqs
mem_csv = snakemake.input.members
out_mem = snakemake.output.member
out_clu = snakemake.output.cluster
params = snakemake.params
enabled = bool(params.enabled)
anno_dir = os.path.join(os.path.dirname(out_mem), "anno")
os.makedirs(anno_dir, exist_ok=True)
log_path = os.path.join(anno_dir, "annotate_tools.log")

ACCRE = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")


def acc_of(value):
    match = ACCRE.search(str(value))
    return match.group(0) if match else str(value).split()[0]


def run_tool(cmd, name, *, timeout=5400, cwd=None, stdout_path=None, allow_output_on_error=None):
    with open(log_path, "a", encoding="utf-8") as log:
        log.write("\n$ " + " ".join(map(str, cmd)) + "\n")
        result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout)
        log.write(result.stdout)
        log.write(result.stderr)
    if stdout_path is not None:
        Path(stdout_path).write_text(result.stdout, encoding="utf-8")
    if result.returncode != 0 and not (allow_output_on_error and os.path.exists(allow_output_on_error)):
        raise RuntimeError(f"{name} failed with exit code {result.returncode}; see {log_path}")
    return result


members = pd.read_csv(mem_csv)
fam_of = dict(zip(members.acc.astype(str), members.family.astype(str)))
accs = []
with open(seqs_fa, encoding="utf-8", errors="replace") as handle:
    for line in handle:
        if line.startswith(">"):
            accs.append(acc_of(line[1:]))
accs = sorted(set(accs))

clean_fa = os.path.join(anno_dir, "clean.fasta")
with open(seqs_fa, encoding="utf-8", errors="replace") as source, open(clean_fa, "w", encoding="utf-8") as target:
    for line in source:
        target.write(line if line.startswith(">") else line.replace("*", ""))

statuses = {
    "interpro_status": "not_run",
    "foldseek_pdb_status": "not_run",
    "foldseek_afdb_status": "not_run",
    "effectorp_status": "not_run",
    "deeptmhmm_status": "not_run",
    "afdb_name_status": "not_run",
}

ipr_tsv = os.path.join(anno_dir, "interpro.tsv")
fs_pdb_tsv = os.path.join(anno_dir, "foldseek_pdb100.tsv")
fs_afdb_tsv = os.path.join(anno_dir, "foldseek_afdbsp.tsv")
effp_out = os.path.join(anno_dir, "effectorp.txt")
tm_gff = os.path.join(anno_dir, "tmhmm_TMRs.gff3")

if enabled:
    ips_raw = str(params.get("ips", "") or "").strip()
    if ips_raw:
        ips = resolve_executable(ips_raw, "InterProScan")
        if os.path.exists(ipr_tsv):
            os.remove(ipr_tsv)
        run_tool([ips, "-i", clean_fa, "-f", "TSV", "-o", ipr_tsv,
                  "-appl", "Pfam,CDD,Gene3D", "--cpu", "24", "--goterms"], "InterProScan")
        if not os.path.isfile(ipr_tsv):
            raise RuntimeError("InterProScan completed without interpro.tsv")
        statuses["interpro_status"] = "complete"

    fs_raw = str(params.get("fs", "") or "").strip()
    fsdb_raw = str(params.get("fsdb", "") or "").strip()
    if fs_raw or fsdb_raw:
        if not (fs_raw and fsdb_raw):
            raise ValueError("Foldseek annotation requires both tools.foldseek and tools.foldseek_db_dir")
        foldseek = resolve_executable(fs_raw, "Foldseek annotation")
        fsdb = str(Path(fsdb_raw).expanduser().resolve())
        with tempfile.TemporaryDirectory(prefix="suss_anno_") as work:
            query_db = os.path.join(work, "query")
            run_tool([foldseek, "createdb", pdb_dir, query_db], "Foldseek createdb", timeout=1800)
            for db_name, output, status_key in (
                ("pdb100", fs_pdb_tsv, "foldseek_pdb_status"),
                ("afdb_swissprot", fs_afdb_tsv, "foldseek_afdb_status"),
            ):
                database = os.path.join(fsdb, db_name)
                if not (os.path.exists(database) or os.path.exists(database + ".dbtype")):
                    raise FileNotFoundError(f"Configured Foldseek database not found: {database}")
                if os.path.exists(output):
                    os.remove(output)
                fmt = "query,target,theader,evalue,alntmscore,fident,alnlen,prob,qstart,qend,tstart,tend"
                run_tool([foldseek, "easy-search", query_db, database, output,
                          os.path.join(work, "tmp_" + db_name), "--format-output", fmt,
                          "-e", "0.01", "--max-seqs", "10", "--threads", "24"],
                         f"Foldseek {db_name}")
                if not os.path.isfile(output):
                    raise RuntimeError(f"Foldseek {db_name} completed without an output file")
                statuses[status_key] = "complete"
                if db_name == "afdb_swissprot":
                    statuses["afdb_name_status"] = "complete"

    effp_raw = str(params.get("effp", "") or "").strip()
    if effp_raw:
        effp = resolve_file(effp_raw, "EffectorP")
        effp_python = resolve_executable(params.get("effp_python", "python3"), "EffectorP Python")
        run_tool([effp_python, effp, "-i", clean_fa], "EffectorP", timeout=2400, stdout_path=effp_out)
        statuses["effectorp_status"] = "complete"

    tmhmm_raw = str(params.get("tmhmm", "") or "").strip()
    if tmhmm_raw:
        tmhmm = resolve_file(tmhmm_raw, "DeepTMHMM")
        tmhmm_python = resolve_executable(params.get("tmhmm_python", "python3"), "DeepTMHMM Python")
        tmout = os.path.join(anno_dir, "tmhmm")
        shutil.rmtree(tmout, ignore_errors=True)
        expected = os.path.join(tmout, "TMRs.gff3")
        run_tool([tmhmm_python, os.path.basename(tmhmm), "--fasta", os.path.abspath(clean_fa),
                  "--output-dir", os.path.abspath(tmout)], "DeepTMHMM", cwd=os.path.dirname(tmhmm),
                 timeout=3600, allow_output_on_error=expected)
        if not os.path.isfile(expected):
            raise RuntimeError("DeepTMHMM completed without TMRs.gff3")
        shutil.copy(expected, tm_gff)
        statuses["deeptmhmm_status"] = "complete"

name_map = {}
name_map_raw = str(params.get("afdb_name_map", "") or "").strip()
if enabled and name_map_raw:
    name_map_path = resolve_file(name_map_raw, "AFDB Swiss-Prot name map")
    sep = "\t" if name_map_path.endswith((".tsv", ".txt")) else ","
    mapping = pd.read_csv(name_map_path, sep=sep)
    required = {"target", "name"}
    if not required.issubset(mapping.columns):
        raise ValueError(f"AFDB name map must contain columns {sorted(required)}")
    name_map = dict(zip(mapping.target.astype(str), mapping.name.astype(str)))
    statuses["afdb_name_status"] = "complete"

member_domains = defaultdict(list)
ipr_entries = defaultdict(set)
if statuses["interpro_status"] == "complete":
    ipr = pd.read_csv(ipr_tsv, sep="\t", header=None,
                      names=["acc", "md5", "len", "db", "sig", "desc", "start", "end",
                             "evalue", "status", "date", "ipr", "ipr_desc", "go"],
                      usecols=range(14), on_bad_lines="skip")
    ipr["acc"] = ipr["acc"].map(acc_of)
    for _, row in ipr[ipr.db == "Pfam"].iterrows():
        try:
            member_domains[row.acc].append({"desc": row.desc, "start": int(row.start), "end": int(row.end)})
        except (TypeError, ValueError):
            pass
    for _, row in ipr.iterrows():
        if isinstance(row.ipr_desc, str) and row.ipr_desc not in ("-", ""):
            ipr_entries[row.acc].add(row.ipr_desc)

fscols = ["query", "target", "theader", "evalue", "tm", "fident", "alnlen", "prob", "qs", "qe", "ts", "te"]


def best_fs(path, status):
    if status != "complete":
        return {}
    frame = pd.read_csv(path, sep="\t", header=None, names=fscols)
    if frame.empty:
        return {}
    frame["query"] = frame["query"].map(acc_of)
    return frame.loc[frame.groupby("query").evalue.idxmin()].set_index("query").to_dict("index")


def protein_name(hit):
    if not hit:
        return ""
    target = str(hit["target"])
    header = str(hit.get("theader", "") or "").strip()
    if header.startswith(target):
        header = header[len(target):].strip()
    return name_map.get(target, header)


fs_pdb = best_fs(fs_pdb_tsv, statuses["foldseek_pdb_status"])
fs_afdb = best_fs(fs_afdb_tsv, statuses["foldseek_afdb_status"])

eff = {}
if statuses["effectorp_status"] == "complete":
    with open(effp_out, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            parts = [item.strip() for item in line.split("\t")]
            if parts and ACCRE.search(parts[0]):
                eff[acc_of(parts[0])] = parts[-1]

tm_regions = defaultdict(int)
if statuses["deeptmhmm_status"] == "complete":
    with open(tm_gff, encoding="utf-8", errors="replace") as handle:
        for line in handle:
            match = re.match(r"#\s*(\S+)\s+Number of predicted TMRs:\s*(\d+)", line)
            if match:
                tm_regions[acc_of(match.group(1))] = int(match.group(2))


def domain_string(acc):
    domains = sorted(member_domains.get(acc, []), key=lambda item: item["start"])
    return " | ".join(f"{item['desc']}({item['start']}-{item['end']})" for item in domains)


required_evidence = ("interpro_status", "foldseek_pdb_status", "foldseek_afdb_status")
annotation_status = analysis_status(enabled, statuses, required_evidence)
evidence_complete = annotation_status == "complete"
rows = []
novel_values = []
for acc in accs:
    pdb_hit = fs_pdb.get(acc)
    afdb_hit = fs_afdb.get(acc)
    domain_count = len(member_domains.get(acc, []))
    has_domain = domain_count > 0 or bool(ipr_entries.get(acc))
    has_fold = bool(pdb_hit) or bool(afdb_hit)
    target = afdb_hit["target"] if afdb_hit else ""
    call = eff.get(acc, "NA")
    rows.append({
        "acc": acc, "family": fam_of.get(acc, "unassigned"), "annotation_status": annotation_status,
        **statuses, "effectorp": call,
        "is_effector": "effector" in call.lower() and "non" not in call.lower(),
        "n_TMR": tm_regions.get(acc, 0), "pfam_domains": domain_string(acc),
        "n_pfam_dom": domain_count, "multi_domain": domain_count >= 2,
        "interpro_entries": "; ".join(sorted(ipr_entries.get(acc, []))),
        "pdb_hit": pdb_hit["target"] if pdb_hit else "",
        "pdb_tm": round(float(pdb_hit["tm"]), 3) if pdb_hit else np.nan,
        "pdb_fident": round(float(pdb_hit["fident"]), 3) if pdb_hit else np.nan,
        "afdbsp_hit": target, "afdbsp_name": protein_name(afdb_hit),
        "afdbsp_tm": round(float(afdb_hit["tm"]), 3) if afdb_hit else np.nan,
        "has_known_fold": has_fold, "has_any_domain": has_domain,
    })
    novel = novel_call(has_domain, has_fold, evidence_complete)
    novel_values.append(pd.NA if novel is None else novel)

member = pd.DataFrame(rows)
member["novel"] = pd.array(novel_values, dtype="boolean")
os.makedirs(os.path.dirname(out_mem), exist_ok=True)
member.to_csv(out_mem, index=False)

cluster_rows = []
clustered = member[~member.family.isin(["unassigned", "singleton"])]
for family, group in clustered.groupby("family"):
    domains = [domain for value in group.pfam_domains if value for domain in value.split(" | ")]
    top = pd.Series([re.sub(r"\(.*", "", domain) for domain in domains], dtype="object").value_counts()
    known_novel = group.novel.dropna()
    cluster_rows.append({
        "family": family, "n": len(group), "annotation_status": annotation_status,
        "n_effector": int(group.is_effector.sum()) if statuses["effectorp_status"] == "complete" else np.nan,
        "n_novel": int(known_novel.sum()) if len(known_novel) else np.nan,
        "pct_novel": round(100 * known_novel.mean(), 1) if len(known_novel) else np.nan,
        "consensus_domain": top.index[0] if len(top) else "unknown",
        "n_with_pdb": int((group.pdb_hit != "").sum()),
    })
pd.DataFrame(cluster_rows, columns=["family", "n", "annotation_status", "n_effector", "n_novel",
                                    "pct_novel", "consensus_domain", "n_with_pdb"]).to_csv(out_clu, index=False)
print(f"annotate: {len(member)} proteins; status={annotation_status}; {len(cluster_rows)} clusters")
