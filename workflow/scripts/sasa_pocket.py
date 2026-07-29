"""SASA and pocket detection for every QC-passing protein.

Pocket predictions are keyed by accession and full-family keys alias the relevant
hub record. Full-length and domain workbenches therefore reuse one prediction
instead of rerunning a surface method on cropped domain coordinates.
FreeSASA runs in the workflow environment; P2Rank may use a separate Java
environment. Emits sasa_all.csv + pockets.json.
"""
import os, glob, re, json, shutil, subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np, pandas as pd
import freesasa
from runtime_utils import resolve_executable
freesasa.setVerbosity(freesasa.silent)

pdb_dir  = snakemake.input.pdb_dir
qc_csv   = snakemake.input.qc
out_sasa = snakemake.output.sasa
out_pock = snakemake.output.pockets
p2rank   = snakemake.params.p2rank
p2rank_profile = str(snakemake.params.p2rank_profile or "").strip().lower()
fpocket  = snakemake.params.fpocket
java_env = snakemake.params.java_env
enabled  = bool(snakemake.params.enabled)

MAXASA = {'A':129.,'R':274.,'N':195.,'D':193.,'C':167.,'E':223.,'Q':225.,'G':104.,'H':224.,
          'I':197.,'L':201.,'K':236.,'M':224.,'F':240.,'P':159.,'S':155.,'T':172.,'W':285.,'Y':263.,'V':174.}
T2O = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I',
       'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(fn): 
    m = accre.search(os.path.basename(fn)); return m.group(0) if m else os.path.basename(fn)[:-4]

acc2pdb = {acc_of(p): p for p in glob.glob(os.path.join(pdb_dir, "*.pdb"))}
qc = pd.read_csv(qc_csv); keep = set(qc[qc["pass"]].acc) if "pass" in qc.columns else set(qc.acc)

rows = []
for acc, p in acc2pdb.items():
    if acc not in keep: continue
    st = freesasa.Structure(p); ra = freesasa.calc(st).residueAreas()
    for ch in ra:
        for resi in ra[ch]:
            r = ra[ch][resi]; aa = T2O.get(r.residueType, 'X'); mx = MAXASA.get(aa, np.nan)
            rows.append(dict(acc=acc, resi=int(resi), aa=aa, sasa=round(r.total, 2),
                             rel_sasa=round(r.total/mx, 4) if mx else ""))
os.makedirs(os.path.dirname(out_sasa), exist_ok=True)
pd.DataFrame(rows).to_csv(out_sasa, index=False)
print(f"SASA: {len({r['acc'] for r in rows})} proteins, {len(rows)} residues")

# --- pockets on every protein; family keys alias their reference protein ---
famdir = os.path.join(os.path.dirname(out_sasa), "families")
members = pd.read_csv(snakemake.input.members)
P2RANK = resolve_executable(p2rank, "P2Rank") if enabled and str(p2rank).strip() else None
FPOCKET = resolve_executable(fpocket, "fpocket") if enabled and str(fpocket).strip() else None
CONDA = resolve_executable(snakemake.params.conda, "conda") if enabled and java_env else None
if enabled and not (P2RANK or FPOCKET):
    raise ValueError("steps.pocket is enabled but neither P2Rank nor fpocket is configured")
pockets = {}
previous_pockets = {}
if os.path.isfile(out_pock):
    try:
        previous_pockets = json.load(open(out_pock))
    except (OSError, ValueError):
        previous_pockets = {}
family_refs = {}
for ff in sorted(glob.glob(os.path.join(famdir, "*.members.txt"))):
    fam = os.path.basename(ff).split(".")[0]
    m = accre.search(open(ff).readline())
    ref = m.group(0) if m else None
    if ref:
        family_refs[fam] = ref
targets = sorted(keep)
ref_storage = {
    ref: family
    for family, ref in family_refs.items()
    if (
        os.path.isdir(os.path.join(os.path.dirname(out_sasa), "p2rank", family))
        or os.path.isdir(os.path.join(os.path.dirname(out_sasa), "fpocket", family))
    )
}
cached_profiles = {}
for key, value in previous_pockets.items():
    if not isinstance(value, dict):
        continue
    storage_key = str(value.get("storage_key") or key)
    profile = str(value.get("p2rank_profile") or "").strip().lower()
    if profile:
        cached_profiles[storage_key] = profile

def predict_pockets(ref):
    storage_key = ref_storage.get(ref, ref)
    src = acc2pdb.get(ref)
    if not src:
        return ref, {
            "error": "no ref pdb",
            "ref": ref,
            "storage_key": storage_key,
        }
    entry = {
        "ref": ref,
        "storage_key": storage_key,
        "pocket_status": "not_run" if not enabled else ("complete" if P2RANK and FPOCKET else "partial"),
        "p2rank_status": "pending" if P2RANK else "not_run",
        "p2rank_profile": p2rank_profile or "default",
        "fpocket_status": "pending" if FPOCKET else "not_run",
    }
    if not enabled:
        return ref, entry
    # P2Rank
    if P2RANK:
        wd = os.path.join(os.path.dirname(out_sasa), "p2rank", storage_key)
        out_dir = os.path.join(wd, "out")
        profile_marker = os.path.join(wd, "profile.txt")
        os.makedirs(wd, exist_ok=True)
        pcsv = sorted(glob.glob(os.path.join(out_dir, "*_predictions.csv")))
        marker_profile = (
            open(profile_marker).read().strip().lower()
            if os.path.isfile(profile_marker)
            else cached_profiles.get(storage_key, "")
        )
        if marker_profile != (p2rank_profile or "default"):
            pcsv = []
        if not pcsv:
            shutil.rmtree(out_dir, ignore_errors=True)
            ds = os.path.join(wd, f"{storage_key}.ds")
            open(ds, "w").write(os.path.abspath(src) + "\n")
            p2cmd = [P2RANK, "predict", ds, "-o", out_dir]
            if p2rank_profile and p2rank_profile != "default":
                p2cmd.extend(["-c", p2rank_profile])
            if CONDA:
                p2cmd = [CONDA, "run", "-n", java_env, *p2cmd]
            subprocess.run(
                p2cmd, capture_output=True, text=True, timeout=600, check=True
            )
            pcsv = sorted(glob.glob(os.path.join(out_dir, "*_predictions.csv")))
        open(profile_marker, "w").write((p2rank_profile or "default") + "\n")
        exact = [
            path for path in pcsv
            if (
                os.path.basename(path).endswith(f"{ref}.pdb_predictions.csv")
                or os.path.basename(path).endswith(f"{ref}_predictions.csv")
                or os.path.basename(path).endswith(
                    f"{storage_key}.pdb_predictions.csv"
                )
                or os.path.basename(path).endswith(
                    f"{storage_key}_predictions.csv"
                )
            )
        ]
        if len(exact) != 1:
            raise RuntimeError(
                f"{ref}: expected one P2Rank predictions CSV for {ref}, found {len(exact)}"
            )
        pp = pd.read_csv(exact[0]); pp.columns = [c.strip() for c in pp.columns]
        entry["p2rank_status"] = "complete"
        if len(pp):
            all_pockets = []
            for idx, pred in pp.iterrows():
                rid = str(pred.get("residue_ids", ""))
                resis = sorted({int(x.split("_")[-1]) for x in rid.split()
                                if x.split("_")[-1].isdigit()})
                all_pockets.append({
                    "pocket_id": int(pred.get("rank", idx + 1)),
                    "score": float(pred.get("score", 0)),
                    "probability": (
                        float(pred.get("probability"))
                        if pd.notna(pred.get("probability"))
                        else None
                    ),
                    "lining_residues": resis,
                })
            top = max(all_pockets, key=lambda p: p["score"])
            entry["p2rank"] = {"top_score": top["score"],
                               "top_probability": top["probability"],
                               "n_pockets": len(all_pockets),
                               "lining_residues": top["lining_residues"], "pockets": all_pockets}
    # fpocket (local only)
    if FPOCKET:
        fwd = os.path.join(os.path.dirname(out_sasa), "fpocket", storage_key); os.makedirs(fwd, exist_ok=True)
        existing_info = sorted(glob.glob(os.path.join(fwd, "*_out", "*_info.txt")))
        if existing_info:
            info = existing_info[0]
            output_stem = os.path.basename(info).removesuffix("_info.txt")
        else:
            tgt = os.path.join(fwd, f"{ref}.pdb"); shutil.copy(src, tgt)
            subprocess.run(
                [FPOCKET, "-f", tgt],
                capture_output=True,
                text=True,
                timeout=600,
                check=True,
            )
            info = os.path.join(fwd, f"{ref}_out", f"{ref}_info.txt")
            output_stem = ref
        if os.path.exists(info):
            txt = open(info).read()
            blk = re.findall(r"Pocket\s+(\d+)\s*:\s*\n\s*Score\s*:\s*([\-\d.]+)", txt)
            sc = {int(n): float(s) for n, s in blk}
            if sc:
                topn = max(sc, key=sc.get)
                all_pockets = []
                for pocket_id, score in sorted(sc.items()):
                    atm = os.path.join(
                        fwd,
                        f"{output_stem}_out",
                        "pockets",
                        f"pocket{pocket_id}_atm.pdb",
                    )
                    resis = set()
                    if os.path.exists(atm):
                        for line in open(atm):
                            if line.startswith(("ATOM", "HETATM")):
                                try: resis.add(int(line[22:26]))
                                except (ValueError, IndexError): pass
                    all_pockets.append({"pocket_id": pocket_id, "score": round(score, 3),
                                        "lining_residues": sorted(resis)})
                top = next(p for p in all_pockets if p["pocket_id"] == topn)
                entry["fpocket"] = {"top_score": round(sc[topn], 3), "n_pockets": len(sc),
                                    "lining_residues": top["lining_residues"], "pockets": all_pockets}
        entry["fpocket_status"] = "complete"
    return ref, entry


workers = max(1, int(getattr(snakemake, "threads", 1)))
with ThreadPoolExecutor(max_workers=workers) as executor:
    futures = {executor.submit(predict_pockets, ref): ref for ref in targets}
    for future in as_completed(futures):
        ref, entry = future.result()
        pockets[ref] = entry
for family, ref in family_refs.items():
    if ref in pockets:
        pockets[family] = dict(pockets[ref])
json.dump(pockets, open(out_pock, "w"))
singletons = int((members.family == "singleton").sum())
complete = len([
    accession for accession in targets
    if "p2rank" in pockets.get(accession, {})
    or "fpocket" in pockets.get(accession, {})
])
print(
    f"pockets: {complete}/{len(targets)} proteins "
    f"({len(family_refs)} full-family aliases, {singletons} protein singletons)"
)
