"""sasa_pocket rule — per-residue relative SASA for ALL proteins (incl. singletons)
via freesasa, plus pocket detection on each family reference (fpocket local + P2Rank).
freesasa runs everywhere; fpocket is LOCAL only (not installed on 4070); P2Rank needs
the java17 conda env. Emits sasa_all.csv + pockets.json. (recovered from sasa_all.py +
fpocket_all.py; hub-ref aware.)
"""
import os, glob, re, json, shutil, subprocess
import numpy as np, pandas as pd
import freesasa
freesasa.setVerbosity(freesasa.silent)

pdb_dir  = snakemake.input.pdb_dir
qc_csv   = snakemake.input.qc
out_sasa = snakemake.output.sasa
out_pock = snakemake.output.pockets
p2rank   = snakemake.params.p2rank
java_env = snakemake.params.java_env

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

# --- pockets on each family reference (P2Rank on 4070/java17; fpocket falls to local) ---
famdir = os.path.join(os.path.dirname(out_sasa), "families")
FPOCKET = shutil.which("fpocket")
pockets = {}
for ff in sorted(glob.glob(os.path.join(famdir, "*.members.txt"))):
    fam = os.path.basename(ff).split(".")[0]
    ref = None
    m = accre.search(open(ff).readline())
    ref = m.group(0) if m else None
    src = acc2pdb.get(ref)
    if not src: pockets[fam] = {"error": "no ref pdb", "ref": ref}; continue
    entry = {"ref": ref}
    # P2Rank
    wd = os.path.join(os.path.dirname(out_sasa), "p2rank", fam); os.makedirs(wd, exist_ok=True)
    ds = os.path.join(wd, f"{fam}.ds"); open(ds, "w").write(os.path.abspath(src) + "\n")
    try:
        env = dict(os.environ)
        jvm = f"/home/claude/.conda/envs/{java_env}/lib/jvm"
        if os.path.isdir(jvm): env["JAVA_HOME"] = jvm; env["PATH"] = jvm + "/bin:" + env["PATH"]
        subprocess.run([p2rank, "predict", ds, "-o", os.path.join(wd, "out")],
                       capture_output=True, env=env, timeout=600)
        pcsv = glob.glob(os.path.join(wd, "out", "*_predictions.csv"))
        if pcsv:
            pp = pd.read_csv(pcsv[0]); pp.columns = [c.strip() for c in pp.columns]
            if len(pp):
                top = pp.iloc[0]; rid = str(top.get("residue_ids", ""))
                resis = sorted({int(x.split("_")[-1]) for x in rid.split() if x.split("_")[-1].isdigit()})
                entry["p2rank"] = {"top_score": float(top.get("score", 0)), "n_pockets": len(pp),
                                   "lining_residues": resis}
    except Exception as e:
        entry["p2rank_error"] = str(e)[:100]
    # fpocket (local only)
    if FPOCKET:
        fwd = os.path.join(os.path.dirname(out_sasa), "fpocket", fam); os.makedirs(fwd, exist_ok=True)
        tgt = os.path.join(fwd, f"{ref}.pdb"); shutil.copy(src, tgt)
        subprocess.run([FPOCKET, "-f", tgt], capture_output=True)
        info = os.path.join(fwd, f"{ref}_out", f"{ref}_info.txt")
        if os.path.exists(info):
            txt = open(info).read()
            blk = re.findall(r"Pocket\s+(\d+)\s*:\s*\n\s*Score\s*:\s*([\-\d.]+)", txt)
            sc = {int(n): float(s) for n, s in blk}
            if sc:
                topn = max(sc, key=sc.get)
                atm = os.path.join(fwd, f"{ref}_out", "pockets", f"pocket{topn}_atm.pdb")
                resis = set()
                if os.path.exists(atm):
                    for line in open(atm):
                        if line.startswith(("ATOM", "HETATM")):
                            try: resis.add(int(line[22:26]))
                            except: pass
                entry["fpocket"] = {"top_score": round(sc[topn], 3), "n_pockets": len(sc),
                                    "lining_residues": sorted(resis)}
    pockets[fam] = entry
json.dump(pockets, open(out_pock, "w"))
print(f"pockets: {len([f for f in pockets if 'p2rank' in pockets[f] or 'fpocket' in pockets[f]])} families")

