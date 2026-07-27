"""Merge sequence evolution, structural conservation, and SASA on the family hub."""
import json, os, re
import numpy as np, pandas as pd

r4s_path = snakemake.input.r4s
structural_path = snakemake.input.structural
status_path = snakemake.input.status
sasa_csv = snakemake.input.sasa
ref_path = snakemake.input.ref            # {fam}.ref — the EXACT protein r4s numbered by
pdb_dir  = getattr(snakemake.input, "pdb_dir", None)
out_csv  = snakemake.output[0]
out_pdb  = snakemake.output[1] if len(snakemake.output) > 1 else None
strain   = snakemake.params.get("strain", "") if hasattr(snakemake.params, "get") else ""
fam      = snakemake.wildcards.fam
# reference accession = whatever the conservation rule pinned (MSA first header), so the
# r4s numbering and the SASA rows we filter are guaranteed to be the same protein.
ref = None
if os.path.exists(ref_path):
    ref = open(ref_path).read().strip() or None
if ref is None:
    famfile = os.path.join(os.path.dirname(os.path.dirname(r4s_path)), f"{fam}.members.txt")
    if os.path.exists(famfile):
        m = re.search(r"[A-Z]{2,3}\d{4,}\.\d+", os.path.basename(open(famfile).readline().strip()))
        ref = m.group(0) if m else None

def read_r4s(path):
    rows = []
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["resi", "aa", "cons"])
    for line in open(path):
        m = re.match(r"^\s*(\d+)\s+(\S)\s+(-?[\d.]+)", line.strip())
        if m: rows.append(dict(resi=int(m.group(1)), aa=m.group(2), cons=float(m.group(3))))
    return pd.DataFrame(rows)

r4s = read_r4s(r4s_path)
if "aa" in r4s:
    r4s = r4s.rename(columns={"aa": "aa_r4s"})
structural = (
    pd.read_csv(structural_path)
    if os.path.exists(structural_path) and os.path.getsize(structural_path) > 0
    else pd.DataFrame(columns=["resi", "aa"])
)
if "aa" in structural:
    structural = structural.rename(columns={"aa": "aa_struct"})
status = json.load(open(status_path)) if os.path.exists(status_path) else {}
sasa = pd.read_csv(sasa_csv)
s = sasa[sasa.acc == ref][["resi","aa","rel_sasa"]].copy()
s = s.rename(columns={"aa": "aa_sasa"})
s["rel_sasa"] = pd.to_numeric(s["rel_sasa"], errors="coerce")
if len(structural):
    m = structural.merge(s, on="resi", how="left", suffixes=("_struct", "_sasa"))
else:
    m = s.copy()
if len(r4s):
    m = m.merge(r4s, on="resi", how="left")
else:
    m["cons"] = np.nan
    m["aa_r4s"] = pd.NA
m["sequence_conservation"] = -pd.to_numeric(m["cons"], errors="coerce")
# Backward-compatible column name. It is explicitly sequence evolutionary
# conservation and remains NA when Rate4Site is not scientifically applicable.
m["conservation"] = m["sequence_conservation"]
m["rate4site_status"] = status.get("rate4site_status", "status_unavailable")
valid = m.dropna(subset=["rel_sasa","conservation"])
r = float(np.corrcoef(valid.conservation, valid.rel_sasa)[0,1]) if len(valid) > 2 else float("nan")
if "aa_r4s" in m and "aa_sasa" in m:
    compared = m.dropna(subset=["aa_r4s", "aa_sasa"])
    aa_match = float((compared.aa_r4s == compared.aa_sasa).mean()) if len(compared) else float("nan")
else:
    aa_match = float("nan")
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
m.to_csv(out_csv, index=False)

# conservation-colored reference PDB: write conservation into the B-factor column so the
# viewer (3Dmol) can color by conservation. resi -> conservation from the signature.
if out_pdb and pdb_dir and ref:
    import glob
    cons_by_resi = dict(zip(m.resi, m.sequence_conservation))
    cand = ([os.path.join(pdb_dir, f"{strain}_{ref}.pdb")] if strain else []) + \
           glob.glob(os.path.join(pdb_dir, f"*{ref}*.pdb"))
    src = next((p for p in cand if os.path.exists(p)), None)
    if src:
        out_lines = []
        for line in open(src):
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                try: ri = int(line[22:26])
                except ValueError: out_lines.append(line); continue
                b = cons_by_resi.get(ri, 0.0)
                out_lines.append(f"{line[:60]}{b:6.2f}{line[66:]}")
            else:
                out_lines.append(line)
        with open(out_pdb, "w") as fh: fh.writelines(out_lines)
    else:
        open(out_pdb, "w").write("")   # empty placeholder so the output exists
print(f"{fam} signature: ref={ref} n_aligned={len(m)} aa_match={aa_match:.3f} "
      f"cons_sasa_r={r:.3f} (n_valid={len(valid)}, "
      f"rate4site={status.get('rate4site_status', 'unknown')})")
