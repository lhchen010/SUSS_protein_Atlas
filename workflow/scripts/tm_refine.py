"""tm_refine rule — per-family US-align TM matrix, an ALGORITHM-INDEPENDENT cross-check of
the Foldseek TM used for clustering. Foldseek (3Di + its own TMalign module) builds the
families; US-align (rigid-body superposition, the TM-align successor, v2024xxxx) then
recomputes the true symmetric TM within each family. We do NOT re-cluster on this — it is a
validation/robustness layer: the card and family_summary show both, and disagreeing pairs
(|Foldseek TM − US-align TM| large) are flagged.

Runs US-align in -dir batch mode (one process does the family's all-vs-all), so cost is
O(m^2) within a family only (families are small; whole-set all-vs-all stays Foldseek's job).
Symmetric TM follows the same configured min/max/mean convention as clustering.

Output {fam}_TM_usalign.csv mirrors {fam}_TM.csv (labeled square matrix, diag=1.0).
Missing tools, structures, failed commands, and incomplete pair output block the rule.
"""
import os, re, glob, subprocess, tempfile, shutil
import numpy as np, pandas as pd
from runtime_utils import resolve_executable, symmetric_tm

famfile  = snakemake.input.famfile
pdb_dir  = snakemake.input.pdb_dir
out_tm   = snakemake.output.tm
fam      = snakemake.wildcards.fam
usalign  = snakemake.params.usalign          # absolute path to USalign binary
strain   = snakemake.params.get("code", "")  # filename prefix <code>_<acc>.pdb
sym_mode = snakemake.params.sym

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(s):
    m = accre.search(str(s)); return m.group(0) if m else str(s)

members = [acc_of(l) for l in open(famfile) if l.strip()]
members = list(dict.fromkeys(members))        # unique, keep order (member[0] = ref)

# resolve each member to a real PDB file (try <code>_<acc>.pdb, then *<acc>*.pdb)
def find_pdb(acc):
    cands = []
    if strain:
        cands.append(os.path.join(pdb_dir, f"{strain}_{acc}.pdb"))
    cands += glob.glob(os.path.join(pdb_dir, f"*{acc}*.pdb"))
    for c in cands:
        if os.path.exists(c):
            return c
    return None

pdbs = {a: find_pdb(a) for a in members}
have = [a for a in members if pdbs[a]]
missing = [a for a in members if not pdbs[a]]
if missing:
    raise FileNotFoundError(f"{fam}: missing PDB files for {', '.join(missing)}")
usalign = resolve_executable(usalign, "US-align")

TM = pd.DataFrame(0.0, index=members, columns=members)
for m in members:
    TM.loc[m, m] = 1.0

note = ""
if len(have) >= 2:
    # stage the family's structures into a temp dir with predictable names <acc>.pdb,
    # then run US-align -dir all-vs-all in one process.
    tmp = tempfile.mkdtemp(prefix=f"usaln_{fam}_")
    try:
        acc_by_stem = {}
        for a in have:
            stem = a                                   # temp filename stem = accession
            shutil.copy(pdbs[a], os.path.join(tmp, stem + ".pdb"))
            acc_by_stem[stem] = a
        listf = os.path.join(tmp, "list.txt")
        with open(listf, "w") as fh:
            fh.write("\n".join(acc_by_stem) + "\n")
        # -outfmt 2 = tab table: PDBchain1 PDBchain2 TM1 TM2 RMSD ID1 ID2 IDali L1 L2 Lali
        cmd = [usalign, "-dir", tmp + "/", listf, "-suffix", ".pdb", "-outfmt", "2"]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=3600, check=True)
        stem_re = re.compile(r"^(.*?)\.pdb(?::\S+)?$")
        for line in res.stdout.splitlines():
            if not line or line.startswith("#"):
                continue
            f = line.split("\t")
            if len(f) < 4:
                continue
            m1 = stem_re.match(f[0]); m2 = stem_re.match(f[1])
            if not (m1 and m2):
                continue
            a1 = acc_by_stem.get(m1.group(1)); a2 = acc_by_stem.get(m2.group(1))
            if a1 is None or a2 is None or a1 == a2:
                continue
            try:
                tm1, tm2 = float(f[2]), float(f[3])
            except ValueError:
                continue
            sym = symmetric_tm(tm1, tm2, sym_mode)
            TM.loc[a1, a2] = max(TM.loc[a1, a2], sym)
            TM.loc[a2, a1] = TM.loc[a1, a2]
        expected = len(have) * (len(have) - 1) // 2
        observed = int((TM.values[np.triu_indices(len(TM), 1)] > 0).sum())
        if observed != expected:
            raise RuntimeError(f"{fam}: US-align returned {observed}/{expected} expected pairs")
        note = f"US-align on {len(have)}/{len(members)} members ({sym_mode} symmetry)"
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
else:
    note = f"only {len(have)} member(s) — nothing to align"

os.makedirs(os.path.dirname(out_tm), exist_ok=True)
TM.round(3).to_csv(out_tm)
print(f"tm_refine {fam}: {note}; matrix {len(members)}x{len(members)} -> {out_tm}")
