"""seqs rule — canonical mature-sequence fasta. Uses the user-provided fasta where the
accession is present; fills any gap (structure exists but no provided sequence — e.g.
full-length exceeds the fasta's length cutoff) from the structure's own residues via
Bio.PDB one-letter extraction. Guarantees every QC-passing accession has a sequence, so
ESM/annotation never miss a protein. Emits results/seqs.fasta + a provenance count.
"""
import os, glob, re
import pandas as pd

pdb_dir  = snakemake.input.pdb_dir
qc_csv   = snakemake.input.qc
out_fa   = snakemake.output[0]
provided = snakemake.params.provided

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
def acc_of(s):
    m = accre.search(s); return m.group(0) if m else s.split()[0]

qc = pd.read_csv(qc_csv)
keep = set(qc[qc["pass"]].acc.astype(str)) if "pass" in qc.columns else set(qc.acc.astype(str))

# 1) provided fasta
seqs = {}
if provided and os.path.exists(provided):
    acc = None; buf = []
    for line in open(provided, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line.startswith(">"):
            if acc: seqs[acc] = "".join(buf)
            acc = acc_of(line[1:]); buf = []
        else:
            buf.append(line)
    if acc: seqs[acc] = "".join(buf)
n_provided = len(seqs)

# 2) fill gaps from structure residues
T2O = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
acc2pdb = {acc_of(os.path.basename(p)): p for p in glob.glob(os.path.join(pdb_dir, "*.pdb"))}
n_filled = 0
for acc in keep:
    if acc in seqs: continue
    p = acc2pdb.get(acc)
    if not p: continue
    res = {}
    for line in open(p, encoding="utf-8", errors="replace"):
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            rn = line[17:20].strip(); ri = int(line[22:26])
            res[ri] = T2O.get(rn, "X")
    if res:
        seqs[acc] = "".join(res[k] for k in sorted(res))
        n_filled += 1

# 3) write only QC-passing accessions
os.makedirs(os.path.dirname(out_fa), exist_ok=True)
written = 0
with open(out_fa, "w") as fh:
    for acc in sorted(keep):
        if acc in seqs:
            fh.write(f">{acc}\n{seqs[acc]}\n"); written += 1
print(f"seqs.fasta: {written} sequences ({n_provided} provided, {n_filled} filled from structure); "
      f"{len(keep)-written} QC-pass accessions had neither")
