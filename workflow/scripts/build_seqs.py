"""seqs rule — canonical mature-sequence fasta. Uses the user-provided fasta where the
accession is present; fills any gap from the structure's own residues (PDB or mmCIF).
Guarantees every QC-passing accession has a sequence, so ESM/annotation never miss a protein.
"""
import os
import pandas as pd

from id_utils import extract_accession
from structure_io import list_structure_files, read_ca_sequence

pdb_dir = snakemake.input.pdb_dir
qc_csv = snakemake.input.qc
out_fa = snakemake.output[0]
provided = snakemake.params.provided
strain = snakemake.config.get("strain", {}).get("code", "")
custom_pattern = (
    str(snakemake.config.get("input", {}).get("accession_pattern", "") or "").strip() or None
)


def acc_of(s, *, strip_prefix=False):
    return extract_accession(
        s,
        strain=strain,
        custom_pattern=custom_pattern,
        strip_strain_prefix=strip_prefix,
    )


qc = pd.read_csv(qc_csv)
keep = set(qc[qc["pass"]].acc.astype(str)) if "pass" in qc.columns else set(qc.acc.astype(str))

# 1) provided fasta
seqs = {}
if provided and os.path.exists(provided):
    acc = None
    buf = []
    with open(provided, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if acc:
                    seqs[acc] = "".join(buf)
                acc = acc_of(line[1:], strip_prefix=False)
                buf = []
            else:
                buf.append(line)
    if acc:
        seqs[acc] = "".join(buf)
n_provided = len(seqs)

# 2) fill gaps from structure residues
acc2struct = {acc_of(p.name, strip_prefix=True): p for p in list_structure_files(pdb_dir)}
n_filled = 0
for acc in keep:
    if acc in seqs:
        continue
    p = acc2struct.get(acc)
    if not p:
        continue
    seq = read_ca_sequence(p)
    if seq:
        seqs[acc] = seq
        n_filled += 1

# 3) write only QC-passing accessions
os.makedirs(os.path.dirname(out_fa), exist_ok=True)
written = 0
with open(out_fa, "w", encoding="utf-8") as fh:
    for acc in sorted(keep):
        if acc in seqs:
            fh.write(f">{acc}\n{seqs[acc]}\n")
            written += 1
print(
    f"seqs.fasta: {written} sequences ({n_provided} provided, {n_filled} filled from structure); "
    f"{len(keep) - written} QC-pass accessions had neither"
)
