"""rnaseq_family rule — slice ONE family's expression matrix out of the aggregate.

Root fix for the F0-clobber race: the old single rnaseq rule wrote every family's
{fam}_expression.csv as an UNDECLARED side effect, so Snakemake tracked no dependency
edge for them and a parallel rule touching a family dir could clobber one (F0 lost its
file this way). Here each {fam}_expression.csv is a DECLARED per-family output, so
Snakemake owns it, guarantees it before cards/atlas, and never lets another job race it.

Reads the aggregate results/rnaseq_expression.csv (acc + conditions, already in
biological order) and results/members.csv, emits members x conditions for this family.
A family whose members have no RNAseq locus gets a header-only file (card omits panel).
"""
import os
import pandas as pd

agg_path = snakemake.input.agg
members  = snakemake.input.members
fam      = snakemake.wildcards.fam
out      = snakemake.output[0]

agg = pd.read_csv(agg_path)
conds = [c for c in agg.columns if c != "acc"]
mem = pd.read_csv(members)
fam_accs = set(mem[mem.family == fam].acc.astype(str))

sub = agg[agg["acc"].astype(str).isin(fam_accs)] if len(agg) else agg.iloc[0:0]
os.makedirs(os.path.dirname(out), exist_ok=True)
sub.to_csv(out, index=False, columns=["acc"] + conds)
print(f"{fam} expression: {len(sub)}/{len(fam_accs)} members mapped; conditions={conds}")
