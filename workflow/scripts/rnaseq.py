"""rnaseq rule — per-family expression matrices from a strain RNAseq workbook.
Input xlsx has two sheets in the lab-standard format the preflight validator enforces:
'id_mapping' (headers: protein_accession, gene_id) and 'expression' (first column
gene_id, then one column per RNAseq sample). For each family, emit {fam}_expression.csv
= members x conditions
(replicates averaged per condition), member accession in the first column. Families whose
members have no RNAseq locus get an empty file (rule still succeeds; the card simply
omits the RNAseq panel). Produces results/families/{fam}/{fam}_expression.csv for all fams
plus results/rnaseq_expression.csv (all mapped members)."""
import os, re
import numpy as np, pandas as pd

xlsx     = snakemake.input.xlsx
members  = snakemake.input.members          # results/members.csv (acc,family,...)
out_all  = snakemake.output.all
fam_dir  = snakemake.params.fam_dir

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
xl = pd.ExcelFile(xlsx)
expr = xl.parse("expression")
# lab-standard sheet is 'id_mapping' with headers [protein_accession, gene_id];
# tolerate the older 'mapping'/[accession, locus] names for backward compatibility.
map_sheet = "id_mapping" if "id_mapping" in xl.sheet_names else "mapping"
mapping = xl.parse(map_sheet)
expr = expr.rename(columns={expr.columns[0]: "locus"})
acc_col = next((c for c in ("protein_accession", "accession") if c in mapping.columns), mapping.columns[0])
loc_col = next((c for c in ("gene_id", "locus") if c in mapping.columns), mapping.columns[1])
acc2loc = dict(zip(mapping[acc_col].astype(str), mapping[loc_col].astype(str)))

# collapse replicates to conditions: strip trailing .N / _N
sample_cols = [c for c in expr.columns if c != "locus"]
cond_of = {c: re.sub(r"[._]\d+$", "", str(c)) for c in sample_cols}
uniq = list(dict.fromkeys(cond_of.values()))   # order as they appear in the sheet
# biological ordering hint: control/mock -> in-vitro -> time course (1<3<7 DPI/HPI).
# falls back to sheet order for any label the pattern doesn't recognise.
def _cond_key(c):
    s = str(c).lower()
    base = 0
    if re.search(r"(^|[_\- ])(c|ctrl|control|mock)($|[_\- ])", s): base = 0
    elif "vh" in s or "vitro" in s or "vegetat" in s: base = 1
    else: base = 2
    m = re.search(r"(\d+)\s*(dpi|hpi|dai|h|d)\b", s) or re.search(r"(\d+)", s)
    t = int(m.group(1)) if m else 0
    return (base, t, str(c))
conds = sorted(uniq, key=_cond_key)
cond_mean = pd.DataFrame({cd: expr[[c for c in sample_cols if cond_of[c] == cd]].mean(axis=1)
                          for cd in conds})
cond_mean.insert(0, "locus", expr["locus"].astype(str))
loc2row = cond_mean.set_index("locus")

def resolve(loc):
    # RNAseq count tables often carry a transcript suffix (.t1) the mapping lacks
    if not loc:
        return None
    if loc in loc2row.index:
        return loc
    for suf in (".t1", ".1", ".t1.1"):
        if loc + suf in loc2row.index:
            return loc + suf
    return None

mem = pd.read_csv(members)
# Build ONLY the aggregate (acc + conditions for every mapped member of any family). The
# per-family {fam}_expression.csv files are produced by the separate rnaseq_family rule,
# which declares each as a tracked output (prevents the clobber race). We no longer write
# per-family files here as undeclared side effects.
rows_all = []
# EVERY mapped member — clustered families AND singletons (singletons are one of the two
# downloadable summary lists and must carry RNAseq too).
for a in mem.acc.astype(str).tolist():
    loc = resolve(acc2loc.get(a))
    if loc is not None:
        r = loc2row.loc[loc]
        rows_all.append(dict(acc=a, **{cd: float(r[cd]) for cd in conds}))

os.makedirs(os.path.dirname(out_all), exist_ok=True)
# de-dup on acc (a member appears once) preserving first occurrence
seen = set(); dedup = []
for row in rows_all:
    if row["acc"] not in seen:
        seen.add(row["acc"]); dedup.append(row)
pd.DataFrame(dedup, columns=["acc"] + conds).to_csv(out_all, index=False)
n_fam = mem[mem.family != "singleton"].family.nunique()
print(f"rnaseq aggregate: {n_fam} families + singletons, {len(dedup)} members with expression, conditions={conds}")
