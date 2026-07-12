"""composition rule — cluster-composition Excel. Two sheets:
cluster_summary (one row/family, with per-strain member counts n_<strain>, is_core,
is_lineage_specific placeholders for multi-strain) and member_composition (one row/
member with strain/species/host/pLDDT/length/n_neighbors). Strain parsed from the
<strain>_<accession> filename convention; single-strain here fills one strain column.
"""
import os, re
import pandas as pd

mem_csv  = snakemake.input.members
anno_csv = snakemake.input.annotation
out_xlsx = snakemake.output[0]

mem = pd.read_csv(mem_csv)
anno = pd.read_csv(anno_csv)
strain = snakemake.config["strain"]["code"]
species = snakemake.config["strain"].get("species", "")
host = snakemake.config["strain"].get("host", "")

mem = mem.merge(anno[["acc","is_effector","novel","pdb_hit"]], on="acc", how="left")
mem["strain"] = strain; mem["species"] = species; mem["host"] = host

# member_composition sheet
member_comp = mem[["acc","family","strain","species","host","plddt","length","deg"]].rename(
    columns={"plddt":"pLDDT","deg":"n_neighbors"})

# cluster_summary sheet (n_<strain> expands per strain when multi-strain)
fams = mem[mem.family != "singleton"]
rows = []
for f, g in fams.groupby("family"):
    row = dict(family=f, n_total=len(g))
    for s, gs in g.groupby("strain"): row[f"n_{s}"] = len(gs)
    row.update(n_species=g.species.nunique(),
               is_core=(g.strain.nunique() > 1),           # present in >1 strain (multi-strain)
               is_lineage_specific=(g.strain.nunique() == 1),
               host_range=";".join(sorted(set(g.host.dropna()))),
               n_effector=int(g.is_effector.fillna(False).sum()),
               n_novel=int(g.novel.fillna(False).sum()))
    rows.append(row)
cluster_summary = pd.DataFrame(rows).sort_values("n_total", ascending=False)

os.makedirs(os.path.dirname(out_xlsx), exist_ok=True)
with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xl:
    cluster_summary.to_excel(xl, sheet_name="cluster_summary", index=False)
    member_comp.to_excel(xl, sheet_name="member_composition", index=False)
print(f"composition: {len(cluster_summary)} clusters, {len(member_comp)} members -> {out_xlsx}")

