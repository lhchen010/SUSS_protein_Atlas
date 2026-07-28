"""Infer sequence-homologous subgroups inside each structure-defined family."""

import os

import igraph as ig
import pandas as pd

from v3_utils import coverage, protein_id, select_blast_relationships


columns = [
    "query",
    "target",
    "pident",
    "alnlen",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]
blast = pd.read_csv(snakemake.input.blastp, sep="\t", names=columns)
blast["q"] = blast["query"].map(protein_id)
blast["t"] = blast["target"].map(protein_id)
blast = blast[blast.q != blast.t].copy()
blast["qcov"] = [
    coverage(aligned, length) for aligned, length in zip(blast.alnlen, blast.qlen)
]
blast["scov"] = [
    coverage(aligned, length) for aligned, length in zip(blast.alnlen, blast.slen)
]
blast["min_coverage"] = blast[["qcov", "scov"]].min(axis=1)
blast = select_blast_relationships(
    blast,
    evalue_threshold=float(snakemake.params.evalue),
    coverage_threshold=float(snakemake.params.coverage),
)

members = pd.read_csv(snakemake.input.members)
rows = []
edge_rows = []
for family, group in members[members.family != "singleton"].groupby("family"):
    accessions = sorted(group.acc.astype(str).unique())
    graph = ig.Graph()
    graph.add_vertices(accessions)
    selected = blast[
        blast.q.isin(accessions) & blast.t.isin(accessions)
    ].copy()
    if len(selected):
        graph.add_edges(list(zip(selected.q, selected.t)))
    components = graph.connected_components()
    ranked = sorted(
        components,
        key=lambda component: (-len(component), min(graph.vs[component]["name"])),
    )
    subgroup_for = {}
    for index, component in enumerate(ranked):
        subgroup = f"{family}.S{index}"
        for vertex in component:
            subgroup_for[graph.vs[vertex]["name"]] = subgroup
    for accession in accessions:
        subgroup = subgroup_for[accession]
        subgroup_size = sum(value == subgroup for value in subgroup_for.values())
        rows.append(
            {
                "acc": accession,
                "family": family,
                "sequence_subgroup": subgroup,
                "n_members": subgroup_size,
                "homology_status": (
                    "homologous_group" if subgroup_size >= 2 else "sequence_singleton"
                ),
            }
        )
    for row in selected.itertuples():
        edge_rows.append(
            {
                "family": family,
                "sequence_subgroup": subgroup_for[row.q],
                "q": row.q,
                "t": row.t,
                "pident": row.pident,
                "alnlen": row.alnlen,
                "qcov": row.qcov,
                "scov": row.scov,
                "min_coverage": row.min_coverage,
                "evalue": row.evalue,
                "bitscore": row.bitscore,
            }
        )

output = pd.DataFrame(
    rows,
    columns=[
        "acc",
        "family",
        "sequence_subgroup",
        "n_members",
        "homology_status",
    ],
)
edge_output = pd.DataFrame(
    edge_rows,
    columns=[
        "family",
        "sequence_subgroup",
        "q",
        "t",
        "pident",
        "alnlen",
        "qcov",
        "scov",
        "min_coverage",
        "evalue",
        "bitscore",
    ],
)
os.makedirs(os.path.dirname(snakemake.output.members), exist_ok=True)
output.to_csv(snakemake.output.members, index=False)
edge_output.to_csv(snakemake.output.edges, index=False)
print(
    f"sequence subgroups: {output.sequence_subgroup.nunique()} across "
    f"{output.family.nunique()} structural families"
)
