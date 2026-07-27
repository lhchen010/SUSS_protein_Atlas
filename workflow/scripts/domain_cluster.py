"""Build domain-segment families from local Foldseek 3Di+AA hits."""

import os

import igraph as ig
import leidenalg
import pandas as pd

from v3_utils import aggregate_domain_bridges, domain_segments


columns = [
    "query",
    "target",
    "qstart",
    "qend",
    "tstart",
    "tend",
    "alnlen",
    "qlen",
    "tlen",
    "evalue",
    "prob",
    "bits",
    "lddt",
    "fident",
]
try:
    hits = pd.read_csv(snakemake.input.tsv, sep="\t", names=columns)
except pd.errors.EmptyDataError:
    hits = pd.DataFrame(columns=columns)
segments, edges = domain_segments(
    hits,
    evalue_threshold=float(snakemake.params.evalue),
    probability_threshold=float(snakemake.params.probability),
    min_aligned_residues=int(snakemake.params.min_aligned),
    min_shorter_coverage=float(snakemake.params.min_coverage),
    interval_overlap=float(snakemake.params.interval_overlap),
)

if len(segments):
    graph = ig.Graph()
    graph.add_vertices(segments.segment_id.tolist())
    if len(edges):
        graph.add_edges(list(zip(edges.source, edges.target)))
        graph.es["weight"] = edges.prob.fillna(0).astype(float).tolist()
    partition = leidenalg.find_partition(
        graph,
        leidenalg.RBConfigurationVertexPartition,
        weights="weight" if len(edges) else None,
        resolution_parameter=float(snakemake.params.resolution),
        seed=int(snakemake.params.seed),
    )
    membership = dict(zip(graph.vs["name"], partition.membership))
    segments["community"] = segments.segment_id.map(membership)
    sizes = segments.community.value_counts()
    eligible = sizes[sizes >= int(snakemake.params.min_size)].index.tolist()
    ranked = sizes.loc[eligible].sort_values(ascending=False)
    labels = {community: f"D{index}" for index, community in enumerate(ranked.index)}
    segments["domain_family"] = segments.community.map(labels)
    segments = segments[segments.domain_family.notna()].copy()
    all_retained_edges = edges[
        edges.source.isin(set(segments.segment_id))
        & edges.target.isin(set(segments.segment_id))
    ].copy()
    segment_family = dict(zip(segments.segment_id, segments.domain_family))
    bridges = aggregate_domain_bridges(all_retained_edges, segment_family)
    all_retained_edges["domain_family"] = all_retained_edges.source.map(segment_family)
    edges = all_retained_edges[
        all_retained_edges.domain_family
        == all_retained_edges.target.map(segment_family)
    ].copy()
else:
    segments["community"] = pd.Series(dtype=int)
    segments["domain_family"] = pd.Series(dtype=str)
    edges["domain_family"] = pd.Series(dtype=str)
    bridges = aggregate_domain_bridges(edges, {})

families = []
for family, group in segments.groupby("domain_family"):
    family_edges = edges[edges.domain_family == family]
    families.append(
        {
            "domain_family": family,
            "n_segments": len(group),
            "n_proteins": group.acc.nunique(),
            "n_edges": len(family_edges),
            "mean_probability": family_edges.prob.mean() if len(family_edges) else None,
            "mean_lddt": family_edges.lddt.mean() if len(family_edges) else None,
            "mean_aligned_residues": (
                family_edges.alnlen.mean() if len(family_edges) else None
            ),
        }
    )
families = pd.DataFrame(
    families,
    columns=[
        "domain_family",
        "n_segments",
        "n_proteins",
        "n_edges",
        "mean_probability",
        "mean_lddt",
        "mean_aligned_residues",
    ],
)

for output in (
    snakemake.output.families,
    snakemake.output.members,
    snakemake.output.edges,
    snakemake.output.bridges,
):
    os.makedirs(os.path.dirname(output), exist_ok=True)
families.to_csv(snakemake.output.families, index=False)
segments[
    ["domain_family", "segment_id", "acc", "start", "end", "length", "community"]
].to_csv(snakemake.output.members, index=False)
edges.to_csv(snakemake.output.edges, index=False)
bridges.to_csv(snakemake.output.bridges, index=False)
print(
    f"domain families: {len(families)} families, {len(segments)} segments, "
    f"{segments.acc.nunique() if len(segments) else 0} proteins, "
    f"{len(bridges)} cross-family bridges"
)
