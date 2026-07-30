from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from v3_utils import (
    aggregate_domain_bridges,
    blast_identity_matrix,
    coverage,
    domain_sequence_records,
    domain_segments,
    foldmason_column_scores,
    merge_intervals,
    relationship_components,
    select_blast_relationships,
    whole_fold_edges,
)


def test_whole_fold_edges_require_tm_and_reciprocal_coverage():
    table = pd.DataFrame(
        [
            {
                "query": "A.pdb",
                "target": "B.pdb",
                "alntmscore": 0.7,
                "qtmscore": 0.8,
                "ttmscore": 0.7,
                "lddt": 0.8,
                "fident": 0.2,
                "alnlen": 80,
                "qlen": 100,
                "tlen": 100,
                "evalue": 1e-6,
                "bits": 90,
            },
            {
                "query": "A.pdb",
                "target": "C.pdb",
                "alntmscore": 0.8,
                "qtmscore": 0.8,
                "ttmscore": 0.8,
                "lddt": 0.8,
                "fident": 0.2,
                "alnlen": 30,
                "qlen": 100,
                "tlen": 100,
                "evalue": 1e-6,
                "bits": 90,
            },
        ]
    )

    edges = whole_fold_edges(
        table, tm_threshold=0.5, coverage_threshold=0.5, symmetry="min"
    )

    assert list(zip(edges.q, edges.t)) == [("A", "B")]
    assert edges.iloc[0].min_coverage == 0.8


def test_blast_relationship_uses_best_hit_that_passes_coverage():
    hits = pd.DataFrame([
        {
            "q": "A", "t": "B", "evalue": 1e-30, "bitscore": 200,
            "min_coverage": 0.2, "pident": 60,
        },
        {
            "q": "A", "t": "B", "evalue": 1e-20, "bitscore": 180,
            "min_coverage": 0.8, "pident": 40,
        },
    ])

    selected = select_blast_relationships(
        hits, evalue_threshold=1e-3, coverage_threshold=0.5
    )

    assert len(selected) == 1
    assert selected.iloc[0].evalue == 1e-20
    assert selected.iloc[0].pident == 40


def test_alignment_coverage_is_bounded_to_probability_range():
    assert coverage(110, 100) == 1.0
    assert coverage(-5, 100) == 0.0
    assert pd.isna(coverage(10, 0))


def test_interval_merging_preserves_distinct_domains():
    merged = merge_intervals([(10, 80), (15, 75), (180, 250), (190, 245)])

    assert len(merged) == 2
    assert merged[0][1] < merged[1][0]


def test_domain_segments_allow_a_local_domain_in_long_proteins():
    table = pd.DataFrame(
        [
            {
                "query": "A.pdb",
                "target": "B.pdb",
                "qstart": 20,
                "qend": 89,
                "tstart": 300,
                "tend": 369,
                "alnlen": 70,
                "qlen": 500,
                "tlen": 600,
                "evalue": 1e-8,
                "prob": 0.95,
                "bits": 100,
                "lddt": 0.75,
                "fident": 0.1,
            }
        ]
    )

    segments, edges = domain_segments(
        table,
        evalue_threshold=1e-3,
        probability_threshold=0.5,
        min_aligned_residues=40,
        min_shorter_coverage=0.0,
    )

    assert set(segments.acc) == {"A", "B"}
    assert len(edges) == 1
    assert edges.iloc[0].alnlen == 70


def test_domain_segments_reject_low_local_lddt_even_with_high_probability():
    table = pd.DataFrame(
        [
            {
                "query": "A.pdb",
                "target": "B.pdb",
                "qstart": 20,
                "qend": 89,
                "tstart": 300,
                "tend": 369,
                "alnlen": 70,
                "qlen": 500,
                "tlen": 600,
                "evalue": 1e-12,
                "prob": 0.99,
                "bits": 120,
                "alntmscore": 0.4,
                "lddt": 0.35,
                "fident": 0.1,
            }
        ]
    )

    segments, edges = domain_segments(
        table,
        evalue_threshold=1e-3,
        probability_threshold=0.5,
        min_aligned_residues=40,
        min_shorter_coverage=0.0,
        min_lddt=0.5,
    )

    assert segments.empty
    assert edges.empty


def test_domain_coverage_is_bounded_and_cross_family_bridges_are_aggregated():
    edges = pd.DataFrame([
        {
            "source": "A:1-50", "target": "B:1-50", "prob": 0.9,
            "lddt": 0.7, "alnlen": 55,
        },
        {
            "source": "A:1-50", "target": "C:1-50", "prob": 0.8,
            "lddt": 0.6, "alnlen": 45,
        },
    ])
    bridges = aggregate_domain_bridges(
        edges, {"A:1-50": "D0", "B:1-50": "D1", "C:1-50": "D1"}
    )

    assert len(bridges) == 1
    assert bridges.iloc[0].source_family == "D0"
    assert bridges.iloc[0].target_family == "D1"
    assert bridges.iloc[0].n_edges == 2
    assert round(bridges.iloc[0].mean_lddt, 6) == 0.65


def test_foldmason_column_scores_turn_unscored_columns_into_missing_values():
    scores = foldmason_column_scores(
        {"scores": [0.8, -1, None, 0.2]}, length=4
    )

    assert scores[0] == 0.8
    assert pd.isna(scores[1])
    assert pd.isna(scores[2])
    assert scores[3] == 0.2


def test_domain_sequences_are_cropped_and_keep_exact_segment_ids():
    members = pd.DataFrame(
        [
            {
                "domain_family": "D0",
                "segment_id": "A:2-5",
                "acc": "A",
                "start": 2,
                "end": 5,
            },
            {
                "domain_family": "D0",
                "segment_id": "B:1-3",
                "acc": "B",
                "start": 1,
                "end": 3,
            },
        ]
    )

    records, manifest = domain_sequence_records(
        members, {"A": "ABCDEFG", "B": "MNOP"}
    )

    assert records == {"A:2-5": "BCDE", "B:1-3": "MNO"}
    assert manifest.set_index("segment_id").loc["A:2-5", "parent_length"] == 7


def test_same_parent_domain_segments_are_not_related_without_segment_hit():
    components = relationship_components(
        ["A:1-50", "A:100-150", "B:5-55"],
        pd.DataFrame([{"q": "A:1-50", "t": "B:5-55"}]),
    )

    assert components == [["A:1-50", "B:5-55"], ["A:100-150"]]


def test_domain_blast_identity_uses_exact_segments_not_parent_accessions():
    matrix = blast_identity_matrix(
        pd.DataFrame(
            [
                {
                    "q": "A:1-50",
                    "t": "B:5-55",
                    "pident": 42.0,
                },
                {
                    "q": "A:100-150",
                    "t": "B:5-55",
                    "pident": 18.0,
                },
            ]
        ),
        ["A:1-50", "A:100-150", "B:5-55"],
    )

    assert matrix == [
        [1.0, 0.0, 0.42],
        [0.0, 1.0, 0.18],
        [0.42, 0.18, 1.0],
    ]
