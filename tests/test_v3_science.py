from pathlib import Path
import sys

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "workflow" / "scripts"))

from v3_utils import domain_segments, merge_intervals, whole_fold_edges


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
