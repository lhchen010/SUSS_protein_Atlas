"""Write an auditable best-local-hit report without changing D-family membership."""

import os

import pandas as pd

from v3_utils import domain_match_diagnostics


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
    "alntmscore",
    "lddt",
    "fident",
    "qcov",
    "tcov",
]
try:
    hits = pd.read_csv(snakemake.input.hits, sep="\t", names=columns)
except pd.errors.EmptyDataError:
    hits = pd.DataFrame(columns=columns)
members = pd.read_csv(snakemake.input.members)
try:
    domain_members = pd.read_csv(snakemake.input.domain_members)
except pd.errors.EmptyDataError:
    domain_members = pd.DataFrame(columns=["acc"])

diagnostics = domain_match_diagnostics(
    hits,
    members["acc"].astype(str).tolist(),
    set(domain_members.get("acc", pd.Series(dtype=str)).astype(str)),
    evalue_threshold=float(snakemake.params.evalue),
    probability_threshold=float(snakemake.params.probability),
    min_aligned_residues=int(snakemake.params.min_aligned),
    min_shorter_coverage=float(snakemake.params.min_coverage),
    min_lddt=float(snakemake.params.min_lddt),
    min_alntm=float(snakemake.params.min_alntm),
    borderline_lddt_margin=float(snakemake.params.borderline_lddt_margin),
)
os.makedirs(os.path.dirname(snakemake.output.csv), exist_ok=True)
diagnostics.to_csv(snakemake.output.csv, index=False)
print(
    "domain diagnostics: "
    + ", ".join(
        f"{status}={count}"
        for status, count in diagnostics.status.value_counts().items()
    )
)
