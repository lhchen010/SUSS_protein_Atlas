# Foldseek parameters in SUSS Protein Atlas v3

SUSS Atlas uses Foldseek for two different scientific questions. Their scores and thresholds
must not be interpreted as interchangeable.

## 1. Full-length fold families (`F`)

This layer asks whether two complete proteins share a comparable overall fold.

| Setting | Default | Meaning |
|---|---:|---|
| `clustering.alignment_type` | `1` | TM-align-style global structural alignment |
| `clustering.foldseek_tm` | `0.5` | Minimum symmetric TM score |
| `clustering.tm_symmetric` | `min` | Combine query- and target-normalized TM scores conservatively |
| `clustering.whole_fold_min_coverage` | `0.5` | Both proteins must align over at least half their length |
| `clustering.exhaustive_search` | `true` | Accurate all-vs-all search for atlas-scale datasets |
| `clustering.leiden_resolution` | `1.0` | Community granularity after structural edges are accepted |

The default edge is:

```text
min(query TM, target TM) >= 0.5
AND query coverage >= 0.5
AND target coverage >= 0.5
```

Raise the TM threshold to `0.7` for tighter fold families. Lowering coverage permits partial
matches to influence whole-protein families and should therefore be reported explicitly.

## 2. Local domain families (`D`)

This layer asks whether proteins share a local structural region, even when their complete
architectures differ. It uses local Foldseek 3Di+AA evidence and constructs a graph of matched
segments rather than whole proteins.

| Setting | Default | Meaning |
|---|---:|---|
| `domain_clustering.alignment_type` | `2` | Local 3Di+AA alignment |
| `domain_clustering.evalue` | `1e-3` | Statistical significance threshold |
| `domain_clustering.sensitivity` | `9.5` | Sensitive Foldseek search |
| `domain_clustering.min_probability` | `0.5` | Minimum Foldseek hit probability |
| `domain_clustering.min_aligned_residues` | `40` | Reject short motifs and noisy fragments |
| `domain_clustering.min_shorter_coverage` | `0.0` | Allow a domain embedded in two long proteins |
| `domain_clustering.interval_overlap` | `0.5` | Merge substantially overlapping hits into one segment |

Example: protein A has two domains and protein B contains only one of them. They may share a
`D` family for the matching region but remain in different `F` families. A protein may belong to
multiple `D` families.

## 3. TM score, 3Di+AA score, and coverage

- TM score describes global geometric similarity after normalization by protein length.
- Foldseek 3Di+AA search uses structural-alphabet and amino-acid evidence; e-value, probability,
  bit score, lDDT, and alignment length describe a local hit.
- Coverage is not a similarity score. It states how much of each protein participates in the
  alignment.
- A coverage value of `0.00` means no minimum coverage filter is requested. It does not mean the
  observed alignment covered zero residues; observed query and target coverage remain in result
  tables.

## 4. Recommended presets

| Goal | Full-length TM | Reciprocal coverage | Domain layer |
|---|---:|---:|---|
| General fungal secretome atlas | `0.5` | `0.5` | On, defaults |
| Tight fold families | `0.7` | `0.7` | On, defaults |
| Explore remote shared domains | `0.5` | `0.5` | On; inspect `D` hits carefully |
| Reproduce a paper | Paper value | Paper value | Record every changed field |

The portal exposes common F-family settings directly and keeps D-family controls under
**Advanced**. The effective config and tool versions are saved with every run.

## 5. Structural conservation coverage

`structural_conservation.pair_threshold` is passed to FoldMason `msa2lddtjson`. The default
`0.5` requires at least half of structure-pair subalignments to provide LDDT information for an
MSA column. Columns below that support level are retained in the CSV with a missing
`structural_lddt` value and displayed in grey. This prevents a residue present in only one or a
few structures from appearing strongly conserved.

This threshold is separate from F-family reciprocal coverage and from D-family local-alignment
coverage.

## 6. Scaling

For one strain or approximately 1,500 proteins, exhaustive all-vs-all search is reasonable.
For tens of thousands of proteins, use a separate discovery workflow based on Foldseek
prefiltering or clustering, then run the detailed atlas analyses only for selected groups.
The v3 `F`, `D`, and `S` identifiers and structure-search index are designed to support that
future multi-species workflow without redefining the current scientific layers.
