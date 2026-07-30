# SUSS Protein Atlas v5.1.3 reviewer-response handoff

## Decision

Release candidate status: **PASS after independent acceptance rerun**.

This release responds to the independent v5.1.2 review performed on the
300-protein *Magnaporthe oryzae* 70-15 cohort. The review independently
reproduced the whole-fold partition exactly (ARI 1.0) and identified latent
correctness, provenance, and workflow-ordering defects. v5.1.3 repairs those
defects without changing the validated clustering thresholds or the resulting
31-family/118-singleton partition.

## Acceptance dataset

The acceptance rerun reuses the reviewer's exact input bundle and effective
configuration:

```text
/home/claude/data/suss_reviewer_fix_clean_20260730
```

Final atlas:

```text
/home/claude/data/suss_reviewer_fix_clean_20260730/results/mor_suss_atlas.html
```

Execution log:

```text
/home/claude/suss_audits/20260730-v5.1.3-review-fixes/SUSS_v5.1.3_MOR_REBUILD.log
```

## Reviewer findings and resolution

| Finding | v5.1.3 resolution | Acceptance evidence |
|---|---|---|
| CR-1: missing Foldseek pairs stored as zero and included in `tm_cons_r` | TM matrices use `NA`; agreement uses only mutually finite Foldseek/US-align pairs; pair denominator is shown | F0 has 55/55 measured pairs; synthetic missing-pair tests pass |
| CR-9: full-length Foldseek lacked `--max-seqs` | Added configurable `clustering.max_seqs`, default 100000, to full-length Foldseek | Effective configs and parameter documentation updated |
| CR-3: empty evolutionary conservation advertised as valid | Conservation PDB/reference/button now require finite scores and bounds | Only 4/31 valid review families expose the mode; 27 legitimate `not_applicable` families do not |
| CR-5: unscored structural positions written as zero | Unscored residues use `-1.00` plus a PDB remark and are grey in the viewer; CA/alignment mismatch fails explicitly | F0 masked positions are `-1.00`; no false low-conservation zeros |
| CR-2: undeclared summary and pocket data channels | Split family CSVs are declared outputs; pocket result trees have declared JSON inventories; atlas assembly depends on all five outputs | Clean DAG rebuild completed 31/31 jobs |
| CR-4: ESM correlations joined by row position | ESM, conservation, and SASA are joined by residue ID | Offset-residue regression tests pass |
| CR-6: duplicated accession parsing | Five workflow consumers now use `v3_utils.protein_id` | Unit and full-cohort reruns pass |
| CR-7: hub penalised missing Foldseek pairs as zero | Hub maximises measured pair coverage, then mean measured TM, with deterministic tie-breaking | Measured/expected pair count is displayed; synthetic missingness test passes |
| Provenance and documentation gaps | VERSION/expected header/I/O contract updated; CI parses the Snakefile; release package carries `GIT_COMMIT` | Release archive validation records exact commit and version |

## Scientific result

The acceptance rerun retained the reviewer's independently reproduced result:

| Result | Value |
|---|---:|
| Input proteins | 300 |
| QC-passing proteins | 232 |
| Full-length families | 31 |
| Singletons | 118 |
| F0 members | 11 |
| F0 possible/measured TM pairs | 55/55 |
| Proteins with parsed pockets | 232/232 |
| Valid sequence-evolution conservation families | 4/31 |
| Invalid `nan` B-factor tokens in final atlas | 0 |

The low number of sequence-evolution conservation families is expected. The
analysis runs only when a sequence-related subgroup reaches the configured
minimum size; it is not inferred from structural-family membership.

## Interpretation retained from the review

Two reviewer observations are scientific limitations rather than software
defects:

1. Reciprocal coverage was non-binding in this cohort: every pair passing the
   TM threshold also passed coverage. Coverage remains an explicit configurable
   safeguard, but this dataset does not demonstrate its effect.
2. Foldseek and US-align agree strongly overall but differ near a hard threshold.
   A marginal family assignment should not be described as tool-independent
   merely because the family-level correlation is high.

## Test evidence

Local and 4070 validation:

```text
pytest: 65 passed
Python compile: passed
JavaScript syntax: passed
Snakefile parse/list-rules: passed
git diff --check: passed
M. oryzae reconstruction: 31/31 jobs, exit 0
```

The real 4070 rerun also caught two missing `re` imports in fpocket and
Rate4Site parsers during development. Both were repaired before the successful
continuation, which is recorded in the execution log.

Pocket inventories from the successful run:

```text
P2Rank files: 2320
fpocket files: 9014
proteins represented in merged pocket JSON: 232/232
```

## Claude acceptance checklist

1. Run `pytest -q`; expect 65 passing tests.
2. Parse `workflow/Snakefile` with the template configuration.
3. Inspect a family with `rate4site_status: not_applicable`; no evolutionary
   conservation control should be displayed.
4. Inspect one of the four valid families; the control and finite colour bounds
   should be present.
5. Download structural-conservation PDB data and confirm masked residues are
   `-1.00`, not `0.00`.
6. Inspect structural similarity: the Foldseek/US-align statistic must state
   how many mutually measured pairs were compared.
7. Confirm `family_summary_clustered.csv`,
   `family_summary_singletons.csv`, `p2rank.complete.json`, and
   `fpocket.complete.json` are declared workflow products.
8. Compare the deployed package `GIT_COMMIT` with the v5.1.3 GitHub release
   commit.

## Scope

This response validates the full-length reviewer findings and the shared
viewer/data infrastructure. It does not claim a new independent derivation of
the D-family partition; the original reviewer explicitly left that layer
outside its independent reimplementation.
