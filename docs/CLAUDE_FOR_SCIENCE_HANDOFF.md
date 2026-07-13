# SUSS Protein Atlas v1.0.1-rc1: Claude for Science handoff

Date: 2026-07-14 (Asia/Taipei)<br>
Branch: `fix/release-hardening`<br>
Production status: **not deployed; existing portal on port 8600 is unchanged**

## Purpose

This release candidate hardens the existing v1.0 workflow against silent tool failures,
incorrect optional-step states, annotation false positives, non-portable FoldTree execution,
and unsafe portal uploads. It deliberately preserves the validated clustering definition and
the 100-protein bundled-example family membership.

## Changes Claude should understand

1. External tools now accept either an absolute path or a command on `PATH`. A configured tool
   that is missing, exits unsuccessfully, or omits expected output fails the rule. Blank optional
   tools are recorded as `not_run`; failures are never represented by zeros or empty files.
2. US-align now verifies every expected within-family pair and uses the configured
   `clustering.tm_symmetric` mode. The same symmetric-TM function is used by matrices and
   cross-family edges.
3. FoldTree no longer contains `/home/user` runtime assumptions. Its Python, Snakemake,
   Foldseek, repository, and PATH prefix are configured. Families with fewer than three members
   are skipped instead of receiving zero-byte Newick files.
4. Annotation has component statuses and a tri-state `novel` field. Novelty is called only when
   InterProScan, pdb100 Foldseek, and AFDB Swiss-Prot Foldseek evidence all completed. Singleton
   proteins are excluded from cluster consensus. AFDB search headers now populate
   `afdbsp_name` with real protein names.
5. Optional analyses receive explicit enabled flags. Pocket analysis distinguishes P2Rank and
   fpocket status; a configured detector failure is fatal. Atlas renderer fallback is disabled by
   default. `cards` is explicitly required by `atlas`.
6. `results/used_config.yaml` now contains the effective Snakemake config, engine version,
   git commit when available, input SHA-256 values, PDB count, and configured/resolved tool state.
7. Singleton pLDDT and sequence-length summaries now use the current singleton, not the entire
   singleton group. Stale RNA-seq output is no longer reused when RNA-seq is disabled.
8. Portal hardening adds upload and expanded-archive limits, archive file-count limits, safe
   engine extraction, one-active-job default, CSRF-protected delete, UUID job IDs, and explicit
   `used_config.yaml` generation. Uploaded PDB names are normalized from accession to avoid
   duplicate or mismatched strain prefixes.

## 4070 staging evidence

Core staging root:

`/home/claude/suss_atlas_staging/20260714-hardening/repo`

Full-annotation staging root:

`/home/claude/suss_atlas_staging/20260714-annotation/repo`

Portal staging root and URL from the 4070 itself:

`/home/claude/suss_atlas_staging/20260714-portal`<br>
`http://127.0.0.1:8610`

Pre-change production backup:

`/home/claude/suss_backups/20260714-pre-hardening`

### Bundled 100-protein regression

The release candidate reproduced the existing clustering baseline exactly:

| Family | Members | Mean Foldseek TM | Core-SUSS % |
|---|---:|---:|---:|
| F0 | 28 | 0.594 | 92.8 |
| F1 | 26 | 0.669 | 13.7 |
| F2 | 23 | 0.603 | 79.0 |
| F3 | 7 | 0.605 | 54.5 |
| F4 | 2 | 0.642 | 100.0 |
| F5 | 2 | 0.907 | 0.0 |

Additional checks:

- 100 total proteins, 12 singletons, 536 classification rows.
- US-align nonzero pairs: F0 378/378, F1 325/325, F2 253/253, F3 21/21,
  F4 1/1, F5 1/1.
- FoldTree produced 12 nonempty files: three metrics each for F0-F3. F4/F5 were correctly
  skipped because `n=2`.
- P2Rank completed for all six families; fpocket was explicitly `not_run`; pocket status was
  therefore `partial`, not silently complete.
- Annotation-disabled output contained 100 `not_run` rows and zero non-null novelty calls.
- Re-running the completed core target returned `Nothing to be done`.

### Full annotation regression

- 100/100 rows: `annotation_status=complete`.
- InterPro, pdb100 Foldseek, AFDB Foldseek, EffectorP: 100/100 `complete`.
- 57 AFDB hits and 57 nonblank real protein names.
- Novel calls: 32 true, 68 false, 0 missing after complete evidence.
- Cluster annotation has exactly F0-F5; no synthetic `singleton` consensus row.
- DeepTMHMM was intentionally left blank for the successful run and is recorded as `not_run`.

### Portal regression

Successful jobs:

- `20260714-011053-cor-7941cd78`: 100 PDB + 100 FASTA, 6 families, atlas complete.
- `20260714-011259-smoke-bd7f1dcd`: same input with a different strain code; accession
  normalization succeeded, 6 families, atlas complete.

The latter job completed in about 31 seconds with optional expensive analyses disabled and
produced a 12,788,150-byte self-contained HTML atlas. HTTP checks returned 200 for status,
atlas, Excel summary, and config. An oversized request returned 413; invalid delete CSRF
returned 403.

## Required Claude validation

From the repository root on 4070:

```bash
source /home/user/miniconda3/etc/profile.d/conda.sh
conda activate suss
python -m unittest discover -s tests -v
python -m compileall -q workflow/scripts portal tests
snakemake --snakefile workflow/Snakefile \
  --configfile examples/config.example.yaml --cores 8 \
  results/example_suss_atlas.html results/family_summary.xlsx results/used_config.yaml
```

Claude should then verify:

- `members.csv` reproduces the family sizes and singleton count above.
- Every off-diagonal expected US-align pair is nonzero.
- F0-F3 have three nonempty Newick files each; F4/F5 have none.
- Disabled annotation yields `not_run` and nullable novelty, not false novelty.
- Full annotation yields 57 named AFDB hits and 32/68 novelty calls with the current databases.
- `used_config.yaml` reports `engine_version: 1.0.1-rc1`, `pdb_count: 100`, and input hashes.
- Atlas HTML renders normally; no fallback HTML message appears.
- A portal upload whose archived files already have another strain prefix is normalized and
  still matches the FASTA accessions.

## Known blocker: DeepTMHMM environment

The configured 4070 DeepTMHMM installation fails while generating ESM embeddings:

```text
AttributeError: 'Alphabet' object has no attribute 'unique_no_split_tokens'
```

The failure is in the licensed DeepTMHMM model/environment interaction with installed
`esm 2.0.0`, not in SUSS output parsing. The new wrapper correctly stops instead of reporting a
successful annotation. Before production deployment, Claude should repair or isolate the
DeepTMHMM environment, run the tool directly on `results/anno/clean.fasta`, then re-enable
`tools.deeptmhmm`. Do not patch the licensed third-party source inside this PR.

## Deliberately deferred work

This PR does not change the scientific BLAST/SUSS definition. Bidirectional coverage,
global sequence identity, and structural-alignment coverage remain a separate scientific-metric
change requiring a new baseline. It also does not replace the GenBank-specific accession model
with a manifest, enforce single-chain/single-model PDB input, add user authentication/isolation,
or migrate the portal away from Python's deprecated `cgi` module. Snakemake lint still reports
the pre-existing absence of per-rule `log` and `conda` directives.

## Deployment gate

Do not replace `/home/claude/suss_portal` until Claude has signed off on the checks above and
DeepTMHMM is either repaired or intentionally disabled in the production config. Deployment
should preserve `suss_portal_runs`, replace the engine tar and portal script atomically, restart
the service, and run one final upload smoke test on port 8600. Rollback material is in the backup
directory listed above.
