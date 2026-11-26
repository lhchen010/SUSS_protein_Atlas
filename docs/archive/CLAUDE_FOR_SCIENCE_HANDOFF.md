# SUSS Protein Atlas v1.0.1: Claude for Science handoff

Date: 2026-07-14 (Asia/Taipei)<br>
Release: `main`, tag `v1.0.1`<br>
Production status: **deployed and post-deployment validation passed**

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
- The initial P2Rank-only regression completed for all six families and correctly recorded
  fpocket as `not_run`; the final dual-detector regression is documented below.
- Annotation-disabled output contained 100 `not_run` rows and zero non-null novelty calls.
- Re-running the completed core target returned `Nothing to be done`.

### Full annotation regression

- 100/100 rows: `annotation_status=complete`.
- InterPro, pdb100 Foldseek, AFDB Foldseek, EffectorP: 100/100 `complete`.
- 57 AFDB hits and 57 nonblank real protein names.
- Novel calls: 32 true, 68 false, 0 missing after complete evidence.
- Cluster annotation has exactly F0-F5; no synthetic `singleton` consensus row.
- DeepTMHMM completed for 100/100 proteins through the compatibility launcher. It identified
  five proteins with transmembrane regions, 11 regions total, and a maximum of five in one protein.
- P2Rank and fpocket both completed for F0-F5. fpocket found 5, 7, 9, 5, 12, and 8 pockets
  respectively; all six family records have `pocket_status=complete`.

### Portal regression

Successful jobs:

- `20260714-011053-cor-7941cd78`: 100 PDB + 100 FASTA, 6 families, atlas complete.
- `20260714-011259-smoke-bd7f1dcd`: same input with a different strain code; accession
  normalization succeeded, 6 families, atlas complete.
- `20260714-013818-cor-726cbc90`: final v1.0.1 packaged engine with classification,
  conservation, dual-pocket analysis, FoldTree, complete annotation, and DeepTMHMM enabled;
  50 workflow steps completed in 121 seconds and produced a 13,513,416-byte atlas.

The latter job completed in about 31 seconds with optional expensive analyses disabled and
produced a 12,788,150-byte self-contained HTML atlas. HTTP checks returned 200 for status,
atlas, Excel summary, and config. An oversized request returned 413; invalid delete CSRF
returned 403.

Production ESM job `20260714-014753-cor-89d46131` completed 28/28 workflow steps in about
95 seconds. ESM-Scan scored all six family references: 1,028 residue positions total
(67/229/208/88/212/224 for F0-F5), with non-null scores for all 20 amino-acid substitutions at
every position. Its atlas completed successfully.

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
- `used_config.yaml` reports `engine_version: 1.0.1`, `pdb_count: 100`, and input hashes.
- Atlas HTML renders normally; no fallback HTML message appears.
- A portal upload whose archived files already have another strain prefix is normalized and
  still matches the FASTA accessions.

## DeepTMHMM compatibility environment

The original 4070 DeepTMHMM command failed while generating ESM embeddings:

```text
AttributeError: 'Alphabet' object has no attribute 'unique_no_split_tokens'
```

The licensed release requires `fair-esm==0.4.0`; the general user environment provided
`esm 2.0.0`. The 4070 deployment now uses `/home/claude/suss_tools/bin/deeptmhmm-python`,
which isolates `fair-esm 0.4.0` and maps DeepTMHMM's retired Matplotlib style name to its current
equivalent. The licensed third-party source is unchanged. Claude should verify the sample and
100-protein results above before approving a future environment upgrade.

FoldTree also now uses an isolated package work directory per family. This prevents nested
Snakemake metadata collisions during parallel execution and preserves the upstream package's
dangling test-data symlink instead of trying to dereference it. Each family writes
`foldtree_subworkflow.log` for diagnosis.

## Deliberately deferred work

This PR does not change the scientific BLAST/SUSS definition. Bidirectional coverage,
global sequence identity, and structural-alignment coverage remain a separate scientific-metric
change requiring a new baseline. It also does not replace the GenBank-specific accession model
with a manifest, enforce single-chain/single-model PDB input, add user authentication/isolation,
or migrate the portal away from Python's deprecated `cgi` module. Snakemake lint still reports
the pre-existing absence of per-rule `log` and `conda` directives.

## Production deployment result

Production was deployed to `/home/claude/suss_portal` on 2026-07-14. Existing run history was
preserved. The active portal is reachable over Tailscale at:

`http://100.80.77.29:8600`

Deployed identity:

- Engine version: `1.0.1`
- Source: GitHub `main`, tag `v1.0.1`; the packaged tag commit is embedded in `GIT_COMMIT` and
  appears in every new run's `results/used_config.yaml`
- Engine tar SHA-256: `f7b2d88eefa67013add9b2bdcecfd21b4545462471894055804a6a03198f1cb5`
- DeepTMHMM compatibility directory: `/home/claude/suss_tools/deeptmhmm-compat`
- DeepTMHMM launcher: `/home/claude/suss_tools/bin/deeptmhmm-python`

The production engine completed job `20260714-014252-cor-54742061` with 100 structures,
100 sequences, 6 families, and 50/50 workflow steps in about 129 seconds. Classification,
conservation, US-align, FoldTree, annotation, DeepTMHMM, fpocket, and P2Rank were enabled; ESM
was left off for this deployment-focused regression. Results matched the staging metrics above:

- family sizes 28/26/23/7/2/2 plus 12 singletons;
- annotation and DeepTMHMM complete for 100/100 proteins;
- 57 named AFDB hits and 32 novel / 68 non-novel calls;
- five proteins with 11 total predicted transmembrane regions;
- dual pocket status complete for F0-F5;
- 12/12 declared FoldTree outputs nonempty;
- atlas, status, summary, and config endpoints all returned HTTP 200;
- Tailscale access from the deployment client returned HTTP 200.

Immediate pre-deployment rollback material is at:

`/home/claude/suss_backups/20260714-014241-v1.0.1-predeploy`

An earlier untouched backup is also retained at
`/home/claude/suss_backups/20260714-pre-hardening`. To roll back, stop the PID in
`/home/claude/suss_portal/suss_portal.pid`, restore the three deployment files from the backup
(`suss_engine.tar.gz`, `suss_portal.py`, and `launch.sh`), then run `launch.sh`. Do not remove
the production `runs/` directory.
