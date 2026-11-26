# SUSS Protein Atlas v2.0.0: Singleton Workbench

Date: 2026-07-26

## Purpose

Version 2 separates unclustered singleton proteins from the structural-family network. A
singleton is an independent protein that did not enter a Leiden family at the configured
structural threshold. It is not a one-member family and all singletons are not one shared
cluster.

The existing clustering result and SUSS family definitions are unchanged.

## User-visible behavior

The atlas now has two primary views:

1. **Cluster network** retains the comparative family workflow.
2. **Singletons** is a searchable, sortable, paginated table with direct protein evidence.

Selecting a singleton opens:

- its AlphaFold structure colored by pLDDT;
- fpocket and P2Rank predictions and residue downloads;
- ESM-1b tolerance, when enabled;
- RNA-seq expression, when supplied;
- InterPro/Pfam, EffectorP, and DeepTMHMM results;
- Foldseek PDB100 and AFDB/Swiss-Prot hits with TM scores;
- FASTA, PDB, pocket CSV, and evidence workbook downloads.

The singleton view intentionally omits FoldTree, FoldMason/Rate4Site conservation, Foldseek and
US-align within-family matrices, BLAST pair matrices, sequence identity matrices, structural
superposition, and the family network.

## Data contract

- `members.csv` remains the source of singleton membership (`family=singleton`).
- The stable singleton key is the protein accession itself.
- `pockets.json` stores singleton entries under the accession key.
- `esm_all.csv` stores singleton records with the accession in the `family` column for backward
  compatibility with the existing aggregate schema.
- `member_annotation.csv` remains the source for Foldseek database hits, domains, EffectorP,
  DeepTMHMM, component status, and novelty.
- `rnaseq_expression.csv` remains the all-protein expression source.
- The HTML payload contains separate `NET.nodes` and `SINGLETONS` arrays. Singleton accessions are
  never inserted as isolated network nodes.

## Scientific status semantics

Foldseek PDB100 and AFDB/Swiss-Prot searches remain meaningful for a singleton because they compare
the query protein with external structural databases. By contrast, family matrices and FoldTree
require multiple proteins from the submitted dataset and are therefore not calculated.

Novelty remains tri-state. A protein is only called novel when the required domain and structural
searches completed; incomplete evidence remains indeterminate.

## Downloads

The singleton Excel workbook contains direct annotation, Foldseek database hits and TM scores,
pocket summaries and predictions, pocket residues, detector-native pocket tables when available,
and RNA-seq expression. It intentionally excludes family-only sheets:

`foldseek_TM`, `usalign_TM`, `blast_identity`, `blast_pairs`, `foldtree`, `superposition`, and
`per_site`.

## Automated validation

The v2 test suite verifies:

- singleton workbooks include direct evidence and omit family-only analyses;
- PDB100 and AFDB/Swiss-Prot hit names, TM scores, and component states survive payload creation;
- indeterminate novelty remains null;
- a synthetic atlas embeds a singleton as an independent payload with an empty family network;
- the renderer exposes separate singleton search, detail, pLDDT, and download functions;
- all pre-existing family workbook, superposition, search, FoldTree, and runtime tests still pass.

Pre-deployment local result: `24 passed`.

## 4070 staging validation

The existing *C. orbiculare* run was upgraded in an isolated 4070 staging directory before
production deployment. The completed run contains:

- 1,144 proteins;
- 105 structural families;
- 249 singleton proteins;
- 354 pocket targets (105 family representatives plus 249 singletons);
- 354 complete fpocket results and 354 complete P2Rank executions;
- 354 ESM targets, including all 249 singletons;
- direct annotation, EffectorP, and DeepTMHMM records for all 249 singletons;
- 130 singleton PDB100 hits and 144 singleton AFDB/Swiss-Prot hits.

The complete staging workflow finished `112/112` Snakemake jobs. The original v2 HTML was
623 MB because singleton structures were redundantly embedded in several derived representations.
The final renderer stores each singleton PDB once and creates pLDDT, ESM, and pocket-colored
representations in the browser. The resulting self-contained atlas is 436,000,474 bytes.

Chromium acceptance testing confirmed:

- the page completes parsing with no console warnings or errors;
- the singleton button reports `Singletons (249)`;
- the table shows 50 rows per page and `1-50 of 249`;
- accession, annotation, EffectorP, DeepTMHMM, and Foldseek-prefixed searches update the table;
- `annotation:thioredoxin` returns 6 proteins;
- `effectorp:apoplastic` returns 56 proteins;
- `tmr:9` returns 4 proteins;
- `pdb:1xw9` returns 1 protein;
- TDZ13209.1 displays RNA-seq, PDB100/AFDB hits and TM scores, EffectorP, DeepTMHMM,
  fpocket, P2Rank, and ESM evidence;
- pLDDT, pocket, and ESM structure modes render a nonblank 3D canvas;
- a 390 x 844 viewport has no page-level horizontal overflow and stacks the detail panel below
  the singleton table.

## Claude acceptance checklist

1. Open the atlas and confirm only true families appear in **Cluster network**.
2. Switch to **Singletons** and confirm the count matches
   `members.csv[family == "singleton"]`.
3. Search by accession, annotation text, `effectorp:`, `tmr:`, `pdb:`, `afdb:`,
   `foldseek:`, `pocket:`, and `rnaseq:`.
4. Apply the Effector, Novel, Has pocket, and Has TM helix filters.
5. Sort by protein, annotation, Foldseek TM, pocket score, pLDDT, and RNA-seq peak.
6. Open a singleton and verify pLDDT, pocket, and ESM structure modes.
7. Verify the RNA-seq tab displays the selected protein only.
8. Verify annotation shows PDB100 and AFDB/Swiss-Prot hits and their TM scores.
9. Download FASTA, PDB, pocket CSV, filtered singleton CSV, and the singleton Excel workbook.
10. Confirm the singleton detail panel has no FoldTree, pairwise matrix, conservation, or
    superposition controls.
11. Open a normal family and repeat a family smoke test to confirm all v1 comparative functions
    remain available.

## Deployment record

Production deployment completed on 2026-07-26.

| Item | Production value |
|---|---|
| GitHub pull request | `#10` (`SUSS Protein Atlas v2.0.0: singleton workbench`) |
| Runtime merge commit | `43888c9cc69ef9553f31de1ef0509a78a2222e6c` |
| Release tag | `v2.0.0` |
| Engine package SHA-256 | `66a14ab922b51ceef503ad0dbca06f1e201049d6c7d5db459565d4c47244d0f7` |
| Production run | `/home/claude/suss_portal/runs/20260723-003323-cor-7a72a790` |
| Production atlas SHA-256 | `9d45e849a2d98f1f11d48bd97758138930a62721b1549e537b68e16441b88dac` |
| Production atlas size | `436,000,474` bytes |
| Pre-v2 backup | `/home/claude/suss_backups/20260726-2040-pre-v2.0.0` |
| Portal health | HTTP 200 on port 8600 |

The migrated production run reports 1,144 proteins, 105 families, and 249 singletons.
`used_config.yaml` records engine version `2.0.0` and the runtime merge commit above. The final
production workflow consistency pass completed `111/111` jobs after regenerating signatures,
cards, summaries, provenance, and the atlas from the migrated pocket and ESM aggregates.

Production browser validation used the portal's `/atlas` route, not the staging HTTP server. The
page completed parsing, displayed `Singletons (249)`, rendered 50 table rows with
`1-50 of 249`, hid the family network in singleton mode, and emitted no browser warnings or
errors.

To roll back, stop the portal, restore `suss_engine.tar.gz` and the files under `results/` from the
pre-v2 backup to the production paths, restore the backed-up `cards` directory, and run
`/home/claude/suss_portal/launch.sh`. The backup is additive and the original input files and
annotation cache were not deleted.
