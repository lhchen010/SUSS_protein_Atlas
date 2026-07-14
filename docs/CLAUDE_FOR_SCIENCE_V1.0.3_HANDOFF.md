# SUSS Protein Atlas v1.0.3: Claude for Science handoff

## Scope

This release extends the family-level `Download all Fx data` workbook. It does not change family
membership, clustering thresholds, structural alignment, annotation calls, or scientific metrics.

## Behaviour change

Every family workbook now contains an `annotation` sheet with one row per family member. The sheet
is copied from the run's `member_annotation.csv` without dropping or renaming columns, so it retains:

- accession and family identifiers;
- overall annotation status and component execution states;
- InterPro/Pfam domains and InterPro entries;
- PDB100 and AFDB-SwissProt Foldseek hits, scores, identities, and protein name;
- EffectorP call and status;
- DeepTMHMM transmembrane-region count and status;
- domain/fold evidence flags and the tri-state novelty call.

The workbook `README` sheet documents the new table. Existing sheets and download behaviour are
unchanged.

## Validation contract

Automated tests verify that:

1. the `annotation` sheet is present;
2. its accession rows match the requested family members;
3. biological result columns are preserved;
4. InterPro, both Foldseek searches, EffectorP, and DeepTMHMM status columns are preserved.

## 4070 staging regression

The atlas was rebuilt on 4070 from successful 100-protein production job
`20260714-092945-cor-2a801063`. Every embedded workbook was decoded and inspected with pandas and
openpyxl.

| Family | Members | Annotation rows | Missing workbook sheets |
|---|---:|---:|---:|
| F0 | 28 | 28 | 0 |
| F1 | 26 | 26 | 0 |
| F2 | 23 | 23 | 0 |
| F3 | 7 | 7 | 0 |
| F4 | 2 | 2 | 0 |
| F5 | 2 | 2 | 0 |

All six annotation sheets contained the complete 25-column production schema, including every
required annotation value and execution-status field. Existing ZIP structures, transforms,
pockets, FoldTree records, RNA-seq rows, and superposition records also passed the v1.0.2
regression contract.

## Claude acceptance checks

1. Open a family in the deployed atlas and select `Download all Fx data`.
2. Open the downloaded workbook and confirm an `annotation` sheet is present.
3. Confirm the sheet contains one row per member shown in the family card.
4. Confirm `annotation_status`, `interpro_status`, `foldseek_pdb_status`,
   `foldseek_afdb_status`, `effectorp_status`, and `deeptmhmm_status` are present.
5. Confirm domain, structural-hit, EffectorP, TMR, and novelty fields agree with the annotation panel.
6. Run `pytest -q` and confirm all focused tests pass.

Release version: `1.0.3`.
