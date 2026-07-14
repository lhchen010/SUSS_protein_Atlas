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

The 4070 production validation must additionally verify for every family that the number of
`annotation` rows equals the number of family members, and that the complete required annotation
column set is present.

## Claude acceptance checks

1. Open a family in the deployed atlas and select `Download all Fx data`.
2. Open the downloaded workbook and confirm an `annotation` sheet is present.
3. Confirm the sheet contains one row per member shown in the family card.
4. Confirm `annotation_status`, `interpro_status`, `foldseek_pdb_status`,
   `foldseek_afdb_status`, `effectorp_status`, and `deeptmhmm_status` are present.
5. Confirm domain, structural-hit, EffectorP, TMR, and novelty fields agree with the annotation panel.
6. Run `pytest -q` and confirm all focused tests pass.

Release version: `1.0.3`.
