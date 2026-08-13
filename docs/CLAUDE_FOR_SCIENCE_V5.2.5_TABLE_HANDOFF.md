# SUSS Protein Atlas v5.2.5 table handoff

## Scope

This release changes only the atlas presentation layer for two protein lists:

- Full-length analysis > Singletons
- Domain analysis > Unassigned / filtered

No clustering, annotation, pocket, expression, or structure-analysis result is recalculated.

## Changes

1. `EffectorP` and `TMR` are independent columns and sort keys.
2. `pLDDT` and protein `Length` are independent columns and sort keys.
3. Every rendered cell has a stable `data-label` used by the narrow-panel layout.
4. At a table-container width of 900 px or less, rows become complete two-column records.
5. At 520 px or less, rows become one-column records.
6. Annotation and rejection-reason fields span the record width for readability.
7. Narrow layouts do not hide or ellipsize evidence cells and do not require horizontal scrolling.
8. Wide layouts retain the compact sortable table presentation.

## Independent acceptance checks

1. Open Full-length analysis > Singletons and resize the main panel below 900 px.
2. Confirm one record exposes 10 labeled fields: Protein, Annotation, PDB100 hit, AFDB / Swiss-Prot,
   EffectorP, TMR, Pocket, pLDDT, Length, and RNA-seq peak.
3. Open Domain analysis > Unassigned / filtered at the same width.
4. Confirm one record exposes 14 labeled fields, including the same four split evidence fields and
   Domain evidence, Best local match, and Why no D family.
5. Confirm neither list has inaccessible cells or horizontal overflow in narrow mode.
6. Expand the main panel above 900 px and confirm column headers and sorting remain available.

## Automated coverage

`tests/test_atlas_downloads.py::test_protein_list_tables_split_evidence_and_preserve_narrow_fields`
guards the split columns, responsive breakpoints, per-cell labels, and visible overflow behavior.
