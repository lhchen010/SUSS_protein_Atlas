# SUSS Protein Atlas v5.2.4 UI handoff

Date: 2026-08-14  
Scope: domain-unassigned workbench, responsive annotation, and main/detail view controls  
Scientific algorithms changed: no

## Summary

This patch fixes three atlas presentation defects without changing F-family, D-family, S-subgroup,
Foldseek, FoldMason, US-align, pocket, ESM, RNA-seq, or annotation calculations.

1. Domain `Unassigned / filtered` is now a first-class table with search, sorting, filters,
   pagination, rich evidence columns, and a domain-context detail workbench.
2. Selecting an unassigned protein no longer redirects to the full-length singleton/family page.
   An explicit button is still available when the reviewer intentionally wants the parent
   full-length context.
3. Per-member annotation preserves all six evidence fields as the detail panel changes size.
4. `Hide detail view` now works after divider dragging and on narrow layouts; the main view fills
   the released area and all active networks/viewers are resized.

## Root causes

- `#unassigned` was omitted from the fixed-height and scroll-container selectors. The renderer
  inserted every record into one table, so later records were clipped by `#primary`.
- The unassigned row handler called `openProtein()`, whose defined behavior is to switch to the
  full-length family or singleton axis.
- Entering the unassigned mode automatically collapsed the detail panel. This competed with the
  user's own panel state.
- A divider drag writes an inline `width` on `#side`; the old collapsed CSS could not override it.
- The annotation table used a fixed 780 px minimum width inside a 560 px panel.

## Data contract additions

Each `DOMAIN_UNASSIGNED` record now includes direct per-protein fields used elsewhere in the atlas:

- gene/protein name, EffectorP, DeepTMHMM, novelty state
- Pfam and InterPro
- PDB100 and AFDB/Swiss-Prot Foldseek hits and TM values
- annotation and Foldseek execution status
- accession-level P2Rank/fpocket residues and scores
- RNA-seq peak condition/value in addition to the complete expression dictionary

Family-reference pockets are never assigned to another family member. The unassigned structure
viewer reports only accession-level pocket evidence.

## Acceptance evidence

- Local automated tests: `74 passed`.
- Python compilation and whitespace checks: pass.
- Atlas rebuilt from run `20260811-160901-BC-036270db` on the 4070:
  - 124 full-length families
  - 327 full-length singletons
  - 269 domain-unassigned proteins
  - 45,165,487-byte backend atlas preview
- Domain-unassigned table renders 50 rows per page and reports `1-50 of 269`.
- Selecting `ATZ45198.1` remains on `Domain analysis > Unassigned / filtered` and exposes structure,
  domain evidence, RNA-seq, annotation, and sequence tabs.
- At a 560 px detail panel, annotation uses a six-column table with no hidden width.
- At a 420 px detail panel, every member becomes a labeled field block; all six fields remain.
- After dragging the detail panel to 360 px, `Hide detail view` reduces it to 0 px and expands the
  main view to 1280 px.
- At a 680 px browser viewport, hiding detail expands the main view to the full remaining height.

## Claude acceptance checklist

1. Open the deployed BC atlas and select `Domain analysis > Unassigned / filtered`.
2. Confirm the first page contains 50 records and use next/previous pagination.
3. Search for `ATZ45198.1`; select it and confirm the domain axis remains active.
4. Review all five detail tabs and confirm the structure is loaded from the portal artifact route.
5. Open full-length F23 annotation; resize the detail panel through wide and narrow states and
   confirm Protein, EffectorP, DeepTMHMM, Pfam, PDB, and AFDB fields all remain visible.
6. Drag the divider, then hide and show the detail view. Repeat below 700 px viewport width.
7. Confirm F/D memberships and all scientific matrices are unchanged from v5.2.3.

## Atlas download architecture note

The production portal uses a server-backed HTML shell. Large structures, workbooks, and ZIP files
are intentionally fetched only when selected, so downloading the HTML alone is not an offline
copy. For public web publication, use a static atlas bundle (`index.html`, manifest, compressed
data, and structure objects) on a CORS-capable static/object host. Google Drive is suitable for a
single archival ZIP but is not recommended as the live structure endpoint because direct links,
CORS behavior, throttling, and range requests are not stable application contracts.
