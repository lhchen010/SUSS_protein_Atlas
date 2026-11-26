# SUSS Protein Atlas v1.0.2: Claude for Science handoff

## Scope

This release corrects structural superposition and expands family-level downloads. It does not
change Foldseek clustering, Leiden family membership, SUSS classification thresholds, or existing
annotation calls.

## Behaviour changes

### Structural viewer superposition

The old renderer loaded every selected PDB in its original coordinate frame. The structures were
therefore displayed together but were not superposed.

The builder now:

1. Uses the canonical family hub as the fixed reference.
2. Reads FoldMason's amino-acid structural MSA to identify corresponding CA atoms.
3. Computes a proper-rotation Kabsch fit for each member.
4. Iteratively excludes pairs farther than 4 A so flexible loops do not displace the conserved
   structural core.
5. Embeds only the compact rotation and translation matrices. The browser applies those transforms
   to both the 3Dmol view and the downloadable multi-model superposed PDB.

The Excel `superposition` sheet records the reference, correspondence method, retained and total CA
counts, core and all-pair RMSD, rotation matrix, and translation vector for every member.

### Member structure download

`Download members -> All structures` now returns `<family>_member_structures.zip`. The archive
contains one original `<accession>.pdb` per member plus `manifest.tsv`. It no longer combines all
members into one multi-model PDB.

### Complete family workbook

`Download all <family> data` now creates an auditable workbook with these sheets when applicable:

| Sheet | Contents |
|---|---|
| `README` | Sheet-level data dictionary |
| `members` | Family/member mapping |
| `foldseek_TM` | Symmetric Foldseek TM-score matrix |
| `usalign_TM` | Independent US-align TM-score matrix |
| `blast_identity` | BLASTp best-HSP identity matrix |
| `blast_pairs` | Pair-level BLAST fields and SUSS classification |
| `pocket_summary` | fpocket/P2Rank status and top-pocket summaries |
| `pocket_predictions` | One row per detected pocket |
| `pocket_residues` | One row per pocket-lining residue |
| `fpocket_pockets` | All descriptors from the native fpocket info file |
| `p2rank_pockets` | Every column from the native P2Rank predictions CSV |
| `foldtree` | Newick for every configured metric: foldtree, alntmscore, and lddt |
| `RNAseq` | Per-member, replicate-collapsed expression values |
| `per_site` | Conservation, SASA, pocket, and other residue evidence |
| `superposition` | Fit provenance and transforms |

For runs produced before v1.0.2, the builder backfills complete pocket tables from retained raw
fpocket and P2Rank outputs. New runs also write all pocket predictions into `pockets.json`.

## Automated validation

- Local: 12 tests passed under Python 3.11.
- GitHub Actions: `python-tests` passed for push and pull-request events.
- Synthetic rigid transform: recovered with RMSD below 1e-6 A.
- Synthetic FoldMason gap correspondence: correct atoms selected, RMSD below 1e-6 A.
- ZIP contract: one PDB per member plus manifest.
- Workbook contract: complete evidence-sheet set and pocket methods verified.
- Legacy pocket compatibility: full tables reconstructed from raw detector outputs.

## 4070 regression

Candidate atlas was rebuilt from production job `20260714-083204-cor-3113bbf5`:

| Family | Members / ZIP PDBs / transforms | Pocket predictions | FoldTree metrics | Max trimmed-core RMSD |
|---|---:|---:|---:|---:|
| F0 | 28 / 28 / 28 | 5 | 3 | 2.483 A |
| F1 | 26 / 26 / 26 | 7 | 3 | 2.221 A |
| F2 | 23 / 23 / 23 | 11 | 3 | 2.278 A |
| F3 | 7 / 7 / 7 | 6 | 3 | 2.492 A |
| F4 | 2 / 2 / 2 | 14 | 0 | 1.701 A |
| F5 | 2 / 2 / 2 | 10 | 0 | 1.050 A |

F4 and F5 correctly have no FoldTree output because FoldTree requires at least three structures.
F0 and F1 have complete P2Rank status but no P2Rank pocket predictions; their workbook still
contains the native empty P2Rank table and complete fpocket results.

The F0 28-member superposition was also inspected in the deployed 3Dmol renderer. The structural
core is co-located and variable loops remain visibly splayed, as expected.

## Claude acceptance checks

1. Open F0 and switch to `Superpose selected`; confirm the legend says hub-referenced
   FoldMason/Kabsch superposition and the core is aligned.
2. Toggle a tree tip or clade and confirm the viewer updates only the selected members.
3. Download `All structures (ZIP)` and confirm 28 PDB files plus `manifest.tsv` for F0.
4. Download the F0 workbook and inspect all sheets listed above.
5. Confirm F0 has 28 RNA-seq rows, both TM matrices, BLAST pair data, five fpocket pockets, and
   three FoldTree Newick records.
6. Confirm F2 includes native outputs for both fpocket and P2Rank.
7. Confirm F4/F5 have empty FoldTree sheets rather than fabricated trees.
8. Run `pytest -q`; all 12 tests should pass.

## Files changed

- `workflow/builders/html_builder.py`
- `workflow/builders/template/renderer.js`
- `workflow/scripts/sasa_pocket.py`
- `workflow/Snakefile`
- `tests/test_atlas_downloads.py`
- `.github/workflows/ci.yml`

Release version: `1.0.2`.
