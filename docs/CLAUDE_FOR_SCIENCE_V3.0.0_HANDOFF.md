# SUSS Protein Atlas v3.0.0 handoff for Claude for Science

## Purpose

This document records the v3 scientific redesign, implementation, validation, and deployment
state so Claude for Science can independently inspect and test the release.

## Scientific model

- `F` families are whole-protein fold communities from global Foldseek TM similarity,
  reciprocal structural coverage, and Leiden clustering.
- `D` families are local structural-segment communities from Foldseek 3Di+AA hits.
- `S` subgroups are reciprocal-coverage-controlled BLAST connected components inside F families.
- FoldMason AA/3Di alignments cover the structural family and support superposition and
  structural conservation.
- MAFFT sequence alignments cover eligible S subgroups.
- Rate4Site is run only on a sufficiently large homologous S subgroup; unavailable results remain
  null and are not presented as evolutionary conservation.
- FoldTree is labeled as a structural relationship tree. FastTree on the MAFFT alignment is
  labeled as an exploratory sequence relationship tree.

## Major implementation changes

- Added configurable reciprocal structural coverage to global F-family edges.
- Added local segment extraction, overlap consolidation, Leiden D-family clustering, and CSV
  outputs.
- Added BLAST reciprocal coverage, S-subgroup tables, MAFFT L-INS-i, FastTree WAG, and typed
  sequence-analysis status.
- Added FoldMason-derived per-site occupancy, amino-acid entropy, 3Di entropy, and hub-relative
  local-distance agreement.
- Added an atlas Foldseek database and protein-to-F/D/S/singleton index.
- Added a portal PDB/mmCIF structure-search endpoint.
- Split atlas sequence and structural MSA viewers and downloads.
- Split structural conservation from sequence evolutionary conservation.
- Added D-family atlas mode and links from segments to F-family or singleton evidence.
- Expanded per-family Excel exports with annotation, pockets, Foldseek, US-align, BLAST,
  S subgroups, D memberships, both MSA types, trees, structural conservation, Rate4Site status,
  and RNA-seq.

## Configuration defaults

See `docs/FOLDSEEK_PARAMETERS.md`. Important defaults:

```yaml
clustering:
  alignment_type: 1
  foldseek_tm: 0.5
  whole_fold_min_coverage: 0.5
  tm_symmetric: min

domain_clustering:
  alignment_type: 2
  evalue: 0.001
  min_probability: 0.5
  min_aligned_residues: 40
  min_shorter_coverage: 0.0

classification:
  blast_evalue: 0.001
  min_reciprocal_coverage: 0.5
```

## Acceptance checklist for Claude

1. Confirm an F edge requires both the configured symmetric TM and reciprocal coverage.
2. Confirm D families use segment nodes and can link proteins from different F families.
3. Confirm a protein may have multiple D-family memberships.
4. Confirm S subgroups never redefine F-family membership.
5. Confirm MAFFT input is restricted to the representative protein's S subgroup.
6. Confirm Rate4Site is not run on unrelated members merely because they share an F family.
7. Confirm structural conservation uses FoldMason AA/3Di correspondence and is labeled separately.
8. Confirm FoldTree and sequence-tree labels do not claim the same evolutionary interpretation.
9. Upload a query structure in the portal and verify hits map to F/D/S/singleton records.
10. Open the atlas and test all three modes, structure superposition, MSA viewers, ZIP structure
    download, and the complete family Excel workbook.

## Validation record

Validation date: 2026-07-28 (Asia/Taipei).

- Local unit/integration tests: 31 passed on Python 3.11.
- Python compile check and `git diff --check`: passed.
- JavaScript syntax check with Node.js on 4070: passed for `databridge.js` and `renderer.js`.
- Snakemake dry run on 4070: parsed the complete workflow successfully.
- Full 4070 staging regression: 320/320 steps completed in about eight minutes with 16 cores.
- Whole-fold regression: all 1,144 proteins retained exactly the same F-family/singleton
  membership as v2.1 (105 F families, 895 clustered proteins, 249 singletons).
- Local-domain result: 106 D families from 1,099 consolidated segments in 1,080 proteins;
  13 proteins have multiple D memberships and 45 D families connect more than one F/singleton
  group.
- Sequence result: 302 S subgroups. MAFFT, FastTree, and Rate4Site completed for 38 eligible
  representative subgroups; 67 families were explicitly `not_applicable`.
- Structural conservation: 105/105 F families produced both CSV and conservation-colored PDB.
- Search database: 1,144/1,144 proteins indexed.
- Domain-aware acceptance example:
  `TDZ22527.1` (full-length singleton, residues 2-127) and `TDZ21966.1`
  (`F23`, residues 1-118) both map to `D20`.
- Portal structure-search smoke test: a `TDZ22527.1` PDB query returned HTTP 200 and included
  `TDZ22527.1`, `TDZ21966.1`, `D20`, `F23`, and singleton mappings.
- Browser acceptance test on the deployed reference atlas:
  - all three atlas modes loaded;
  - search reduced the D-family table to D20 and opened its 13 matched segments;
  - F4 loaded the structure viewer, FoldMason AA/3Di viewers, MAFFT viewer, sequence tree,
    conservation controls, superposition, ZIP, FASTA, MSA, and Excel controls;
  - `TDZ22527.1` loaded as a singleton with structure, pocket, RNA-seq, annotation, sequence,
    and evidence-download controls.
- Atlas build: completed without fallback; generated the atlas HTML, family summary workbook,
  master table, and composition workbook.
- GitHub Actions: pending release push.
- Production deployment: complete at `http://100.80.77.29:8600/`.
  The validated result is available in portal history as `20260728-v3-reference`.
  The previous engine and portal files were retained with `.pre-20260728-v3.0.0` suffixes.
