# SUSS Protein Atlas v3.1.0 handoff for Claude for Science

## Purpose

This release corrects the interpretation and display of structural-conservation scores and
turns Domain families from a flat list into an evidence-linked network workbench. It preserves
the v3 whole-fold (`F`), local-domain (`D`), sequence-subgroup (`S`), and singleton definitions.

## Scientific corrections

### Structural conservation

The previous renderer reversed the intended red/blue interpretation. The previous per-site
calculation could also give a high hub-relative score to a sparsely supported alignment column.

v3.1 now uses FoldMason's official `msa2lddtjson` output:

```text
FoldMason structural MSA
    -> msa2lddtjson --pair-threshold 0.5
    -> map alignment columns to reference residues
    -> leave unsupported columns missing
    -> write CSV plus B-factor-encoded PDB
```

Interpretation in the atlas:

- red: high structural conservation (high FoldMason LDDT);
- white: intermediate structural conservation;
- blue: low structural conservation;
- light gray: no score at the configured pair-support threshold.

The default pair threshold is configurable:

```yaml
structural_conservation:
  pair_threshold: 0.5
```

The pLDDT and ESM gradients were audited at the same time:

- pLDDT: red is higher confidence and blue is lower confidence;
- ESM tolerance: red is mutation-tolerant and blue is constrained.

### Domain-family evidence model

The local-domain clustering calculation is unchanged from v3.0. D families remain connected
components/Leiden communities of consolidated local 3Di+AA Foldseek segments. v3.1 exposes the
underlying evidence rather than presenting only a table.

The atlas now has two domain-network levels:

1. The overview network has one node per D family. An edge means that retained local Foldseek
   relationships bridge segments assigned to the two D communities.
2. Selecting a D family opens a segment network. Nodes are protein segments and edges are the
   retained local Foldseek hits used by the domain analysis.

Each segment detail view provides:

- the full parent structure, with the matched segment in teal and available pocket residues in
  red;
- parent F-family or singleton assignment;
- local Foldseek probability, LDDT, aligned length, identity, and coverage in edge tooltips;
- coordinate-overlapping Pfam/InterPro annotations;
- Foldseek PDB100 and AFDB/Swiss-Prot annotation evidence;
- EffectorP and DeepTMHMM results;
- RNA-seq values;
- segment CSV, local-edge CSV, segment FASTA, and selected-segment PDB downloads.

Search now matches D-family ID, accession, segment coordinates, gene/annotation text,
EffectorP/DeepTMHMM evidence, parent family, and Foldseek-hit labels. A unique match can be opened
directly with Enter.

## Data-contract changes

`domain_cluster` now produces two additional declared outputs:

```text
results/domain_edges.csv
results/domain_cross_edges.csv
```

`domain_edges.csv` contains retained segment-level relationships and their Foldseek statistics.
`domain_cross_edges.csv` aggregates bridges between D communities for the overview network.
Existing `domain_families.csv` and `domain_members.csv` remain available.

The atlas payload adds:

```text
DOMAIN_EDGES
DNET.nodes
DNET.edges
```

Domain members are enriched at build time from existing annotation, pocket, RNA-seq, F/singleton,
sequence, and embedded-structure payloads. No new annotation method is inferred from the network.

## Configuration

The v3 Foldseek defaults remain:

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

structural_conservation:
  pair_threshold: 0.5
```

`docs/FOLDSEEK_PARAMETERS.md` explains why whole-fold TM thresholds and local-domain
3Di+AA thresholds answer different questions.

## Acceptance checklist for Claude

1. Confirm `structural_conservation.py` calls FoldMason `msa2lddtjson` and converts `-1` or absent
   scores to missing values.
2. Confirm only `structural_scored_resi` receives the red-white-blue color scheme; unsupported
   residues remain light gray.
3. Confirm structural conservation is red for high LDDT and blue for low LDDT.
4. Confirm pLDDT is red for high confidence and ESM is red for mutation tolerance.
5. Open Domain families and confirm the first view is a D-family overview network.
6. Search `TDZ22527.1`, open D20, and confirm the segment network contains
   `TDZ22527.1:2-127` and `TDZ21966.1:1-118`.
7. Select the direct edge between those two segments and confirm it reports approximately:
   probability 1.0, local LDDT 0.7476, alignment length 132, and shorter coverage 1.0.
8. Select a segment and confirm the 3D viewer shows the complete parent protein, highlights only
   the matched interval in teal, and displays available pocket residues in red.
9. Confirm overlapping domain annotations, parent F/singleton assignment, database hits,
   EffectorP/DeepTMHMM, and RNA-seq are visible without switching to the full-length mode.
10. Exercise the domain CSV, local-edge CSV, segment FASTA, and segment-PDB downloads.
11. Confirm Full-length families and Singletons still load and their existing downloads work.
12. Run the Python test suite and a Snakemake dry run.

## Validation record

Validation date: 2026-07-28 (Asia/Taipei).

- Local Python tests: 35 passed.
- Python compile, JavaScript syntax, and `git diff --check`: passed.
- GitHub Actions branch and pull-request `python-tests`: passed.
- Domain regression: 106 D families, 1,099 consolidated segments, 1,080 proteins, and
  49 cross-family overview bridges.
- D20 acceptance example: 13 segments and 31 retained local edges. The two requested proteins
  have a direct edge with probability 1.0, local LDDT 0.7476, 132 aligned residues, and
  shorter-segment coverage 1.0.
- F4 targeted conservation smoke test: 194 of 231 reference residues received official
  FoldMason LDDT scores; all scored columns met the 0.5 pair-support threshold and the PDB
  B-factor round trip differed by at most 0.005.
- ESM reuse was checksum-equivalent rather than existence-based: all 1,144 current sequences
  exactly matched the validated source run, and the cache covered all 354 requested references.
- Full 4070 rebuild: 535/535 Snakemake steps completed. A follow-up dry run reported that all
  requested files were present and up to date.
- Whole-fold regression: 1,144 proteins, 105 F families, 895 clustered proteins, and
  249 singletons.
- Domain outputs: 8,807 retained segment edges in addition to the 49 aggregated D-family
  bridges.
- Structural conservation: 105/105 F families produced CSV and B-factor-encoded PDB outputs.
- Annotation: 1,144/1,144 proteins reported `annotation_status=complete`.
- Pocket analysis: 354/354 references reported complete P2Rank results and 354/354 reported
  complete fpocket results. SASA covered 456,195 residues from all 1,144 proteins.
- Atlas output: 474,392,011 bytes (452.42 MiB), SHA-256
  `025af08977a4f3c92f061a5c3d125df78c086e2eba9af7ca390e84a23acc5ba7`.
- Remote Python compile and JavaScript syntax checks passed. The production Python environment
  intentionally does not include pytest; the 35-test suite passed locally and in GitHub Actions.
- Browser acceptance on the production portal passed:
  - the 106-node D-family overview and D20's 13-segment network rendered at non-zero canvas sizes;
  - the D20 workbench exposed structure, local edges, parent placement, annotations, pockets,
    EffectorP/DeepTMHMM, RNA-seq, and downloads;
  - F4 displayed the corrected `structurally variable` blue,
    `structurally conserved` red, and `insufficient pair coverage` gray legend, with the official
    FoldMason LDDT label and pair threshold 0.50;
  - the `TDZ22527.1` singleton retained its structure, pLDDT/pocket/ESM modes, annotation,
    RNA-seq, FASTA/PDB downloads, and complete evidence workbook.

## Deployment and rollback

- Production portal: `http://100.80.77.29:8600/`
- Validated reference atlas:
  `http://100.80.77.29:8600/atlas?id=20260728-v3.1-reference`
- Staging engine: `/home/claude/suss_atlas_staging/20260728-v3/engine`
- Staging backup: `/home/claude/suss_atlas_staging/20260728-v3/backups/20260728-030852`
- Production engine backup:
  `/home/claude/suss_portal/suss_engine.tar.gz.pre-20260728-v3.1.0`
- GitHub pull request: `https://github.com/lhchen010/SUSS_protein_Atlas/pull/16`
- GitHub release:
  `https://github.com/lhchen010/SUSS_protein_Atlas/releases/tag/v3.1.0`

For rollback, restore the backed-up `suss_engine.tar.gz`, restart the portal only if the portal
process itself was changed, and select the previous validated reference run from portal history.
