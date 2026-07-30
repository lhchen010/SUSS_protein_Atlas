# SUSS Protein Atlas v5.1.0 handoff

This release corrects the scientific ownership of Domain-family sequence
analysis. D families remain structurally defined by local Foldseek 3Di+AA
links, but their sequence subgroups, sequence trees, and Rate4Site analyses now
use independently extracted domain segments instead of whole-parent BLAST
relationships.

## Why this correction was required

The v5.0.2 Domain workbench cropped sequences before MAFFT, but decided which
segments belonged to a sequence subgroup using `blastp_allvsall.tsv`, whose
queries and coverage values described complete parent proteins. It also joined
two D segments automatically when they came from the same parent accession.

Those behaviors could:

- connect D segments because another part of the parent proteins was homologous;
- miss a homologous D segment because its parent-level reciprocal coverage was
  too low;
- merge distinct repeated domains from one parent without segment evidence;
- make a Domain sequence tree or Rate4Site status depend on Full-analysis data.

v5.1.0 removes all four paths.

## Independent D-family sequence workflow

The new declared Snakemake data flow is:

1. `domain_sequences` extracts every `accession:start-end` sequence and writes
   `results/domain_segments.fasta`.
2. It also writes `results/domain_sequence_manifest.csv`, preserving the D
   family, parent accession, coordinates, segment length, and parent length.
3. `domain_blastp` searches that segment FASTA against itself.
4. D sequence relationships require the configured BLAST e-value and reciprocal
   coverage on the cropped segments.
5. Connected components become D-specific sequence subgroups.
6. MAFFT, FastTree, and Rate4Site run independently inside eligible D
   subgroups.
7. Distinct segments from one parent are not joined unless their segment
   sequences have a qualifying hit.

BLAST coverage is clamped to the valid 0-1 range. Gapped alignments can otherwise
report an alignment length slightly greater than the ungapped query length.

## Structure and sequence evidence remain separate

Each D family now exposes two explicitly different amino-acid identity views:

- Domain-segment BLASTp best-HSP identity: independent sequence-search evidence.
- FoldMason-aligned amino-acid identity: identity over structurally
  corresponding columns.

The following D-family analyses use cropped D structures and remain independent
of F families:

- local Foldseek network;
- all-pairs US-align TM matrix and superposition transforms;
- FoldMason AA/3Di structural MSA;
- FoldMason per-column lDDT structural conservation;
- FoldMason structural guide tree;
- optional FoldTree structural trees.

## Parent-level evidence

Some measurements belong to a protein or gene, not to a family. Recomputing the
same value when a protein appears in both an F and D view would create duplicate
evidence rather than an independent family analysis.

The following are therefore calculated once and mapped by accession or residue
coordinates:

- RNA-seq expression;
- EffectorP and DeepTMHMM;
- parent annotation, with Pfam/InterPro calls filtered by coordinate overlap in
  the D view;
- fpocket and P2Rank predictions on the complete parent, with D-overlapping
  lining residues identified;
- ESM-1b parent-context residue scores.

ESM uses a symmetric representative policy by default:

- F workbench: the independently selected F representative;
- D workbench: the independently selected D structural hub;
- Full singleton: that singleton protein.

The D viewer opens on its D hub and labels members with available ESM data.
Selecting an unscanned non-hub reports that state explicitly instead of showing
a red/blue legend over a grey structure. Exhaustive `domain_members` and
`all_proteins` ESM scopes remain available at substantially higher cost.

## Reference regression

Deployed reference: `20260730-v5.1.0-reference`

Portal URL:
`http://100.80.77.29:8600/atlas?id=20260730-v5.1.0-reference`

Domain workbench schema: `4`.

Inputs and structural clustering are unchanged from v5.0.2:

| Artifact | v5.0.2 vs v5.1.0 |
|---|---|
| `members.csv` | identical SHA256 |
| `domain_members.csv` | identical SHA256 |
| `domain_edges.csv` | identical SHA256 |
| `foldseek_allvsall.tsv` | identical SHA256 |
| F families | 105, unchanged |
| D families | 172, unchanged |
| D segments | 1,014, unchanged |

D-specific sequence-analysis changes:

| Metric | v5.0.2 | v5.1.0 |
|---|---:|---:|
| D Rate4Site-eligible subgroups | 69 | 73 |
| D segments with Rate4Site scores | 545 | 596 |
| D families with complete FoldTree outputs | not enabled | 110 |

Focused examples:

- D40 remains one six-segment sequence subgroup and has D-specific MAFFT,
  FastTree, Rate4Site, FoldMason, US-align, and three FoldTree metrics.
- D70 changes from parent-BLAST `2 + 1 + 1` to domain-BLAST `3 + 1`.
  Its three-member subgroup remains below the default four-sequence Rate4Site
  minimum and is reported as `3 / 4`, not as a missing renderer.
- D94 remains one three-member domain-sequence subgroup and is explicitly
  reported as `3 / 4` for Rate4Site.

## Downloads and interface

Every D-family Excel workbook and complete ZIP now contains:

- cropped domain sequences and complete parent sequences;
- domain BLASTp hits;
- domain BLASTp identity matrix;
- separate FoldMason-aligned AA identity matrix;
- every D sequence-subgroup MAFFT alignment and FastTree;
- FoldMason AA/3Di MSA;
- US-align matrix and transforms;
- FoldTree and FoldMason structural trees;
- structural and sequence conservation tables;
- mapped pockets, ESM, RNA-seq, and annotation.

The Conservation panel always reports each D subgroup's segment count,
relationship count, Rate4Site minimum, and exact status. Unavailable evidence is
shown as unavailable rather than hidden.

## Automated validation

- Local Python tests: 55 passed.
- 4070 Python tests: 55 passed.
- GitHub Actions: passed.
- Domain workflow dry-run: passed.
- D workbench rebuild: 172/172 families completed.
- ESM representative scan: 453/453 proteins completed, including 172/172
  independent D-family hubs.
- D40/D70/D94 focused regression: passed.
- Download regression: D40, D70, and D94 Excel and complete ZIP packages passed.
- Portal artifact smoke test: atlas, D40/D70 Excel, and D40 complete ZIP returned
  HTTP 200.
- Browser regression: Domain overview and segment networks rendered without
  console errors; D40 structure/ESM/heatmaps, D70 and D94 Rate4Site states, and
  the F8 RNA-seq/FoldTree/sequence/structure views passed.

## Claude acceptance checklist

1. Open the v5.1.0 reference and confirm that Domain and Full analysis remain
   separate navigation axes.
2. Open D40, then Sequence + MSA. Confirm that the first heatmap is labeled
   Domain-segment BLASTp identity and the second is FoldMason-aligned identity.
3. Confirm that D40 exposes one six-segment MAFFT subgroup, a sequence tree,
   Rate4Site scores, FoldMason evidence, US-align, and FoldTree structural trees.
4. Open D70 and confirm two sequence subgroups with sizes 3 and 1.
5. Open D70 Conservation + Pockets and confirm Rate4Site reports `3 / 4`
   rather than silently disappearing.
6. Open D94 and confirm one three-member subgroup with Rate4Site status `3 / 4`.
7. In each D family, confirm the default selected member is the D hub.
8. Select ESM on the D hub and confirm parent-context residue coloring appears.
9. Select a non-hub without ESM and confirm the viewer explicitly reports that
   it was not scanned.
10. Download a D-family Excel and ZIP. Confirm both identity matrices,
    `domain_blastp_hits`, sequence-subgroup MAFFT/tree files, structures,
    conservation, pockets, annotation, and RNA-seq are present.
11. Run a new portal job and confirm `config/config.yaml` is initialized even
    when the engine archive contains only the 4070 reference config.

## Intentional semantics

- F and D memberships are independent and linked only through a crosswalk.
- F-family sequence evolution uses full-protein sequence evidence.
- D-family sequence evolution uses cropped D-segment sequence evidence.
- Structural similarity does not imply sequence homology.
- Rate4Site is not forced across sequence-unrelated structural members.
- RNA-seq and protein predictions are parent-level measurements mapped into
  the selected family view.
