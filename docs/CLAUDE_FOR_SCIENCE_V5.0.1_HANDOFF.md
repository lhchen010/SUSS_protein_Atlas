# SUSS Protein Atlas v5.0.1 handoff

This maintenance release restores Full-family evidence that was absent from the
v5.0.0 reference build and completes visual and analytical parity between the
Full and Domain workbenches. It is intended for review and regression testing by
Claude for Science.

## Root cause and workflow correction

The `cluster` checkpoint previously declared `results/families/` as a directory
output. Snakemake owns directory outputs and may replace the complete directory
when the checkpoint is rerun. Because downstream FoldMason, FoldTree, BLAST,
US-align, and conservation files were also stored under that directory, a
re-clustering run could remove completed Full-family analyses while preserving a
valid family summary. The resulting atlas therefore looked complete but lacked
Full-family heatmaps, FoldTree figures, and conservation assets.

v5.0.1 separates ownership:

- `results/family_members/` contains checkpoint-generated membership files.
- `results/families/` contains downstream family analyses and downloads.
- per-family rules read membership from the first directory and write analyses
  only to the second.

This is a workflow durability fix. It does not change the reference F-family
membership or scientific clustering thresholds.

## Full-family restoration

The v5.0.1 reference contains the complete Full-family assets:

- Foldseek TM-score heatmap;
- independent US-align TM-score heatmap;
- BLASTp best-HSP sequence-identity heatmap;
- FoldTree structural relationship figure and interactive tree;
- FoldMason structural conservation and eligible Rate4Site sequence
  conservation;
- MAFFT sequence MSA and FoldMason AA/3Di structural MSA.

Missing optional files now produce an explicit unavailable message instead of a
broken image.

The Full structure viewer also adds white and black backgrounds while retaining
cartoon, surface, stick, sphere, and line representations.

## Domain-family parity

Each D family now exposes the same evidence organization where the analysis is
scientifically applicable:

- white/black viewer backgrounds;
- cartoon, surface, stick, sphere, and line representations;
- full-parent fpocket and P2Rank controls using their reported lining residues;
- FoldMason-aligned domain sequence-identity heatmap;
- FoldMason structural conservation across the complete D family;
- Rate4Site sequence conservation only within eligible reciprocal-coverage
  BLAST/MAFFT subgroups;
- per-segment sequence-conservation status and residue scores;
- workbook and package exports for sequence identity and sequence conservation.

Domain sequence identity is the pairwise amino-acid identity over FoldMason AA
alignment columns where both segments contain a residue. It describes sequence
identity within the structurally aligned domain region; it is not a whole-protein
BLAST identity.

Rate4Site uses MAFFT alignments from BLAST-supported sequence subgroups with at
least the configured minimum number of sequences. A D family may therefore have
complete structural conservation but unavailable sequence conservation. This is
intentional and prevents unrelated proteins from being interpreted as a
sequence-homologous group.

## Reference validation

Reference run: `20260730-v5.0.1-reference`

- 1,144 QC-passing proteins.
- 105 Full families and 249 protein singletons.
- F-family membership is identical to the previous validated reference: 1,144
  accession-to-family assignments compared equal.
- 105/105 Full-family directories contain Foldseek TM, BLAST identity,
  US-align TM, and structural-conservation outputs.
- Eligible Full families contain FoldTree Newick output.
- 172/172 D families contain a sequence-identity matrix.
- Every D sequence-identity matrix is symmetric with a diagonal of 1.
- 69 D sequence subgroups completed Rate4Site.
- 545 D segments contain sequence-conservation mappings.
- 187,333 D residue-level sequence-conservation scores were exported.
- 172/172 D workbooks opened successfully with all 16 required sheets.
- 172/172 D complete-package ZIP files passed archive checks.

The 16 D workbook sheets are:

`members`, `foldseek_local_links`, `usalign_TM`, `sequence_identity`,
`sequence_MSA`, `foldmason_AA`, `foldmason_3Di`, `sequence_subgroups`,
`trees`, `superposition`, `RNAseq_parent_proteins`,
`pockets_parent_mapped`, `annotation`, `structural_conservation`,
`sequence_conservation`, and `README`.

## Browser acceptance checklist

1. Open `20260730-v5.0.1-reference`, then open D1.
2. In **Structure + Network**, switch White/Black and each structure style.
3. Select **Pocket**, switch between fpocket and P2Rank, and confirm the complete
   parent protein remains visible with pocket-lining residues highlighted.
4. Select **Sequence conservation** and confirm eligible residues are colored
   without replacing structural-conservation semantics.
5. Open **Sequence + MSA** and confirm the domain sequence-identity heatmap,
   MAFFT MSA, FoldMason AA MSA, and FoldMason 3Di MSA are distinct.
6. Open **Conservation + Pockets** and confirm structural-conservation,
   Rate4Site, ESM, fpocket, and P2Rank counts are reported separately.
7. Download a D workbook and confirm `sequence_identity` and
   `sequence_conservation` sheets.
8. Download a D complete-package ZIP and confirm
   `tables/domain_sequence_identity.csv` and
   `tables/rate4site_sequence_conservation.csv`.
9. Open F43 and confirm White/Black backgrounds.
10. Open **FoldTree (figure)** and confirm the structural and sequence trees are
    both rendered and correctly labelled.
11. Open **Struct sim (TM)** and confirm both Foldseek and US-align heatmaps.
12. Open **Sequence + MSA** and confirm the BLASTp sequence-identity heatmap and
    the three alignment views.

## Automated validation

- Python source compilation: passed.
- Renderer JavaScript syntax check: passed.
- Test suite: 48 passed.
- Artifact audit: passed for all 105 F families and 172 D families.
- Staging browser acceptance: passed for D1 and F43.

## Intentional limitations

- Sequence conservation is unavailable when no eligible BLAST-supported subgroup
  reaches the configured minimum size.
- Domain pocket evidence is computed on the complete parent protein and mapped
  to domain coordinates; cropped segments are not submitted to pocket detectors.
- Domain sequence identity and Full BLAST identity answer different questions
  and must not be compared as the same metric.
- FoldTree is a structural relationship tree, not a sequence phylogeny.
