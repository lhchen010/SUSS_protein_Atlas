# SUSS Protein Atlas v5.0.0 handoff

This document describes the v5.0.0 changes for Claude for Science and provides a
focused validation checklist. The release keeps the existing full-length family
workflow while adding a separate domain-family evidence workspace.

## Analysis scopes

- **Full**: full-length structural families (`F`), plus protein singletons.
- **Domain**: retained local structural-domain families (`D`), plus proteins with
  no retained D-family match.
- **Combined**: runs both scopes and builds one atlas linking `F`, `D`, and
  singleton records.

An `F` family and a `D` family are not interchangeable. An `F` family describes
whole-protein structural similarity. A `D` family describes a retained local
structural match and can connect proteins whose remaining regions differ.

## Domain-family evidence

Each D family now has the same evidence-oriented organization used by the full
view:

- Structure and local Foldseek network.
- Explicitly separated FoldMason structural guide tree, optional FoldTree, and
  MAFFT/FastTree sequence trees.
- Independent all-pairs US-align matrix and hub transforms.
- Sequence downloads, FoldMason AA/3Di alignments, and MAFFT alignments for
  reciprocal-coverage sequence subgroups.
- FoldMason residue-level structural conservation.
- RNA-seq, annotation, P2Rank, fpocket, and available ESM evidence.
- Links between parent F families and every matching D family.

The structure selector has one interaction model. Selecting one segment shows its
complete parent protein and highlights the matching domain. Selecting two or more
segments and pressing **Superpose selected** aligns the complete parent proteins
using the retained domain transforms. This makes the aligned region and the
unaligned parent context visible together.

## Trees and conservation

- FoldMason guide trees are labelled as structural guide trees.
- MAFFT/FastTree trees are labelled as sequence-relationship trees.
- FoldTree is never presented as a sequence tree. Domain FoldTree is optional and
  defaults to disabled.
- Structural conservation uses FoldMason lDDT. Higher conservation is red and
  lower conservation is blue.
- P2Rank and fpocket run on complete parent proteins, then residues are mapped to
  domain coordinates. Cropped domain structures are not used for pocket
  prediction.

## ESM scope

The default is `signals.esm_scope: representatives`, matching the representative
structure used by each F family plus all protein singletons. For the reference
dataset this is 354 proteins (105 F representatives and 249 singletons).

`signals.esm_scope: all_proteins` remains available for users who explicitly want
masked-marginal scans for every QC-passing protein and accept the additional
runtime. The atlas labels domain ESM evidence as available only when the complete
parent protein was scanned.

## Output organization

Large reference PDB payloads are served lazily by the portal rather than embedded
repeatedly in the HTML. The atlas remains a single entry page, while detailed
artifacts live under `results/downloads/`.

Each D family has:

- cropped-domain structure ZIP;
- complete-parent structure ZIP;
- domain and parent FASTA;
- FoldMason AA and 3Di MSA;
- sequence subgroup MAFFT MSA and tree;
- US-align matrix and transforms;
- local Foldseek edges;
- conservation tables;
- Excel workbook;
- complete evidence package ZIP.

The D-family package includes all available files above. Missing optional analyses
are represented by explicit status metadata, not empty success files.

Each D-family Excel workbook contains these sheets:
`members`, `foldseek_local_links`, `usalign_TM`, `sequence_MSA`,
`foldmason_AA`, `foldmason_3Di`, `sequence_subgroups`, `trees`,
`superposition`, `RNAseq_parent_proteins`, `pockets_parent_mapped`,
`annotation`, `structural_conservation`, and `README`.

## Portal controls

The portal exposes **Full**, **Domain**, and **Combined** analysis scopes. Advanced
options include domain FoldTree. Full-only sequence MSA, conservation, and
full-family FoldTree options are disabled when Domain-only is selected.

The result page has four primary views:

- Full-length families
- Domain families
- No D-family match
- Protein singletons

The global search continues to cover gene/accession, annotation, EffectorP, TM
results, database hits, and novelty labels. Search results highlight and open the
corresponding F family, D family, unmatched protein, or singleton.

## Reference validation on 4070

Reference run: `20260729-v5-reference`

- 1,144 QC-passing proteins.
- 105 full-length families and 249 protein singletons.
- 172 domain families containing 1,014 retained segments.
- 149 proteins with no retained D-family match.
- FoldMason completed for 172/172 D families.
- 842 non-hub domain transforms completed.
- Median domain-fit RMSD: 2.245 A.
- Median symmetric domain TM-score: 0.74525.
- No missing cells in retained D-family US-align matrices.
- P2Rank and fpocket completed for 1,144/1,144 proteins.
- P2Rank reported one or more pockets for 1,072 proteins; 72 valid runs reported
  no pocket.
- Default representative ESM evidence completed for 354/354 proteins.
- 172/172 D-family package ZIP files passed archive integrity checks.
- 172/172 D-family Excel files opened successfully and contained the required
  sheets.

## Claude acceptance checklist

1. Open `20260729-v5-reference` and confirm the four primary result views.
2. Open D1 and verify that one selected segment shows a complete parent protein
   with its matching region highlighted.
3. Select two or more D1 members and verify that **Superpose selected** aligns the
   matching domains while retaining full parent context.
4. Confirm that Trees presents FoldMason, FoldTree status, and MAFFT/FastTree as
   distinct analyses.
5. Confirm that Structural similarity contains both the independent US-align
   matrix and retained Foldseek edge table.
6. Confirm that RNA, Annotation, Sequence/MSA, and Conservation/Pockets display
   domain-member data rather than full-family summary rows.
7. Search `NLP` and confirm D35 (Necrosis inducing protein / NPP1) is found.
8. Open an F family from a D-family link and navigate back through the D-family
   link.
9. Open a protein singleton and verify structure, RNA-seq, annotation, pocket,
   sequence, PDB, and Excel downloads.
10. Download a D-family structure ZIP and complete package ZIP, then confirm that
    each contains multiple member files.
11. Confirm the red/blue structural-conservation legend reads high/red and
    low/blue.
12. Submit one Full, one Domain, and one Combined portal job and confirm disabled
    options and generated views match the selected scope.

## Known intentional behavior

- Domain FoldTree is disabled by default because FoldMason and independent
  US-align already supply structural-family evidence; enabling it adds cost and
  should be treated as an optional cross-check.
- A protein with no retained D-family match is not asserted to be domain-free.
- Representative-level ESM is the default performance policy, not missing data.
- D-family numbering is independent of F-family numbering. Cross-links, not
  matching numeric IDs, define their relationship.
