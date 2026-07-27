# SUSS Protein Atlas v2.1.0: Sequence and MSA Viewer

Date: 2026-07-27

## Scope

Version 2.1 adds sequence inspection and download controls without changing clustering,
annotation, structure analysis, or SUSS scientific definitions.

## Cluster behavior

The existing **Seq identity** tab is renamed **Sequence + MSA** and now contains:

- the existing BLASTp identity matrix;
- a compact member-sequence viewer;
- a FoldMason structure-guided MSA viewer;
- all mature member sequences as unaligned FASTA;
- the FoldMason MSA as aligned FASTA;
- the existing per-member FASTA download in the Structure tab.

The MSA payload comes from `results/families/<family>/<family>.aln`, which is already used to
establish residue correspondence for structural superposition. Alignment gaps and member order
are preserved. The filename is `<family>_FoldMason_MSA_<n>seqs.fasta`.

## Singleton behavior

A singleton receives a dedicated **Sequence** tab containing its mature amino-acid sequence,
length, and FASTA download. It does not receive an MSA button or synthetic one-sequence alignment.

## UI and payload impact

The viewer uses a fixed-height monospace panel with horizontal and vertical scrolling. It is
created only when the Sequence tab is opened, so it does not add hundreds of sequence elements to
the initial DOM. Only aligned sequence strings are added to the data payload; PDB, Excel, matrix,
and image data are not duplicated.

## Automated validation

The test suite verifies:

- FoldMason headers map back to accessions;
- alignment gaps are preserved;
- family MSA and sequence download functions exist;
- the compact viewer and CSS are present;
- singleton payloads explicitly contain no MSA;
- singleton UI describes MSA as not applicable.

Local result: `26 passed`.

## 4070 staging validation

- Full atlas rebuild completed successfully.
- Atlas size: `436,784,865` bytes, an increase of `784,391` bytes (`~0.18%`)
  compared with v2.0.
- Validated `105` FoldMason MSA files containing `895` aligned records.
- Every cluster MSA contains at least two records and has one consistent aligned length.
- F0 contains `62` aligned records with `1,582` alignment columns.
- Desktop and `390 px` mobile browser checks passed with no horizontal page overflow.
- No browser console errors were recorded.

## Production deployment

The merged runtime commit is `de6b6d1bebf805c315a7a1e888be66c4ca5dc136`.

- Portal: `http://100.80.77.29:8600/`
- Existing production result:
  `http://100.80.77.29:8600/atlas?id=20260723-003323-cor-7a72a790`
- Engine package SHA-256:
  `090373afc1323ffddb23b22dafecb8c5714067787a9a0df2d99fbec119c2a58a`
- Production atlas SHA-256:
  `fdfa9ff85477b0a51747273eb3d7b6a92a89e288aa070189e6524c20641f95fd`
- Production provenance reports `engine_version: 2.1.0` and the runtime commit above.
- Portal restart and HTTP health checks passed.
- Production browser smoke tests passed for F0 member/MSA views and singleton
  `TDZ13209.1`.

The pre-v2.1 rollback bundle is stored at:

`/home/claude/suss_backups/20260727-pre-v2.1.0`

## Claude acceptance checklist

1. Open a multi-member cluster and select **Sequence + MSA**.
2. Confirm the member view shows the selected accession and mature sequence length.
3. Change the member selector and confirm the displayed sequence changes.
4. Select **MSA** and confirm at least two aligned rows and gap characters are visible.
5. Download **All sequences (FASTA)** and confirm sequences are unaligned.
6. Download **MSA (FASTA)** and confirm equal aligned lengths and preserved gaps.
7. Open a singleton and select **Sequence**.
8. Confirm its sequence and FASTA download are available.
9. Confirm no MSA control is shown for the singleton.
10. Repeat a structure, RNA-seq, and annotation smoke test to confirm those views are unchanged.
