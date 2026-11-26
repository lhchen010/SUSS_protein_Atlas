# SUSS Protein Atlas — Methodology Overview

## 1. Conceptual Framework: Separating F, D, and S

Secreted pathogen effectors frequently retain conserved tertiary folds or reusable structural
domains despite extreme primary-sequence divergence (Sequence-Unrelated, Structurally Similar;
**SUSS**). SUSS Protein Atlas explicitly separates three orthogonal evolutionary groupings:

| Layer | Label | Biological Question | Method |
|---|---|---|---|
| Full-length fold family | `F` | Which proteins share an overall global fold? | All-vs-all Foldseek TM-align (`--alignment-type 1`) + reciprocal coverage + Leiden clustering |
| Local structural domain family | `D` | Which proteins share a reusable local structural module? | All-vs-all local Foldseek 3Di+AA (`--alignment-type 2`) + segment-graph clustering |
| Sequence-homologous subgroup | `S` | Which structural relatives retain detectable sequence homology? | Reciprocal-coverage-controlled BLASTp within each `F` (full-length) or `D` (cropped segment) family |

A structural edge is labeled **SUSS** when it passes the configured Foldseek TM and coverage
thresholds without a reciprocal-coverage BLASTp link at the configured e-value. BLASTp
annotates the sequence-divergence spectrum; it never splits a structure-defined family.

## 2. Structure QC and Format Normalization

Input structures (`.pdb`, `.cif`, `.mmcif`, plain or `.gz`) are parsed at preflight and QC
(`workflow/scripts/qc.py`) to compute CA residue count, mean pLDDT, and fraction of confident
residues (pLDDT ≥ 70). Structures passing configured length and pLDDT bounds are staged into
`results/qc_pdb/` (`workflow/scripts/stage_qc_structures.py`), automatically converting mmCIF
inputs into standardized PDB coordinates so every downstream tool operates on an identical,
QC-verified structure set.

## 3. Independent Structural Validation & Alignment

- **Global fold validation (`F`)**: Within every `F` family, pairwise TM-scores are independently
  recomputed with **US-align** (`O(m^2)` per family) and cross-checked against Foldseek TM-scores
  to quantify agreement and flag missing or discordant pairs. Missing pairs remain `NaN` and
  never contribute zeros to family means or hub selection.
- **Hub selection**: Each family hub is chosen as the member maximizing mean measured within-family
  TM-score, and placed first so FoldMason alignments, conservation coordinates, pocket references,
  and superposition frames share one consistent reference protein.
- **Domain workbench (`D`)**: For each `D` family, matched residue intervals are cropped into
  domain structures and sequences for independent **FoldMason** AA/3Di MSAs, **US-align**
  matrices, **BLASTp** subgroups, **MAFFT** MSAs, **FastTree** phylogenies, and optional
  **FoldTree** structural guide trees, without borrowing `F`-family evolutionary results.

## 4. Conservation, Pockets, Mutational Tolerance & Annotation

- **Structural conservation**: Official **FoldMason** per-column LDDT scores; alignment columns
  lacking the configured fraction of supporting structure pairs (`pair_threshold`, default `0.5`)
  remain explicitly unscored (grey).
- **Sequence evolutionary conservation**: **Rate4Site** is executed strictly on eligible
  reciprocal-coverage `S` subgroups (`min_rate4site_sequences`, default `4`). Sequence-unrelated
  members are never forced into a single sequence alignment.
- **Surface accessibility & pockets**: **FreeSASA**, **fpocket**, and **P2Rank** (defaulting to
  its AlphaFold profile so pLDDT B-factors are not treated as crystallographic B-factors) run on
  complete parent structures; `D` workbenches map parent pocket residues onto domain intervals.
- **Mutational tolerance**: Optional **ESM-1b** (`ESM-Scan`) scores run in full parent-protein
  context.
- **Tri-state annotation & novelty**: Functional annotation integrates **InterProScan**,
  **Foldseek PDB100**, **Foldseek AFDB/Swiss-Prot**, **EffectorP**, and optional **DeepTMHMM**.
  Novelty (`novel=True`) is called only when all required domain and fold searches complete
  (`annotation_status == "complete"`); incomplete evidence yields `novel = null`.
