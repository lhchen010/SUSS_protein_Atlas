# SUSS Atlas Engine — Pipeline I/O Contract

Single-strain Snakemake engine rule-by-rule input/output/tool/config parameter specification. This is the rule script specification.

## Inputs (user-provided, via upload portal or config)
- `input_dir/` — mature AF2 PDB, filename `<strain>_<accession>.pdb` (or accession-only; strain supplied by config)
- `input_seqs.fasta` — mature sequences (optional; if absent, extracted from structure)
- `rnaseq.xlsx` — two-sheet standard format (optional; see design document §2.1)
- `strain_meta.yaml` — strain codes/species/host/phylogeny (see §2.0)

## Config Parameters (config/config.yaml — lab default)
```
clustering.foldseek_tm=0.5, tm_symmetric=min, leiden_resolution=1.0, min_family_size=2
clustering.whole_fold_min_coverage=0.5
domain_clustering.evalue=1e-3, min_probability=0.5, min_aligned_residues=40, min_lddt=0.5
classification.blast_evalue=1e-3, min_reciprocal_coverage=0.5
qc.min_plddt=50, max_length=1000, min_length=50
signals.foldtree_metrics=[foldtree,alntmscore,lddt], esm_model=esm1b, plm_strategy=masked_marginals
```

## Core Rule Contract

| # | rule | Input | Output | Tool | config Parameter | Scope |
|---|---|---|---|---|---|---|
| 1 | qc | input_dir/*.pdb | qc.csv (acc,length,mean_plddt,frac_conf,pass) | custom (CA B-factor parsing) | min_plddt,max_length,min_length | all |
| 2 | foldseek | QC-passing PDB | foldseek_allvsall.tsv | Foldseek easy-search, global alignment | alignment_type, exhaustive_search | all |
| 3 | cluster | Foldseek TSV + qc.csv | families.csv, members.csv, edges.csv, family_members/ | igraph + Leiden | foldseek_tm, reciprocal coverage, tm_symmetric, resolution, min size | all |
| 4 | foldseek_domains | QC-passing PDB | foldseek_domains.tsv | Foldseek easy-search, local 3Di+AA | e-value, sensitivity | all |
| 5 | domain_cluster | local Foldseek TSV | domain_families/members/edges/cross_edges.csv | interval consolidation + Leiden | probability, aligned length, local lDDT, optional alignment TM/coverage | all |
| 6 | domain_workbench | D members/edges + sequences/PDB | domain_workbench.json, domain identity matrices, structural/sequence conservation | FoldMason + US-align + MAFFT + FastTree + Rate4Site | configured tools and minimum subgroup size | D families |
| 7 | classify | BLAST TSV + retained F edges | classification.csv (core_SUSS/…) | qualifying-hit merge | e-value + reciprocal sequence coverage | F edges |
| 8 | sequence_subgroups | BLAST TSV + F members | sequence_subgroups/edges.csv | connected components | e-value + reciprocal sequence coverage | F families |
| 9 | msa | F-family member PDB | {fam}.aln, {fam}.fasta | FoldMason | — | F families n≥2 |
| 10 | structural_conservation | FoldMason AA/3Di MSA | structural_conservation.csv/PDB | FoldMason msa2lddtjson | pair_threshold | F families n≥2 |
| 11 | sequence_analysis | S subgroup sequences | MAFFT MSA, FastTree, Rate4Site, status | MAFFT + FastTree + Rate4Site | min subgroup size | eligible S subgroups |
| 12 | sasa_pocket | per-protein PDB, members.csv | sasa_all.csv, pockets.json keyed by family or singleton accession | freesasa+fpocket+P2Rank (java17) | — | all proteins; pocket detection on family references and every singleton |
| 13 | esm | canonical sequences, members.csv | esm_all.csv keyed by family or singleton accession | ESM-1b (esmscan.py) | esm_model,plm_strategy | family references and every singleton |
| 14 | foldtree | family members PDB per family | {fam}_{metric}.nwk ×3 | FoldTree (snakemake pipeline) | foldtree_metrics | families n≥2 |
| 15 | annotate | all protein PDB+seq (including singletons) | member_annotation.csv, cluster_annotation.csv | Foldseek (pdb100/afdb)+InterProScan+EffectorP+DeepTMHMM | — | all |
| 16 | rnaseq | rnaseq.xlsx + members | {fam}_expression.csv | pandas (log2 CPM) | — | all (if data present) |
| → | signature | r4s + sasa + pockets | {fam}_signature.csv, B-factor PDB | custom | — | families n≥2 |
| → | cards | all upstream | co_card_{fam}.png | matplotlib | — | per-family |
| → | structure_db | QC-passing PDB + F/D/S assignments | structure_db/atlas*, structure_search_index.csv | Foldseek createdb | — | all |
| → | assemble | master, F/D/S assignments, domain workbench, annotation, pockets, optional ESM/RNA, family assets | atlas.html + optional downloads/ artifacts | custom | html_mode (single/backend) | F, D, and singleton workbenches |

## Key Observations from Recycled Code (must handle during refactor)

- **R1 — BLAST relationship definition (resolved)**: both `classify` and
  `sequence_subgroups` consume the configured e-value and reciprocal-coverage thresholds. They
  filter HSPs first, then choose the best qualifying HSP per undirected pair. The shipped default
  is `1e-3` and `0.5` reciprocal coverage.
- **R2 — QC coupled to Drive download**: Recycled `qc.py` first half is Google Drive MCP download, second half is PDB parsing. → Engine's qc rule performs **parsing+filtering only**; input_dir is prepared by the upload portal; Drive download is the portal's responsibility, not part of the pipeline.
- **R3 — Hardcoded artifact absolute paths**: Recycled code has many hardcoded `/…/artifacts/…/vXXXX_*.csv` paths. → Replace entirely with rule input/output wildcards; zero absolute paths.
- **R4 — min_family_size in recycled code is 3**: `network_families.py` writes `sizes>=3`; lab default changed to 2. → Consume config.clustering.min_family_size.
- **R5 — Symmetric TM**: Recycled code `df[[qtmscore,ttmscore]].min(axis=1)` is correct (= design's tm_symmetric=min); if config provides max/mean, require branching.
- **R6 — Leiden**: `RBConfigurationVertexPartition, resolution_parameter, seed=42`, weight=tm. Retain; resolution consumed from config.
