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
classification.blast_evalue=1.0    # ⚠ see R1 below
qc.min_plddt=50, max_length=1000, min_length=50
signals.foldtree_metrics=[foldtree,alntmscore,lddt], esm_model=esm1b, plm_strategy=masked_marginals
```

## 11 Rule Contract

| # | rule | Input | Output | Tool | config Parameter | Scope |
|---|---|---|---|---|---|---|
| 1 | qc | input_dir/*.pdb | qc.csv (acc,length,mean_plddt,frac_conf,pass) | custom (CA B-factor parsing) | min_plddt,max_length,min_length | all |
| 2 | foldseek | QC-passing PDB | foldseek_allvsall.tsv | foldseek (easy-search) | — | all |
| 3 | cluster | foldseek tsv + qc.csv | families.csv, members.csv, edges.csv | igraph+leidenalg | foldseek_tm,tm_symmetric,leiden_resolution,min_family_size | all |
| 4 | classify | blastp tsv + edges.csv | classification.csv (core_SUSS/…) | blastp + merge | blast_evalue | all |
| 5 | msa | family members PDB per family | {fam}.aln, {fam}.fasta | FoldMason | — | families n≥2 |
| 6 | conservation | {fam}.aln + family tree | {fam}_r4s.res | Rate4Site | — | families n≥2 |
| 7 | sasa_pocket | per-protein PDB (including singletons) | sasa_all.csv, pockets.json | freesasa+fpocket+P2Rank (java17) | — | all |
| 8 | esm | per-sequence (including singletons) | {acc}_esm.csv | ESM-1b (esmscan.py) | esm_model,plm_strategy | all |
| 9 | foldtree | family members PDB per family | {fam}_{metric}.nwk ×3 | FoldTree (snakemake pipeline) | foldtree_metrics | families n≥2 |
| 10 | annotate | all protein PDB+seq (including singletons) | member_annotation.csv, cluster_annotation.csv | Foldseek (pdb100/afdb)+InterProScan+EffectorP+DeepTMHMM | — | all |
| 11 | rnaseq | rnaseq.xlsx + members | {fam}_expression.csv | pandas (log2 CPM) | — | all (if data present) |
| → | signature | r4s + sasa + pockets | {fam}_signature.csv, B-factor PDB | custom | — | families n≥2 |
| → | cards | all upstream | co_card_{fam}.png | matplotlib | — | per-family+singleton |
| → | assemble | cards+structure+matrices+trees | atlas.html, cluster_composition.xlsx, master.csv | custom | html_mode (single/backend) | merge |

## Key Observations from Recycled Code (must handle during refactor)

- **R1 — BLAST threshold inconsistency**: Recycled `classification.py` uses `evalue<=1e-3` for blast_detected judgment, but design lab default is `blast_evalue=1.0`. → Rule script must consume config.classification.blast_evalue, default 1.0 (consistent with design); 1e-3 was the stricter threshold used in original classification. **User confirmation required**: Should core_SUSS determination use 1.0 or 1e-3? (Affects SUSS% — 1.0 more conservatively flags "BLAST-detected" → fewer core_SUSS).
- **R2 — QC coupled to Drive download**: Recycled `qc.py` first half is Google Drive MCP download, second half is PDB parsing. → Engine's qc rule performs **parsing+filtering only**; input_dir is prepared by the upload portal; Drive download is the portal's responsibility, not part of the pipeline.
- **R3 — Hardcoded artifact absolute paths**: Recycled code has many hardcoded `/…/artifacts/…/vXXXX_*.csv` paths. → Replace entirely with rule input/output wildcards; zero absolute paths.
- **R4 — min_family_size in recycled code is 3**: `network_families.py` writes `sizes>=3`; lab default changed to 2. → Consume config.clustering.min_family_size.
- **R5 — Symmetric TM**: Recycled code `df[[qtmscore,ttmscore]].min(axis=1)` is correct (= design's tm_symmetric=min); if config provides max/mean, require branching.
- **R6 — Leiden**: `RBConfigurationVertexPartition, resolution_parameter, seed=42`, weight=tm. Retain; resolution consumed from config.
