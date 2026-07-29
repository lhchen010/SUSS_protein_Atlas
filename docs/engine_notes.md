# SUSS Atlas Engine — Single-Strain Snakemake pipeline

Converge C. orbiculare single-strain SUSS analysis into a config-driven, one-command-run engine.
Given a secreted protein AF2 structure for a strain (+ optional sequences / RNAseq), produce an interactive atlas HTML end-to-end.

## One-command run
```bash
conda activate <env with snakemake+igraph+leidenalg>
# 1. Prepare inputs (see below)  2. Edit strain/input block in config/config.yaml  3.
snakemake --configfile config/config.yaml --cores 16
# Output: results/<atlas_name>.html + cluster_composition.xlsx + all_families_master.csv + used_config.yaml
```
To switch strains, only change contents of `input/pdb` + `config.strain`, then re-run the same command.

## Inputs
| Path | Required | Description |
|---|---|---|
| `input/pdb/<strain>_<accession>.pdb` | Required | Mature AF2 structure (signal peptide removed), pLDDT>50 |
| `input/seqs.fasta` | Recommended | Mature sequences; if absent, extracted from structure |
| `input/rnaseq.xlsx` | Optional | Two-sheet standard format (id_mapping + expression), see suss_pipeline_design.md §2.1 |
| `config.strain` | Required | Strain code / species / host / phylogeny |

## Parameters (config/config.yaml — lab defaults, normally unchanged)
foldseek_tm=0.5 · tm_symmetric=min · leiden_resolution=1.0 · min_family_size=2 ·
blast_evalue=1e-3 (core_SUSS criterion, adjustable) · qc: pLDDT≥50, length 50–1000

## DAG (see engine_rulegraph.svg)
Full-length branch: qc → foldseek → **cluster (checkpoint)** → classify → F-family workbenches.

Domain branch: qc → local Foldseek segments → D-family graph → domain workbench. Each D-family
workbench stores cropped structures, FoldMason AA/3Di MSA, a FoldMason guide tree, reciprocal-
coverage-controlled sequence subgroups with MAFFT/FastTree, all-pairs US-align, complete-parent
superposition transforms, structural conservation, and optional FoldTree output.

Shared per-protein evidence: FreeSASA, P2Rank, fpocket, ESM-Scan, annotation, DeepTMHMM, and
optional RNA-seq. Pocket and ESM records are calculated per protein and reused by F, D, and
singleton views.

Assembly: master → cards → assemble (server-backed atlas HTML) + composition/download workbooks.

Family count is determined after cluster → checkpoint dynamically expands per-family rules.

## Tool allocation (measured, see pipeline_io_contract.md)
- **4070 (heavy work)**: Foldseek, FoldMason, US-align, MAFFT/FastTree, BLASTp, Rate4Site,
  P2Rank, fpocket, ESM-Scan, FoldTree, InterProScan, EffectorP, and DeepTMHMM.
- **Assembly**: Python plotting, workbooks, ZIP packages, and the server-backed HTML data bridge.

## Validation status
Release-specific live-run metrics and acceptance checks are recorded in the corresponding
`CLAUDE_FOR_SCIENCE_V*.md` handoff. The portal must not be promoted from staging until Python
tests, JavaScript syntax, Snakemake target execution, artifact download checks, and browser
acceptance all pass.

## Directory structure
```
config/config.yaml          lab defaults + strain/input/tools
workflow/Snakefile          18 rules + checkpoint
workflow/scripts/*.py       12 rule scripts (all real)
workflow/builders/          card_layout.md + html_builder.py (complete renderer porting point)
pipeline_io_contract.md     each rule I/O + recovered observations R1–R6
resources/recovered_code.tar.gz  original recovered code
```
