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
First half (all): qc → foldseek → **cluster (checkpoint)** → classify
Second half (all including singletons): sasa_pocket · esm · annotate
Second half (per-family, n≥2): msa → conservation → signature; msa → foldtree (×3 metric)
Assembly: master → cards → assemble (atlas HTML) + composition (Excel)

Family count is determined after cluster → checkpoint dynamically expands per-family rules.

## Tool allocation (measured, see pipeline_io_contract.md)
- **4070 (heavy work)**: foldseek, foldmason, rate4site (must include `-a <ref>`), blastp, P2Rank (java17 env),
  ESM-1b (boltzgen GPU env), FoldTree (base env), InterProScan, EffectorP
- **Local (light work)**: freesasa, fpocket (not installed on 4070), signature, plotting, HTML assembly

## Validation status
- First half: live re-execution matches baseline 71 families / 425 proteins / 92 singletons / 1610 edges / core_SUSS 1109 (front_half_validation.md)
- Second half: family F5 full toolchain end-to-end working, signature reproduces aa_match=1.0 cons_sasa_r=−0.443 (back_half_validation.md)
- DAG: 291 jobs end-to-end parse pass (min_size=2 → 71 families expanded)
- **Pending**: complete interactive builder for assemble_html (html_builder.py) — currently falls back to table index;
  full v19 renderer (artifact c1b281c8) awaiting porting into builders/.

## Directory structure
```
config/config.yaml          lab defaults + strain/input/tools
workflow/Snakefile          18 rules + checkpoint
workflow/scripts/*.py       12 rule scripts (all real)
workflow/builders/          card_layout.md + html_builder.py (complete renderer porting point)
pipeline_io_contract.md     each rule I/O + recovered observations R1–R6
resources/recovered_code.tar.gz  original recovered code
```
