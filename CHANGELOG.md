# Changelog

All notable changes to **SUSS Protein Atlas** are documented here. Detailed historical
development and audit notes are preserved under [`docs/archive/`](docs/archive/).

## [5.3.0] — Architecture Hardening, Flexible IDs & mmCIF Support

### Added
- **Centralized flexible protein ID parser** (`workflow/scripts/id_utils.py`): natively supports
  GenBank (`TDZ13209.1`), RefSeq (`XP_001234.1`), AlphaFold DB (`AF-P12345-F1-model_v4`),
  UniProt (`sp|...|` / bare accessions), custom user regex (`input.accession_pattern`), and
  safe filename-stem fallback, eliminating 11 duplicated regexes across the workflow and portal.
- **Native `.cif` / `.mmcif` (+ `.gz`) input support** (`workflow/scripts/structure_io.py`):
  preflight, QC, and portal archive uploads accept AlphaFold 3, Boltz-2, and AFDB v4 mmCIF
  structures and stage standardized `.pdb` files into `results/qc_pdb/` for uniform downstream
  execution across all structure-consuming Snakemake rules.
- **Reproducible container definitions**: added `docker/Dockerfile`, `docker/Apptainer.def`,
  `config/config.yaml.container.example`, and `scripts/build_release_tarball.sh` (macOS
  AppleDouble-free archive builder).
- **Portal hardening & service templates**: added `portal/multipart_form.py` (stdlib-only
  Python 3.11–3.14+ multipart parser), `portal/suss-portal.service` (`systemd --user` unit),
  and Tailscale-interface default binding in `portal/launch_4070.sh`.
- **Methods documentation**: added `docs/METHODS.md` summarizing the F/D/S clustering model,
  validation statistics, conservation masking, and failure semantics.

### Changed
- **Modularized `workflow/builders/html_builder.py`**: extracted `pdb_io.py`, `svg_charts.py`,
  and `export_bundles.py` while preserving backward-compatible re-exports.
- **Testable `domain_workbench.py`**: wrapped Snakemake top-level execution in `main(snakemake)`.
- **Documentation structure**: archived 22 historical handoff notes to `docs/archive/`.

### Fixed
- **`annotate.py` thread oversubscription**: replaced hardcoded `--threads 24` with Snakemake
  `threads` (configurable via `resources.annotate_threads`, capped by `--cores`).
- **Python 3.13+ compatibility**: removed deprecated `cgi.FieldStorage` from `portal/suss_portal.py`
  (zero deprecation warnings under `pytest -W error::DeprecationWarning`).

## [5.2.5] — Responsive Protein Evidence Tables
- Made singleton and domain-unassigned evidence tables responsive on narrow viewports while
  keeping `EffectorP`, `TMR`, `pLDDT`, and length as independent sortable columns.

## [5.2.4] — Domain-Unassigned Workbench
- Added paginated workbench for proteins without a retained D-family assignment, preserving
  full structure, RNA-seq, annotation, pocket, and local-hit diagnostic evidence.

## [5.2.3] — Configurable Structure QC Length & pLDDT
- Exposed `qc.min_length`, `qc.max_length`, and `qc.min_plddt` in the portal and routed both
  global and local Foldseek searches through staged QC-passing structures.

## [5.2.0] — Scientific-Console UI Refresh
- Unified visual design across the portal and linked F/D/singleton atlas workbenches.

## [5.1.0 – 5.1.3] — Independent Domain Evolution & Validation Audits
- Added independent D-segment BLASTp, MAFFT, FastTree, Rate4Site, and borderline/filtered
  local-hit diagnostics; tightened missing-pair statistics and conservation pair-support masking.

## [5.0.0 – 5.0.2] — Dual-Axis F/D Navigation & FoldMason Conservation
- Separated full-length (`F`) and local structural domain (`D`) axes; mapped FoldMason
  per-column LDDT structural conservation onto all aligned segments.

## [1.0.0 – 4.1.0] — Core Snakemake Engine, Portal, Pockets & Validation
- Established global Foldseek + Leiden clustering, US-align cross-validation, FoldTree,
  fpocket/P2Rank pockets, ESM-1b tolerance, InterProScan/EffectorP/DeepTMHMM annotation,
  RNA-seq integration, and server-backed lazy-loading intranet portal.
