# Configuration

Copy the template and edit it for your machine:

```bash
cp config/config.yaml.template config/config.yaml
```

Then edit `config/config.yaml`:

1. **`strain:`** — your organism code (used as the `<code>_<accession>.pdb` filename prefix),
   species name, host, clade.
2. **`input:`** — `pdb_dir` (folder of AF2 `.pdb` files), optional `seqs_fasta`, optional `rnaseq_xlsx`.
3. **`tools:`** — point each tool at where you installed it. Conda tools (see `../environment.yml`)
   work by bare name once `conda activate suss-atlas` is done; others need an absolute path.
4. **`steps:`** — turn analyses on/off. `qc` and `cluster` are always on; the rest are optional
   (no RNAseq? set `rnaseq: false`. macOS? leave `tools.interproscan` blank — see the main README).
   A blank optional tool means `not_run`; a configured tool that is missing or fails stops the rule.
   `cards` is required whenever `atlas` is enabled.

`output.analysis_scope` accepts `full`, `domain`, or `both`. The portal exposes the same choices
and defaults to `both`. A domain-only run retains the lightweight full-family mapping needed to
explain each segment's parent, but hides the F-family overview from the delivered atlas.

`signals.foldtree_rooting` controls structural-tree rooting. Families with at most
`small_family_max` members use midpoint rooting directly. Larger families use MAD first and
fall back to midpoint rooting only when MAD output is missing or invalid. Every family records
the actual method and recovery status in `<family>_foldtree_status.json`.
`steps.domain_foldtree` separately controls FoldTree for every D family. It is disabled by
default because large domain atlases can contain hundreds of families. FoldMason structural
guide trees remain available when this option is off.

`signals.esm_scope` defaults to `family_representatives`: F-family reference proteins,
protein singletons, and each D family's independently selected structural hub. Thus D-family
ESM never borrows an unrelated F representative. `representatives` retains the smaller
Full-only target set, `domain_members` adds every parent represented in a D family, and
`all_proteins` scans the complete input. The exhaustive modes can require orders of
magnitude more GPU time with `masked-marginals`.
`signals.esm_workers` controls concurrent ESM-Scan subprocesses. Keep it at `1` unless
GPU memory has been validated; the 12 GB RTX 4070 reference configuration uses `2`.

`pocket.p2rank_profile` selects the P2Rank feature profile. Use `alphafold` for the
predicted structures expected by this workflow. AlphaFold PDB files store pLDDT in the
B-factor column, so the generic `default` profile is not appropriate unless the input
contains experimental B-factors. Use `default` only for experimental PDB structures.
The effective profile is saved in `results/used_config.yaml`.

`config.yaml.4070.example` is the exact configuration used to produce the published
*C. orbiculare* atlas, for reference (paths are specific to that machine).

Version 5 has three separate relationship layers:

- `clustering:` defines full-length `F` families using global Foldseek TM score and reciprocal
  structural coverage.
- `domain_clustering:` defines local structural-domain `D` families using Foldseek 3Di+AA
  segment hits. The default retains links with local lDDT at least `0.5`; coverage `0.0` means
  no minimum coverage filter, allowing a shared domain inside different full-length proteins.
- `structural_conservation.pair_threshold:` controls the minimum fraction of FoldMason
  structure-pair subalignments that must support a column before it is colored (default `0.5`).
- `classification:` defines sequence-homologous `S` subgroups and the SUSS divergence labels
  using BLAST e-value plus reciprocal sequence coverage.

See [Foldseek parameters](../docs/FOLDSEEK_PARAMETERS.md) before changing these thresholds.
The defaults are a starting point, not a replacement for reporting the effective parameters in a
paper.

`output.html_mode` controls data packaging:

- `single` embeds structures and downloads for a portable offline HTML file.
- `backend` keeps the interactive atlas in HTML and stores large PDB, ZIP, and Excel artifacts
  under `results/downloads/` for the portal to stream on demand. This is recommended for large
  server runs.

In backend mode every D family receives an Excel workbook, cropped-domain structure ZIP,
complete-parent structure ZIP, and complete package ZIP. P2Rank/fpocket and ESM run once per
protein; domain panels map those parent-protein results to segment coordinates rather than
rerunning surface analysis on an artificial cropped boundary.
