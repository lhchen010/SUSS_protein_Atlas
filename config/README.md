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

`signals.foldtree_rooting` controls structural-tree rooting. Families with at most
`small_family_max` members use midpoint rooting directly. Larger families use MAD first and
fall back to midpoint rooting only when MAD output is missing or invalid. Every family records
the actual method and recovery status in `<family>_foldtree_status.json`.

`config.yaml.4070.example` is the exact configuration used to produce the published
*C. orbiculare* atlas, for reference (paths are specific to that machine).

Version 3 has three separate relationship layers:

- `clustering:` defines full-length `F` families using global Foldseek TM score and reciprocal
  structural coverage.
- `domain_clustering:` defines local structural-domain `D` families using Foldseek 3Di+AA
  segment hits.
- `classification:` defines sequence-homologous `S` subgroups and the SUSS divergence labels
  using BLAST e-value plus reciprocal sequence coverage.

See [Foldseek parameters](../docs/FOLDSEEK_PARAMETERS.md) before changing these thresholds.
The defaults are a starting point, not a replacement for reporting the effective parameters in a
paper.
