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

`config.yaml.4070.example` is the exact configuration used to produce the published
*C. orbiculare* atlas, for reference (paths are specific to that machine).

Lab-standard analysis defaults (Foldseek TM ≥ 0.5, Leiden resolution 1.0, min family size 2,
BLAST e-value) live in the `clustering:` / `classification:` blocks and normally don't change —
they define what counts as "structurally similar" and how the SUSS divergence spectrum is labelled.
