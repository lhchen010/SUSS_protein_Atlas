# SUSS Atlas

**A Snakemake pipeline + web portal for building structural-similarity atlases of secreted
protein (effector) repertoires.**

Given a set of AlphaFold structures for one strain's secreted proteins (plus optional sequences
and RNA-seq), the pipeline clusters them into **SUSS families** — *Sequence-Unrelated,
Structurally Similar* groups — and produces a self-contained interactive HTML atlas. Click any
family to see its structures, structural-similarity matrices, sequence-identity spectrum,
per-residue conservation, binding pockets, mutational tolerance, structural phylogeny, and
RNA-seq expression.

> **SUSS** = proteins that share a fold (Foldseek TM ≥ 0.5) even when their sequences have
> diverged past the point BLAST can detect. Sequence-similar members are kept in the same family
> (important residues may have mutated while the fold is conserved); BLAST is used only to *label*
> the divergence spectrum, never to split families.

---

## What you get

- **Interactive atlas HTML** (self-contained, opens offline in any browser) — network of families;
  click a node for an integrated card.
- **Dual structural-similarity** per family: Foldseek TM **and** an independent US-align TM
  cross-check, with the Foldseek↔US-align Pearson r reported (robustness evidence).
- **`family_summary.xlsx`** — one row per cluster: members, annotation, pocket residues, mean
  TM / sequence identity / SUSS %, RNA-seq per condition; split into clustered vs singleton lists.
- **Web portal** — a browser upload page (internal-network): drop in structures + optional
  RNA-seq, pick parameters, get the atlas. Serves history, logs, and downloads.

![pipeline DAG](docs/engine_rulegraph.svg)

---

## Runs on Linux or macOS

Everything runs on **Linux or macOS** (Intel or Apple Silicon). The core tools install from
conda in one command. The **one exception is InterProScan** (Pfam domains), which is 64-bit
Linux only — on macOS the annotation step still gives protein names and effector calls, just not
Pfam-domain labels (or use Docker / the EBI API — see [INSTALL.md](INSTALL.md)).

---

## Quickstart

```bash
# 1. environment (one command; installs Foldseek, FoldMason, BLAST, Rate4Site, fpocket, ...)
conda env create -f environment.yml
conda activate suss-atlas

# 2. tools installed separately (US-align, FoldTree, EffectorP, ESM-Scan, P2Rank) — see INSTALL.md
#    then point config at them.

# 3. external databases for annotation (Foldseek pdb100 + afdb_swissprot)
bash scripts/setup_databases.sh /path/to/foldseek_db

# 4. try the bundled 100-protein example (needs only the conda env + US-align + FoldTree)
snakemake --configfile examples/config.example.yaml --cores 8 results/example_suss_atlas.html
#    compare against examples/EXPECTED.md

# 5. your own data
cp config/config.yaml.template config/config.yaml   # then edit strain / input / tools
#    put AF2 PDBs in input/pdb as <code>_<accession>.pdb
snakemake --configfile config/config.yaml --cores 16
#    -> results/<atlas_name>.html
```

See [INSTALL.md](INSTALL.md) for per-tool install, [config/README.md](config/README.md) for
configuration, and [examples/EXPECTED.md](examples/EXPECTED.md) for the self-check.

---

## Web portal (optional)

```bash
cd portal && python3 suss_portal.py        # serves on http://0.0.0.0:8600
```

An internal-network upload page: users submit structures (+ optional sequences / RNA-seq),
choose the four exposed parameters, and get the atlas — the pipeline stays on your server, only
the generated HTML is shared. Deployment notes in [portal/DEPLOY.md](portal/DEPLOY.md). The
generated atlas HTML is self-contained, so you can also drop it onto any static web host.

---

## Inputs

| Path | Required | Notes |
|---|---|---|
| `input/pdb/<code>_<accession>.pdb` | yes | mature AF2 structures (signal peptide removed, pLDDT > 50) |
| `input/seqs.fasta` | recommended | mature sequences; if absent, extracted from structures |
| `input/rnaseq.xlsx` | optional | two sheets (`id_mapping` + `expression`); template in `templates/` |
| `config.yaml` | yes | strain metadata, tool paths, step toggles |

---

## Method & tool citations

If you use this pipeline, please cite the tools it builds on:

- **Foldseek** — van Kempen et al., *Nat Biotechnol* 2024 (structural search & clustering)
- **FoldMason** — Gilchrist et al., 2024 (multiple structure alignment)
- **US-align** — Zhang et al., *Nat Methods* 2022 (independent TM cross-check)
- **TM-align** — Zhang & Skolnick, *Nucleic Acids Res* 2005 (TM-score)
- **FoldTree** — Moi et al., 2023 (structural phylogeny)
- **Rate4Site** — Pupko et al., *Bioinformatics* 2002 (per-site conservation)
- **fpocket** — Le Guilloux et al., *BMC Bioinformatics* 2009 (pocket detection)
- **P2Rank** — Krivák & Hoksza, *J Cheminform* 2018 (ligand-binding-site prediction)
- **ESM-1b / ESM-Scan** — Rives et al., *PNAS* 2021 (mutational tolerance)
- **InterProScan** — Jones et al., *Bioinformatics* 2014 (Pfam/CDD/Gene3D domains)
- **EffectorP** — Sperschneider & Dodds, 2022 (effector prediction)
- **DeepTMHMM** — Hallgren et al., 2022 (transmembrane regions)
- **Leiden** — Traag et al., *Sci Rep* 2019 (community detection)
- **Snakemake** — Mölder et al., 2021 (workflow engine)

(Each tool carries its own license; obtain the license-gated ones — InterProScan, EffectorP,
DeepTMHMM — from their official sources. See INSTALL.md.)

---

## License

The pipeline code in this repository is released under the **MIT License** (see [LICENSE](LICENSE)).
Third-party tools invoked by the pipeline are **not** covered by this license and retain their
own terms — install them from their official distributions.
