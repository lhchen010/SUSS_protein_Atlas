# SUSS Protein Atlas

**A reproducible Snakemake workflow and intranet portal for discovering sequence-unrelated,
structurally similar protein families in secreted effector repertoires.**

[![CI](https://github.com/lhchen010/SUSS_protein_Atlas/actions/workflows/ci.yml/badge.svg)](https://github.com/lhchen010/SUSS_protein_Atlas/actions/workflows/ci.yml)
[![Release](https://img.shields.io/github/v/release/lhchen010/SUSS_protein_Atlas)](https://github.com/lhchen010/SUSS_protein_Atlas/releases/latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-2f855a.svg)](LICENSE)

SUSS Protein Atlas starts from predicted protein structures for one strain, builds
structure-defined families, labels their sequence-divergence spectrum, and integrates independent
structural validation, conservation, pockets, mutational tolerance, phylogeny, annotation, and
optional expression into a self-contained interactive atlas. Version 2 separates true structural
families from unclustered singleton proteins: families remain in the network, while singletons
receive a searchable evidence table and dedicated single-protein viewer.

> **SUSS = Sequence-Unrelated, Structurally Similar.** A SUSS label identifies a structural edge
> that passes the configured Foldseek TM threshold but is not detected by BLAST at the configured
> e-value. BLAST labels divergence; it never splits a structure-defined family.

## At a glance

| Question | Implementation |
|---|---|
| What defines a family? | Foldseek all-vs-all structural edges followed by Leiden community detection |
| What defines a SUSS relationship? | Structural similarity with no BLAST-detected relationship at the configured threshold |
| How is structure independently checked? | Within-family US-align TM matrices, complete-pair validation, and Foldseek/US-align agreement |
| What biological evidence is integrated? | Rate4Site, fpocket, P2Rank, ESM-Scan, FoldTree, InterProScan, EffectorP, DeepTMHMM, and optional RNA-seq |
| How are singletons handled? | As independent proteins with sequence viewer/download, structure, pockets, ESM, annotation, Foldseek database hits, and optional RNA-seq; no artificial singleton cluster or pairwise family analyses |
| What is delivered? | Offline interactive HTML with family and singleton workspaces, Excel summaries, machine-readable tables, per-protein assets, and run provenance |
| How is it orchestrated? | Snakemake checkpoint expansion after the family count becomes known |

## Workflow

```mermaid
flowchart TB
    subgraph INPUTS["1. Inputs and configuration"]
        direction LR
        PDB["AlphaFold PDB structures"]
        FASTA["Mature protein FASTA"]
        RNA["Optional RNA-seq workbook"]
        CFG["Thresholds, tools, and step toggles"]
    end

    subgraph CORE["2. Preflight and structural family discovery"]
        direction LR
        VALIDATE["Format and integrity preflight"]
        QC["Structure QC<br/>length, pLDDT, confident fraction"]
        FOLDSEEK["Foldseek all-vs-all<br/>directional TM scores"]
        GRAPH["Symmetric TM graph<br/>configured min, mean, or max"]
        LEIDEN["Leiden communities<br/>families and singletons"]
    end

    subgraph SPLIT["3. Analysis scope"]
        direction LR
        FAMILIES["Structural families<br/>2 or more members"]
        SINGLETONS["Singleton proteins<br/>independent records"]
    end

    subgraph FAMILY_EVIDENCE["4A. Family evidence"]
        direction LR
        SEQUENCE["BLAST divergence labels"]
        STRUCTURE["Foldseek + US-align matrices"]
        CONSERVATION["FoldMason + Rate4Site"]
        PHYLOGENY["FoldTree phylogeny"]
        FAMILY_SHARED["Pockets, ESM, annotation,<br/>EffectorP, DeepTMHMM, RNA-seq"]
    end

    subgraph SINGLETON_EVIDENCE["4B. Singleton evidence"]
        direction LR
        SINGLE_STRUCT["AlphaFold structure + pLDDT"]
        SINGLE_POCKET["fpocket + P2Rank"]
        SINGLE_ANNO["InterPro + EffectorP + DeepTMHMM"]
        SINGLE_FOLDSEEK["Foldseek PDB100 + AFDB/Swiss-Prot"]
        SINGLE_EXPR["ESM + optional RNA-seq"]
    end

    subgraph INTEGRATION["5. Typed integration"]
        direction LR
        FAMILY_RECORD["Family records<br/>comparative evidence"]
        SINGLE_RECORD["Singleton records<br/>direct evidence only"]
        STATUS["complete / partial / not_run / failed"]
        FAMILY_RECORD --> STATUS
        SINGLE_RECORD --> STATUS
    end

    subgraph OUTPUTS["6. Research outputs"]
        direction LR
        NETWORK["Family network"]
        WORKBENCH["Searchable singleton workbench"]
        EXCEL["Family summary and composition workbooks"]
        TABLES["CSV, JSON, matrices, trees, PDB, and FASTA"]
        PROV["Effective config, hashes, tools,<br/>engine version, and Git commit"]
    end

    PDB --> VALIDATE
    FASTA --> VALIDATE
    RNA --> VALIDATE
    CFG --> VALIDATE
    VALIDATE --> QC --> FOLDSEEK --> GRAPH --> LEIDEN
    LEIDEN --> FAMILIES
    LEIDEN --> SINGLETONS
    FAMILIES --> SEQUENCE
    FAMILIES --> STRUCTURE
    FAMILIES --> CONSERVATION
    FAMILIES --> PHYLOGENY
    FAMILIES --> FAMILY_SHARED
    SINGLETONS --> SINGLE_STRUCT
    SINGLETONS --> SINGLE_POCKET
    SINGLETONS --> SINGLE_ANNO
    SINGLETONS --> SINGLE_FOLDSEEK
    SINGLETONS --> SINGLE_EXPR
    SEQUENCE --> FAMILY_RECORD
    STRUCTURE --> FAMILY_RECORD
    CONSERVATION --> FAMILY_RECORD
    PHYLOGENY --> FAMILY_RECORD
    FAMILY_SHARED --> FAMILY_RECORD
    SINGLE_STRUCT --> SINGLE_RECORD
    SINGLE_POCKET --> SINGLE_RECORD
    SINGLE_ANNO --> SINGLE_RECORD
    SINGLE_FOLDSEEK --> SINGLE_RECORD
    SINGLE_EXPR --> SINGLE_RECORD
    STATUS --> NETWORK
    STATUS --> WORKBENCH
    STATUS --> EXCEL
    STATUS --> TABLES
    STATUS --> PROV

    classDef input fill:#e8f1fb,stroke:#3974a8,color:#17212b,stroke-width:1.5px;
    classDef core fill:#e8f5ec,stroke:#348357,color:#17212b,stroke-width:1.5px;
    classDef split fill:#f0f2f4,stroke:#667985,color:#17212b,stroke-width:1.5px;
    classDef family fill:#fff4d8,stroke:#a97816,color:#17212b,stroke-width:1.5px;
    classDef singleton fill:#e8f5f2,stroke:#267a6a,color:#17212b,stroke-width:1.5px;
    classDef integration fill:#f7e9ee,stroke:#a34d68,color:#17212b,stroke-width:1.5px;
    classDef output fill:#eeeafa,stroke:#6a57a5,color:#17212b,stroke-width:1.5px;
    style INPUTS fill:#f4f8fc,stroke:#8aabc8,stroke-width:1px;
    style CORE fill:#f2f9f4,stroke:#86b398,stroke-width:1px;
    style SPLIT fill:#f7f8f9,stroke:#a7b1b7,stroke-width:1px;
    style FAMILY_EVIDENCE fill:#fffbef,stroke:#c7aa63,stroke-width:1px;
    style SINGLETON_EVIDENCE fill:#f1faf8,stroke:#7ab1a6,stroke-width:1px;
    style INTEGRATION fill:#fcf4f7,stroke:#bd8295,stroke-width:1px;
    style OUTPUTS fill:#f7f5fc,stroke:#9b90c3,stroke-width:1px;
    class PDB,FASTA,RNA,CFG input;
    class VALIDATE,QC,FOLDSEEK,GRAPH,LEIDEN core;
    class FAMILIES,SINGLETONS split;
    class SEQUENCE,STRUCTURE,CONSERVATION,PHYLOGENY,FAMILY_SHARED family;
    class SINGLE_STRUCT,SINGLE_POCKET,SINGLE_ANNO,SINGLE_FOLDSEEK,SINGLE_EXPR singleton;
    class STATUS,FAMILY_RECORD,SINGLE_RECORD integration;
    class NETWORK,WORKBENCH,EXCEL,TABLES,PROV output;
```

The diagram above shows the scientific data flow. The exact rule-level graph remains available
for workflow development:

<details>
<summary>Show the Snakemake rule DAG</summary>

![Snakemake rule DAG](docs/engine_rulegraph.svg)

</details>

## Outputs

| Output | Purpose |
|---|---|
| `results/<atlas_name>.html` | Self-contained family network plus searchable singleton workbench and integrated evidence panels |
| `results/family_summary.xlsx` | Clustered families and singletons with members, evidence, TM statistics, SUSS labels, pockets, and expression |
| `results/cluster_composition.xlsx` | Family membership and annotation composition |
| `results/all_families_master.csv` | Machine-readable integrated family table |
| `results/member_annotation.csv` | Per-protein annotation values and component execution states |
| `results/families/<family>/` | Downloadable family workbook with per-member annotation, Foldseek and US-align matrices, BLAST similarity, FoldMason MSA, complete pocket outputs, FoldTree trees and rooting status, conservation, structures, and RNA-seq |
| Singleton downloads | Mature-sequence FASTA plus an evidence workbook with annotation, Foldseek PDB100/AFDB hits and TM scores, pockets, and RNA-seq; family-only matrices, MSA, FoldTree, conservation, and superposition are intentionally absent |
| `results/used_config.yaml` | Effective configuration plus input hashes, resolved tools, engine version, and Git commit |

## Interactive atlas search

The atlas has two primary views:

- **Cluster network** searches and highlights structure-defined families.
- **Singletons** provides a sortable, paginated table with filters for effector calls, novelty,
  pockets, and transmembrane helices. Selecting a row opens structure, pocket, ESM, RNA-seq, and
  direct annotation evidence for that protein.

Every cluster exposes all mature member sequences as FASTA plus the FoldMason structure-guided
MSA as aligned FASTA. The compact sequence viewer switches between an individual member and the
MSA only when its tab is opened. A singleton exposes its one mature sequence and FASTA download;
MSA is intentionally omitted because it requires at least two proteins.

Both views search locally without a server round trip. Plain text matches accessions, annotations,
InterPro/Pfam terms, Foldseek PDB100 and AFDB/Swiss-Prot hits, EffectorP calls, DeepTMHMM results,
novelty, pockets, and expression. Field prefixes are available for precise queries:

| Prefix | Example |
|---|---|
| `gene:` / `acc:` | `gene:TDZ13877.1` |
| `annotation:` | `annotation:Peroxidase` |
| `effectorp:` | `effectorp:non-effector` |
| `tmr:` / `deeptmhmm:` | `tmr:1` |
| `pdb:` / `afdb:` / `foldseek:` | `afdb:hydrolase` |
| `pocket:` / `rnaseq:` | `pocket:p2rank` |
| `structtm:` | `structtm:0.65` |
| `family:` / `novel:` / `suss:` | `family:F2` |

Matching clusters retain their scientific color and receive an orange outline; non-matches are
visually muted. Press Enter to open a unique result or fit multiple matches, and Escape to reset.

## Quickstart

### 1. Create the environment

```bash
conda env create -f environment.yml
conda activate suss-atlas
```

Install the separately distributed tools needed for the analyses you plan to enable, then set
their paths in the config. See [INSTALL.md](INSTALL.md) for US-align, FoldTree, P2Rank,
ESM-Scan, EffectorP, DeepTMHMM, InterProScan, and database setup.

### 2. Run the validated example

```bash
snakemake \
  --configfile examples/config.example.yaml \
  --cores 8 \
  results/example_suss_atlas.html
```

The bundled dataset contains 100 *Colletotrichum orbiculare* proteins. Compare family counts,
sizes, US-align completeness, and tree outputs with [examples/EXPECTED.md](examples/EXPECTED.md).

### 3. Run a new strain

```bash
cp config/config.yaml.template config/config.yaml
# Edit strain metadata, input paths, enabled steps, and tool paths.

snakemake --configfile config/config.yaml --cores 16
```

## Inputs

| Input | Required | Notes |
|---|---:|---|
| `input/pdb/<code>_<accession>.pdb` | Yes | Mature AlphaFold structures; signal peptide removed |
| `input/seqs.fasta` | Recommended | Mature sequences; missing entries are derived from structure residues |
| `input/rnaseq.xlsx` | No | Sheets `id_mapping` and `expression`; templates are in `templates/` |
| `config/config.yaml` | Yes | Strain metadata, thresholds, step toggles, tool paths, and output settings |

The preflight blocks empty or malformed structure sets before expensive analyses begin. The
current accession parser is designed for versioned GenBank-style IDs; see the handoff report for
known extension points if you need UniProt, JGI, Ensembl, or arbitrary identifiers.

## Failure and optional-step semantics

SUSS Atlas does not use zero matrices, empty trees, or blank files to represent failed tools.

| Situation | Recorded behavior |
|---|---|
| Step disabled | Typed output with `not_run`; unavailable scientific calls remain null |
| Optional component omitted | Parent analysis may be `partial`; completed components remain usable |
| Configured tool missing or unsuccessful | Rule fails with a diagnostic log |
| Tool reports success but expected output is incomplete | Rule fails validation |
| Annotation evidence incomplete | `novel` remains null rather than becoming a false positive |

This distinction is carried into member tables, family summaries, the atlas, and provenance.

## Web portal

The included portal provides an internal-network upload workflow with persistent history, live
logs, original-input downloads, run parameters, and generated atlas/Excel downloads.

```bash
cd portal
python3 suss_portal.py
```

The portal enforces upload and expanded-archive limits, safe archive extraction, a configurable
active-job limit, CSRF-protected deletion, and UUID job IDs. It is intended for a trusted intranet
or Tailscale network; it does not provide public-service authentication or user isolation. See
[portal/DEPLOY.md](portal/DEPLOY.md).

## Platform support

The core workflow supports Linux and macOS. InterProScan is the main native-platform exception:
its official local distribution is 64-bit Linux. On macOS, leave `tools.interproscan` blank or use
Docker/the EBI service; fold and effector evidence can still run, while novelty remains uncalled
unless all required evidence is complete.

## Documentation

| Document | Contents |
|---|---|
| [INSTALL.md](INSTALL.md) | Environment, external tools, databases, and DeepTMHMM compatibility |
| [config/README.md](config/README.md) | Configuration fields and step behavior |
| [examples/EXPECTED.md](examples/EXPECTED.md) | Reproducible 100-protein acceptance baseline |
| [docs/pipeline_io_contract.md](docs/pipeline_io_contract.md) | Rule inputs, outputs, parameters, and contracts |
| [docs/CLAUDE_FOR_SCIENCE_V2.0.0_HANDOFF.md](docs/CLAUDE_FOR_SCIENCE_V2.0.0_HANDOFF.md) | v2 singleton workbench design, validation, deployment evidence, and Claude acceptance checklist |
| [portal/DEPLOY.md](portal/DEPLOY.md) | Intranet portal deployment and operational scope |

## Citation and licenses

SUSS Atlas integrates Foldseek, FoldMason, US-align, TM-align, FoldTree, Rate4Site, fpocket,
P2Rank, ESM-1b/ESM-Scan, InterProScan, EffectorP, DeepTMHMM, Leiden, and Snakemake. Cite the
underlying methods used in your enabled analysis; installation notes are collected in
[INSTALL.md](INSTALL.md).

<details>
<summary>Methods and tools to cite</summary>

- **Foldseek** — van Kempen et al., *Nature Biotechnology* 2024
- **FoldMason** — Gilchrist et al., 2024
- **US-align** — Zhang et al., *Nature Methods* 2022
- **TM-align** — Zhang and Skolnick, *Nucleic Acids Research* 2005
- **FoldTree** — Moi et al., 2023
- **Rate4Site** — Pupko et al., *Bioinformatics* 2002
- **fpocket** — Le Guilloux et al., *BMC Bioinformatics* 2009
- **P2Rank** — Krivák and Hoksza, *Journal of Cheminformatics* 2018
- **ESM-1b / ESM-Scan** — Rives et al., *PNAS* 2021
- **InterProScan** — Jones et al., *Bioinformatics* 2014
- **EffectorP** — Sperschneider and Dodds, 2022
- **DeepTMHMM** — Hallgren et al., 2022
- **Leiden** — Traag et al., *Scientific Reports* 2019
- **Snakemake** — Mölder et al., 2021

</details>

The pipeline code is released under the [MIT License](LICENSE). Third-party tools, models, and
databases retain their own licenses; obtain license-gated components from their official sources.
