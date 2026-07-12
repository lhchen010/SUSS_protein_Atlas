# Installation

The pipeline runs on **Linux or macOS** (Intel or Apple Silicon). It needs (1) a conda
environment for the core tools, (2) a few tools installed separately, and (3) two external
structural databases. This file covers all three.

---

## 1. Core environment (conda — one command)

```bash
conda env create -f environment.yml     # mamba env create -f environment.yml is faster
conda activate suss-atlas
```

This installs the workflow engine (Snakemake), clustering (igraph/leidenalg), and every
**conda-installable** tool. All of these have macOS builds (Intel + Apple Silicon) and Linux builds:

| Tool | Role | Channel |
|---|---|---|
| **Foldseek** | all-vs-all structural search + clustering (the core) | bioconda |
| **FoldMason** | multiple-structure alignment (MSA / superpose order) | bioconda |
| **Rate4Site** | per-residue conservation | bioconda |
| **BLAST+** | sequence-identity divergence spectrum | bioconda |
| **fpocket** | binding-pocket detection | bioconda |
| **mTM-align** | optional multiple-alignment alternative | bioconda |

---

## 2. Tools installed separately

These are not on conda (GitHub source or license-gated). Install each, then put its path in
`config/config.yaml` under `tools:`. **All are cross-platform except InterProScan.**

### US-align — independent-algorithm TM cross-check  *(macOS + Linux)*
Single C++ file, compiles in seconds.
```bash
git clone https://github.com/pylelab/USalign && cd USalign && make
# or: g++ -O3 -ffast-math -o USalign USalign.cpp
```
Set `tools.usalign: /path/to/USalign/USalign` (or put it on PATH).

### FoldTree — structural phylogeny  *(macOS + Linux)*
```bash
git clone https://github.com/DessimozLab/fold_tree
```
Set `tools.fold_tree: /path/to/fold_tree`. (It bundles QuickTree/FastME/MAD — no extra installs.)

### EffectorP — effector prediction  *(macOS + Linux, Python)*
```bash
git clone https://github.com/JanaSperschneider/EffectorP-3.0
```
Set `tools.effectorp: /path/to/EffectorP.py`.

### ESM-Scan — mutational tolerance  *(macOS + Linux; PyTorch)*
```bash
git clone https://github.com/xuebingwu/ESM-Scan   # upstream: Xuebing Wu lab, Columbia
pip install torch fair-esm          # CPU works; Apple-Silicon MPS or CUDA if available
```
Set `tools.esm_scan: /path/to/ESM-Scan/esmscan.py` and `tools.esm_python:` to a python with torch+esm.
(Optional step — set `steps.esm: false` to skip.)

### P2Rank — ligand-binding-site prediction  *(macOS + Linux; needs Java 17)*
```bash
conda install -c bioconda p2rank      # or download from https://github.com/rdkit/p2rank/releases
conda install -c conda-forge openjdk=17
```
Set `tools.p2rank: prank`. (Optional — fpocket alone covers pockets if P2Rank is absent.)

### InterProScan — Pfam / CDD / Gene3D domains  *(⚠️ 64-bit LINUX ONLY)*
InterProScan ships precompiled Linux binaries and **does not run natively on macOS**.
Three options:

- **Linux:** download from https://www.ebi.ac.uk/interpro/download/InterProScan/ and set
  `tools.interproscan: /path/to/interproscan.sh`.
- **macOS — skip it:** leave `tools.interproscan: ""`. The `annotate` step still runs and
  populates **protein names** (Foldseek vs AFDB-SwissProt) and **effector** calls (EffectorP);
  only Pfam-domain labels (GH16, CFEM, …) are omitted. The engine handles the blank cleanly.
- **macOS — full domains via Docker:** run the official InterProScan Docker image, or submit
  the FASTA to the EBI InterProScan web API, then drop the resulting TSV into
  `results/anno/interproscan.tsv` before the `annotate` step.

### DeepTMHMM — transmembrane regions  *(optional; only fills `n_TMR`)*
Optional and easy to skip — leave `tools.deeptmhmm: ""`. If you want it, install per its
BioLib/DeepTMHMM instructions and set the path plus `tools.deeptmhmm_python:` (a python with
torch+h5py). Nothing else depends on it.

---

## 3. External databases

Foldseek annotation needs two reference databases (PDB100 + AFDB/Swiss-Prot). Download them with:

```bash
bash scripts/setup_databases.sh /path/to/foldseek_db
```

Then set `tools.foldseek_db_dir: /path/to/foldseek_db` in the config. See that script for size
requirements (tens of GB). InterProScan bundles its own data with its release download.

---

## Quick sanity check

```bash
conda activate suss-atlas
foldseek version && foldmason version && blastp -version && rate4site -h 2>&1 | head -1
snakemake --version
```

Then run the bundled example (see `examples/`) to confirm your environment reproduces the
expected output before running your own data.
