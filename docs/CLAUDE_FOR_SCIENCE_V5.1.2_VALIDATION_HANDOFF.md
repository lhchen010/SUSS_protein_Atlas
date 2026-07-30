# SUSS Protein Atlas v5.1.2 validation handoff

## Decision

Release candidate status: **PASS**.

SUSS Protein Atlas was run on a ten-protein subset of the original
*Colletotrichum orbiculare* F23 family. Foldseek, FoldMason, US-align, FoldTree,
fpocket, and P2Rank were then run independently on the same 4070 host and inputs.
The workflow and direct-tool results passed all 29 programmed comparisons.

The audit found and fixed seven implementation defects. The final rerun,
functional checks, local tests, download inspection, and browser acceptance all
passed, subject to the explicitly listed residual visual check.

## Acceptance atlas

Portal run:

```text
20260730-v5.1.2-F23-audit
```

Portal URL:

```text
http://100.80.77.29:8600/atlas?id=20260730-v5.1.2-F23-audit
```

4070 audit workspace:

```text
/home/claude/suss_audits/20260730-v5.1.1-F23
```

The directory name records the release from which the audit started. The final
page title and validated code are v5.1.2.

## Audit cohort

The cohort is the complete ten-member F23 family from the validated v5.1.1
reference atlas:

```text
TDZ16948.1
TDZ18954.1
TDZ19006.1
TDZ19854.1
TDZ19858.1
TDZ20033.1
TDZ21966.1
TDZ23942.1
TDZ24074.1
TDZ24636.1
```

Inputs:

| Input | Source | SHA-256 |
|---|---|---|
| 10 AlphaFold-style PDB files | v5.1.1 reference run | Individual hashes are recorded in the execution log |
| `seqs.fasta` | F23 mature sequences | `10c41ee7c25cd27c1c773bfde31e12a0cbab990d7cfa5330052d55683afbf376` |
| `rnaseq.xlsx` | Original 13,253-gene workbook, 5 conditions x 3 replicates | `b68c8df585036ce15e1bb2c0c613c2366a85a1812f305220c872d5b54bccdaa5` |

The workflow calculated five condition means for all ten selected proteins. All
50 values matched direct means calculated from the raw replicate columns; the
maximum absolute floating-point difference was
`7.28e-12`.

## 4070 environment

| Component | Version or executable |
|---|---|
| Host | `mpfi-linux1`, Linux 6.8.0-124 |
| Workflow Python | 3.11.15 |
| Workflow Snakemake | 9.23.1 |
| FoldTree Snakemake | 8.25.3 |
| Foldseek | 9.427df8a |
| FoldMason | 1.763a428 |
| US-align | 20241108 |
| fpocket | 4.0 |
| P2Rank | 2.5.1, `alphafold` profile in the configured Java 17 environment |

The audit began from Git commit
`ca2cfbd3d4c5ec3386c0580d61503d1cd59bfb62` and applied the v5.1.2 fixes
described below.

## Final atlas results

| Result | Value |
|---|---:|
| QC-passing proteins | 10/10 |
| Full-length families | 1 (`F0`) |
| F0 members | 10 |
| Retained full-length structural edges | 43 |
| `core_SUSS` edges | 42 |
| `moderate_paralog` edges | 1 |
| Domain families | 2 (`D0`, `D1`) |
| Retained domain segments | 7 |
| Domain-unassigned or filtered proteins | 3 |
| RNA-seq mappings | 10/10 |
| Complete annotation records | 10/10 |
| fpocket completed | 10/10 |
| P2Rank completed | 10/10 |
| ESM representatives scanned | 2/2 |

`D0` contains five segments and `D1` contains two. ESM was run on the F-family
representative and both D-family hubs; one protein serves as both F0 and D0 hub,
so two unique parent proteins were scanned.

The largest reciprocal-coverage sequence subgroup in F0 contains only two
proteins. With `min_rate4site_sequences: 4`, sequence Rate4Site is correctly
reported as `not_applicable`; this is not a missing output.

## Independent comparisons

All checks in `controls/comparison_results.json` passed:

| Comparison | Result |
|---|---|
| Direct full-length Foldseek vs workflow | PASS, 100 directed rows after deterministic normalization |
| Direct local/domain Foldseek vs workflow | PASS, 67 directed rows |
| Direct FoldMason AA alignment | PASS, byte-identical |
| Direct FoldMason 3Di alignment | PASS, byte-identical |
| Direct FoldMason guide tree | PASS, byte-identical |
| Direct US-align vs workflow | PASS, all 45 pairs; maximum pre-display rounding delta `0.0005` |
| Direct FoldTree tree | PASS, all 45 pairwise patristic distances identical |
| Direct FoldTree alignment-TM tree | PASS, all 45 pairwise patristic distances identical |
| Direct FoldTree lDDT tree | PASS, all 45 pairwise patristic distances identical |
| Direct fpocket vs workflow | PASS, 10/10 proteins |
| Direct P2Rank vs workflow | PASS, 10/10 proteins, `alphafold` profile |
| Family pocket alias | PASS, F0 resolves to `TDZ20033.1` |
| Raw RNA replicate means vs workflow | PASS, 50/50 values |
| Annotation completion | PASS, 10/10 |
| F23 reference membership and 43 edges | PASS |
| Reference classification, matrices, conservation, alignments, trees, RNA and annotation | PASS |

Foldseek row order is not treated as scientific data. Direct and workflow TSVs
were compared after deterministic row normalization. FoldTree Newick child/root
serialization may differ while representing the same tree; therefore every
pairwise patristic distance was compared instead of requiring byte-identical
Newick text.

## Defects found and fixed

### 1. Family analysis inputs used an obsolete directory

The checkpoint correctly writes member lists to:

```text
results/family_members/F*.members.txt
```

Structural conservation, signature fallback, pocket family aliases, and ESM
family aliases still looked under the analysis directory in some paths. This
caused a missing-file failure or silently missing aliases.

All four consumers now receive the checkpoint-owned member directory or member
file as an explicit Snakemake input. This also makes dependency tracking correct.

### 2. Successful zero-pocket P2Rank results were incomplete

P2Rank can complete successfully and report zero pockets. The workflow recorded
`p2rank_status: complete` but omitted the typed `p2rank` result object in this
case, making the atlas display it as unavailable.

Successful zero-pocket results now contain:

```json
{
  "top_score": null,
  "top_probability": null,
  "n_pockets": 0,
  "lining_residues": [],
  "pockets": []
}
```

In this cohort, TDZ19854.1, TDZ20033.1, and TDZ21966.1 legitimately have zero
P2Rank pockets. Direct P2Rank confirms those results.

### 3. Undirected pair metadata depended on input row direction

Protein IDs were sorted into canonical pair order after selecting Foldseek or
BLAST rows, but directional fields were not swapped with them. Query/target TM,
coverage, and length metadata could therefore be attached to the wrong protein,
and reciprocal rows with equal primary scores could depend on input order.

Pair orientation now swaps all directional fields together and uses deterministic
tie-breakers. Added tests reverse reciprocal-row order and verify identical
output.

### 4. Small disconnected domain networks overlapped

A two-to-six-node domain overview with no bridge edges could keep running physics
and draw nodes on top of each other. Such networks now use fixed horizontal
positions and disabled physics. Larger or connected networks retain the normal
interactive layout.

### 5. Structure ZIP creation assumed the strain prefix

The first audit configuration used display code `corF23` while the files were
named `cor_<accession>.pdb`. Analysis tools accepted the directory, but the HTML
download builder used the code as a filename prefix and silently produced empty
ZIPs.

The audit configuration was corrected to `cor`. The builder was also hardened:
it accepts an exact configured prefix, a bare accession filename, or one unique
alternative prefix. Multiple matching structures for one accession now cause a
clear ambiguity error.

Final ZIP inspection:

| Download | PDB files | Manifest |
|---|---:|---|
| `F0_member_structures.zip` | 10 complete parent structures | yes |
| `D0_member_structures.zip` | 5 cropped domain structures | yes |
| `D0_parent_structures.zip` | 5 complete parent structures | yes |
| `D1_member_structures.zip` | 2 cropped domain structures | yes |
| `D1_parent_structures.zip` | 2 complete parent structures | yes |

Every archived PDB was opened and checked for atomic-coordinate records.

## Test evidence

Local:

```text
pytest: 61 passed, 1 Python cgi deprecation warning
Python compileall: passed
JavaScript syntax check: passed
git diff --check: passed
```

4070:

```text
Final Snakemake atlas run: exit 0
Final ESM completion run: exit 0
Structure ZIP audit: PASS
Functional completeness audit: PASS
Independent/reference comparison: 29/29 PASS
```

Browser acceptance confirmed:

- v5.1.2 audit title and F/D navigation;
- separated D0 and D1 overview nodes;
- D0 segment statistics and parent F0 relationship;
- D0 ESM parent-context availability;
- RNA-seq heatmap and structural-similarity heatmaps;
- FoldTree and sequence-analysis status;
- full and domain download controls.

The in-app automation environment did not expose browser `fetch`, so it remained
at `Loading selected full proteins...`. The same structure artifact endpoint
returned HTTP 200 with the expected 77,355-byte TDZ20033.1 PDB. Interactive 3D
rendering in a normal user browser remains a manual acceptance item and is not
counted among the 29 computational comparisons.

## Logs and artifacts

On 4070:

```text
/home/claude/suss_audits/20260730-v5.1.1-F23/reports/SUSS_v5.1.2_F23_4070_EXECUTION.log
/home/claude/suss_audits/20260730-v5.1.1-F23/reports/SUSS_v5.1.2_F23_4070_AUDIT_BUNDLE.tar.gz
```

SHA-256:

```text
c37cf99cb9825928749cd41e9de651da91d22da2e3eae6219d2e2452503d672e  SUSS_v5.1.2_F23_4070_EXECUTION.log
4a283ba95d552d0a765400a6af350b34fb0bcaa5076fee670adcef2bf9f43721  SUSS_v5.1.2_F23_4070_AUDIT_BUNDLE.tar.gz
```

The consolidated log includes failed intermediate attempts. In particular,
`03_suss_atlas_full.log` records the family-path defect that stopped the first
run. Later logs record each fix and the successful final reruns. The bundle also
contains effective configuration, comparison JSON, key result tables, download
archives, and the reusable comparison script. It excludes the original PDB and
RNA workbook inputs.

## Claude acceptance checklist

1. Open `20260730-v5.1.2-F23-audit`.
2. Confirm the header reports 1 family, 0 singletons, and 10 proteins.
3. Open F0 and confirm ten members, RNA-seq, Foldseek and US-align heatmaps, and
   the three FoldTree views.
4. Confirm sequence conservation is explicitly `not applicable`, with the
   subgroup-size reason, rather than silently missing.
5. Download `F0_member_structures.zip`; confirm ten PDBs plus `manifest.tsv`.
6. Open Domain analysis; confirm D0 and D1 are separate nodes.
7. Open D0; confirm five segments, five parent proteins, F0 linkage, ESM
   availability, pocket evidence, RNA-seq, annotation, and independent D-family
   structural/sequence tabs.
8. Download both D0 structure ZIPs and confirm five cropped domains versus five
   complete parents.
9. Open D1 and repeat the structure/ESM/download check for two members.
10. Use a normal browser to verify the selected-member 3D viewer, white/black
    background, display styles, pockets, conservation coloring, and multi-member
    superposition.
11. Compare any disputed result with
    `controls/comparison_results.json` and the corresponding numbered execution
    log.

## Files most relevant to review

- `workflow/Snakefile`
- `workflow/scripts/v3_utils.py`
- `workflow/scripts/signature.py`
- `workflow/scripts/sasa_pocket.py`
- `workflow/scripts/esm_scan.py`
- `workflow/builders/html_builder.py`
- `workflow/builders/template/renderer.js`
- `validation/compare_f23_audit.py`
- `tests/test_v3_science.py`
- `tests/test_atlas_downloads.py`
