# SUSS Protein Atlas v4.0.0 handoff for Claude for Science

## Purpose

Version 4 turns local structural-domain families into complete analysis workbenches, corrects
several metric labels and reference-selection rules, and introduces a server-backed data layout
for large atlases. The existing scientific hierarchy is preserved:

- `F`: full-length structure-defined families;
- `D`: local structure-defined domain families;
- `S`: sequence-homologous subgroups inside an F family;
- singleton: an independent protein without a retained full-length family.

The D-family workbench, MAFFT alignment, tree, and US-align superposition are downstream
descriptions. They do not redefine D-family membership.

## Scientific definitions

### Full-length F family

An undirected structural graph edge is retained when:

```text
min(query TM, target TM) >= clustering.foldseek_tm
AND query coverage >= clustering.whole_fold_min_coverage
AND target coverage >= clustering.whole_fold_min_coverage
```

The default thresholds are TM `0.5` and reciprocal coverage `0.5`. Leiden communities define
the F families.

### Local D family

The local Foldseek search uses alignment type `2` and 3Di+AA evidence. The default retained link
requires:

```text
e-value <= 1e-3
AND probability >= 0.5
AND aligned residues >= 40
AND local lDDT >= 0.5
```

`domain_clustering.min_alntm` is an optional second geometry threshold and defaults to `0.0`.
`domain_clustering.min_shorter_coverage` also defaults to `0.0`, so an embedded domain can be
found even when the full proteins have different architectures. Observed query and target
coverage are still recorded. Leiden edge weights are `probability * local lDDT`.

### Sequence S subgroup and SUSS

BLAST plus reciprocal sequence coverage defines S subgroups and labels structural links. It does
not create or split structural families. MAFFT and Rate4Site are run only for an eligible
sequence-homologous subgroup. FoldMason remains the structural AA/3Di alignment.

## Domain-family workbench

`workflow/scripts/domain_workbench.py` adds the following per D family:

- a deterministic hub selected by summed retained local lDDT;
- mature sequences cropped to the D-segment coordinates;
- MAFFT segment sequence MSA;
- FastTree sequence relationship tree;
- US-align transforms from each segment to the D-family hub;
- per-segment TM and RMSD fit statistics;
- explicit component status values.

The atlas presents the same workbench pattern used for full-length families:

1. **Structure + Network**: full protein with the selected segment highlighted, a stabilized
   segment network, and hub-referenced domain superposition.
2. **Struct sim**: every retained Foldseek local edge with lDDT, alignment TM, probability,
   aligned length, and query/target coverage.
3. **Sequence + MSA**: cropped segment sequences, MAFFT MSA, and sequence tree.
4. **RNA-seq**: expression inherited from each segment's parent protein.
5. **Annotation**: coordinate-overlapping Pfam/InterPro calls separated from parent-protein
   evidence, pockets, EffectorP, DeepTMHMM, and Foldseek database hits.

The domain label and `annotation:`/`domain:` search use only annotations overlapping the segment
coordinates. `gene:`, `acc:`, `accession:`, and `protein:` match accessions, gene names, and
segment IDs. `evidence:` searches the broader parent-protein evidence.

### NLP/NPP regression

All seven true NPP1/NLP proteins are recovered in `D35`:

| Protein | Segment |
|---|---|
| `TDZ17504.1` | `1-263` |
| `TDZ17577.1` | `2-244` |
| `TDZ21662.1` | `6-230` |
| `TDZ24640.1` | `2-240` |
| `TDZ25257.1` | `6-218` |
| `TDZ22296.1` | `376-592` |
| `TDZ25193.1` | `8-244` |

The `TDZ22296.1` C-terminal segment is the important domain-architecture regression: it was
previously missed or represented by the wrong interval.

## Metric and reference corrections

### BLAST relationship selection

`classify` and `sequence_subgroups` now share one relationship selector. It first applies both
the configured e-value and reciprocal-coverage thresholds, then selects the best qualifying HSP
for each undirected pair. A very significant but short HSP can no longer hide a second HSP that
actually passes the relationship definition.

### Explicit denominators

Ambiguous `mean_TM` and `mean_identity` labels are replaced in user-facing tables by:

- `foldseek_TM_all_pairs`: all unique within-family pairs, with undetected pairs represented by
  the matrix value zero;
- `foldseek_TM_detected_pairs`: reported Foldseek pairs only;
- `mean_retained_edge_TM`: graph edges passing clustering thresholds;
- `blast_identity_all_pairs`: all unique pairs, including BLAST-undetected zero entries;
- `blast_identity_detected_pairs`: BLAST-detected pairs only;
- `max_blast_identity`;
- `mean_retained_edge_foldseek_fident`: Foldseek alignment identity on retained structural
  edges, explicitly not BLAST identity.

Singleton pairwise metrics are null with pair count zero. They are not assigned a synthetic TM
value of `1.0`.

### Canonical family hub

The canonical hub now maximizes mean Foldseek TM across every possible within-family partner.
Available below-threshold Foldseek pairs are included and missing pairs contribute zero, matching
the family matrix. Ties are deterministic. This replaces selection based only on retained graph
edges.

The reference change affects 22 families in the *C. orbiculare* reference:

| Family | Previous | v4 hub |
|---|---|---|
| F0 | TDZ18953.1 | TDZ18834.1 |
| F1 | TDZ20858.1 | TDZ15115.1 |
| F2 | TDZ26256.1 | TDZ20663.1 |
| F3 | TDZ18005.1 | TDZ21492.1 |
| F4 | TDZ17686.1 | TDZ22880.1 |
| F5 | TDZ17136.1 | TDZ15166.1 |
| F6 | TDZ14865.1 | TDZ19452.1 |
| F8 | TDZ20849.1 | TDZ19813.1 |
| F9 | TDZ20221.1 | TDZ15165.1 |
| F11 | TDZ16193.1 | TDZ19361.1 |
| F12 | TDZ17640.1 | TDZ20059.1 |
| F13 | TDZ24888.1 | TDZ27012.1 |
| F14 | TDZ15636.1 | TDZ19840.1 |
| F16 | TDZ22883.1 | TDZ14903.1 |
| F19 | TDZ21289.1 | TDZ16976.1 |
| F26 | TDZ14099.1 | TDZ20573.1 |
| F28 | TDZ25288.1 | TDZ20177.1 |
| F35 | TDZ14532.1 | TDZ24120.1 |
| F37 | TDZ17244.1 | TDZ22709.1 |
| F39 | TDZ26080.1 | TDZ16998.1 |
| F49 | TDZ19242.1 | TDZ14381.1 |
| F51 | TDZ25143.1 | TDZ18955.1 |

### P2Rank

P2Rank's ranking `score` and calibrated `probability` are stored separately. Full-family and
singleton atlas labels use `prob` with `top_probability`; Excel retains both fields. fpocket
continues to display its own score.

### Structural conservation color direction

FoldMason column LDDT is stored in the PDB B-factor field on a 0-100 scale. The renderer is
explicitly:

- red: high structural conservation;
- blue: low structural conservation;
- light gray: insufficient FoldMason pair support.

pLDDT red remains higher confidence. ESM red remains greater mutation tolerance, which is a
different biological quantity.

## Server-backed data layout

Portal jobs now default to `output.html_mode: backend`.

- Source PDB files are served lazily through `/artifact?kind=structure`.
- Family and singleton workbooks are stored in `results/downloads/`.
- Family structure bundles remain ZIP files with one original PDB per member plus a manifest.
- The portal streams the atlas and artifacts instead of reading each large file fully into RAM.
- The HTML fetches a structure only when a singleton, domain segment, or full-family
  superposition needs it.
- `output.html_mode: single` remains available for an offline self-contained atlas.

The structure-search result page now links each accession, F family, singleton, and D family
directly to its matching atlas workbench through `protein=` or `open=` query parameters.

## Data-contract changes

New declared output:

```text
results/domain_workbench.json
```

New backend artifact directory:

```text
results/downloads/<F-family>_data.xlsx
results/downloads/<F-family>_member_structures.zip
results/downloads/<singleton>_data.xlsx
```

The atlas payload adds:

```text
DOMAIN_WORKBENCH
BACKEND.enabled
```

`domain_families.csv` adds mean alignment TM and query/target coverage. `domain_edges.csv`
retains alignment TM and both coverage values. Family summary files add the explicit denominator
fields described above and separate P2Rank score/probability.

## Acceptance checklist for Claude

1. Run the Python tests, Python compile checks, JavaScript syntax check, and Snakemake dry run.
2. Confirm all seven NLP/NPP proteins above appear in `D35`, including
   `TDZ22296.1:376-592`.
3. Search `annotation:nlp`; confirm a unique D-family result opens the matching segment rather
   than an unrelated first member.
4. Open D35 and confirm its segment network stops moving after stabilization.
5. Select a segment and confirm the full parent structure is gray, the selected coordinates are
   teal, and pocket residues are red.
6. Turn on D-family superposition and confirm all displayed segments align to the gold hub using
   US-align transforms.
7. Confirm the Struct sim tab reports lDDT, alignment TM, probability, aligned residues, and
   query/target coverage.
8. Confirm the Sequence + MSA tab contains cropped D-segment sequences, MAFFT MSA, and a
   sequence tree.
9. Confirm RNA-seq and Annotation use the same tabbed presentation pattern as full-length
   families.
10. Open F4 and verify all-pair, detected-pair, and retained-edge labels are distinct.
11. Confirm the BLAST qualifying-hit regression passes when the lowest-e-value HSP has inadequate
    coverage but another HSP passes both thresholds.
12. Confirm F4 uses `TDZ22880.1` as its canonical hub.
13. Confirm P2Rank displays probability rather than its ranking score.
14. Confirm structural conservation is red for high LDDT, blue for low LDDT, and gray for
   insufficient support.
15. Open a singleton in backend mode and confirm its structure loads on demand.
16. Download a family workbook and member-structure ZIP; inspect that the ZIP contains multiple
   PDB files and a manifest.
17. Submit a PDB structure search and follow the accession, F-family, and D-family links back
   into the atlas.

## Validated 4070 reference

Validation completed on 2026-07-28 against:

```text
http://100.80.77.29:8600/atlas?id=20260728-v4-reference
```

Reference dataset:

| Check | Validated result |
|---|---:|
| Proteins | 1,144 |
| Full-length F families | 105 |
| Clustered proteins | 895 |
| Singletons | 249 |
| Domain D families | 172 |
| Domain segments | 1,014 |
| Retained domain edges | 3,883 |
| Complete D workbenches | 172 / 172 |
| US-align domain transforms | 1,014 |
| Sequence S subgroups | 302 |
| Complete annotations | 1,144 / 1,144 |
| Complete P2Rank + fpocket targets | 354 / 354 |

The 5,050 retained full-length structural edges classify as:

```text
core_SUSS          2,958
moderate_paralog   1,151
diverged_paralog     740
recent_duplicate     201
```

There are 2,092 sequence-detected structural relationships. For F4, the atlas reports mean
Foldseek TM `0.483` across all 903 pairs, mean retained-edge TM `0.65`, mean BLAST best-HSP
identity `12.7%` across all pairs, `36.3%` across 315 reported pairs, and maximum BLAST identity
`73.8%`. The canonical F4 hub is `TDZ22880.1`.

Generated artifacts:

```text
results/cor_suss_atlas.html
size:   105,406,092 bytes (101 MiB)
sha256: b339566fc7d794784f4d0775e17d9c2098822ed989cad6f652e780aa6822b45e
```

`results/downloads/` contains 354 workbooks and 105 multi-PDB family ZIP archives, 459 files in
total. The HTML is about 78% smaller than the previous approximately 453 MiB self-contained
reference because original structures and download bundles are served lazily.

Automated validation:

```text
pytest:                         41 passed
Python compile check:           passed
renderer.js syntax check:       passed
databridge.js syntax check:     passed
Snakemake final dry run:        nothing to be done
Portal health / atlas / ZIP:    HTTP 200
```

Browser validation confirmed F4 metric labels and hub, the structural-conservation legend,
D35's stabilized network and nonblank structure viewer, US-align family superposition, all five
D-family tabs, `annotation:nlp`, `gene:TDZ22296.1`, lazy singleton structures, and singleton
P2Rank probability. A structure query with `TDZ22880.1.pdb` returned 39 Foldseek hits and direct
links to the accession, F4, and D0 workbenches. The query directory was removed after completion.

Known non-blocking warnings:

- Python 3.11 reports that the standard-library `cgi` module is deprecated for Python 3.13.
- Four baseline files created before Snakemake provenance tracking (`foldseek`, `qc`, `seqs`,
  and `validate`) report missing historical metadata in dry-run output; all requested outputs
  are present and current.

Pre-v4 rollback copies on 4070:

```text
/home/claude/suss_portal/suss_engine.tar.gz.pre-20260728-v4.0.0
/home/claude/suss_portal/suss_portal.py.pre-20260728-v4.0.0
```
