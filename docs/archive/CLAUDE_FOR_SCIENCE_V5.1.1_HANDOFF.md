# SUSS Protein Atlas v5.1.1 handoff

## Purpose

This patch fixes a misleading interpretation in the Domain analysis layer. A protein outside a
retained `D` family was previously shown only as `No retained match`, even when Foldseek had
reported a strong local hit that narrowly missed one configured filter.

The patch adds an auditable diagnostic layer. It does **not** change any Foldseek threshold,
accepted edge, Leiden community, `F` family, or `D` family.

## Scientific behavior

Every protein now receives one of four Domain diagnostic states:

| State | Meaning |
|---|---|
| `retained_family` | At least one segment is in a retained D family |
| `borderline` | The best local hit passes every filter except local lDDT and is within the configured diagnostic margin |
| `filtered_hit` | Foldseek reported a non-self local hit, but the best candidate failed one or more filters |
| `no_raw_hit` | Foldseek reported no non-self local hit under the search settings |

An additional `unclustered_retained_hit` state is supported for a hit that passes all edge
filters but does not survive community/minimum-family-size selection.

`domain_clustering.borderline_lddt_margin` defaults to `0.05`. It is diagnostic only:

```text
borderline =
  passes e-value
  AND passes probability
  AND passes aligned-residue count
  AND passes shorter-protein coverage
  AND passes optional alignment-TM threshold
  AND local lDDT < min_lddt
  AND local lDDT >= min_lddt - borderline_lddt_margin
```

Changing this margin changes only the explanatory label. It never retains an edge or assigns a
protein to a D family.

## New artifact

`results/domain_match_diagnostics.csv` contains one row per protein:

- best local Foldseek match;
- query-oriented and match-oriented coordinates;
- e-value, probability, bit score, alignment TM, local lDDT, sequence identity;
- query/match/shorter coverage and aligned length;
- raw directed-hit count;
- number of filters passed;
- exact failed-filter text;
- diagnostic status and summary.

The ranking prefers a fully passing hit, then a diagnostic-borderline hit, then the candidate
passing the most filters, followed by lDDT, alignment TM, probability, alignment length, and
e-value. Incoming and outgoing directed Foldseek records are reoriented to the protein being
reported, so the displayed coordinates and coverages remain unambiguous.

## Atlas changes

The Domain tab formerly called `No retained match` is now `Unassigned / filtered`.

The table reports:

- protein and full-length family context;
- diagnostic state;
- best local match and both aligned ranges;
- local lDDT, alignment TM, probability, and both observed coverages;
- exact failed threshold;
- annotation.

The table has a `Download diagnostic CSV` action. It opens at full width; selecting a protein
restores the detail panel and opens the existing full-protein/F-family workbench.

Search includes accession, annotation, diagnostic status, best match, filter reason, scores, and
expression fields. Searching for an accession can therefore show both that protein and another
unassigned protein whose best match is the searched accession.

## Reference result

Deployed reference:

```text
20260730-v5.1.1-reference
```

Portal:

```text
http://100.80.77.29:8600/atlas?id=20260730-v5.1.1-reference
```

Reference counts:

| Diagnostic state | Proteins |
|---|---:|
| `retained_family` | 995 |
| `borderline` | 28 |
| `filtered_hit` | 57 |
| `no_raw_hit` | 64 |
| Domain-unassigned total | 149 |

The existing biological results remain:

| Result | v5.1.1 |
|---|---:|
| Input proteins | 1,144 |
| F families | 105 |
| Full-length singletons | 249 |
| D families | 172 |
| Retained D segments | 1,014 |
| Cross-D-family bridges | 6 |

## TDZ19858.1 acceptance case

Expected result:

```text
Full-length context: F23
Domain state: Borderline local hit
Best local match: TDZ23942.1
TDZ19858.1 range: 1-161
TDZ23942.1 range: 8-175
local lDDT: 0.4974
alignment TM: 0.5607
probability: 1.0
observed coverage: 0.982 / 0.960
failed filter: local lDDT 0.497 < 0.5
```

This is the intended interpretation: TDZ19858.1 remains outside a retained D family because the
configured local-lDDT cutoff is `0.5`, but the atlas no longer implies that Foldseek found no
local structural evidence.

## Validation completed

- Local test suite: `56 passed`.
- Python compile check: passed.
- Whitespace/error check: passed.
- Snakemake executed the new diagnostic rule on the 4070 reference dataset.
- Generated CSV counts sum to all 1,144 proteins.
- Existing F/D family counts and memberships were not recomputed or changed.
- Backend HTML rebuilt from the existing validated v5.1.0 result set plus the new diagnostic
  artifact.
- Browser acceptance verified the Domain tab, state counts, search, TDZ19858.1 values,
  full-width layout, and transition back to the F23 workbench.

## Claude acceptance checklist

1. Open `20260730-v5.1.1-reference`.
2. Select `Domain analysis` and `Unassigned / filtered`.
3. Confirm counts `28 borderline`, `57 filtered / unclustered`, and `64 no raw hit`.
4. Search `TDZ19858.1`.
5. Confirm the values listed in the TDZ19858.1 acceptance case.
6. Confirm the text states that the diagnostic label does not change D-family membership.
7. Download `domain_unassigned_match_diagnostics.csv` and compare TDZ19858.1 with
   `results/domain_match_diagnostics.csv`.
8. Select TDZ19858.1 and confirm the existing F23 full-length workbench opens.
9. Regression-check one retained D family and one full-length family to confirm their previous
   workbenches remain unchanged.

## Files most relevant to review

- `workflow/scripts/v3_utils.py`: diagnostic classification and ranking.
- `workflow/scripts/domain_diagnostics.py`: Snakemake artifact writer.
- `workflow/Snakefile`: rule and assemble dependency.
- `workflow/builders/html_builder.py`: CSV-to-atlas data contract.
- `workflow/builders/template/renderer.js`: Domain unassigned table and download.
- `docs/FOLDSEEK_PARAMETERS.md`: parameter semantics.
- `tests/test_v3_science.py`: borderline scientific unit test.
- `tests/test_atlas_downloads.py`: renderer/workflow contract assertions.
