# SUSS Protein Atlas v4.1.0 handoff for Claude for Science

## Purpose

Version 4.1 corrects pocket-reference selection, applies the appropriate P2Rank
profile to predicted structures, restores full-length family RNA-seq plots, and
replaces the all-or-nothing domain superposition with a member-selectable viewer.
It does not change F-family, D-family, S-subgroup, or SUSS classification
thresholds.

## Pocket correctness

### Root cause

Version 4.0 changed the canonical hub of 22 F families. The P2Rank output
directory retained files for both the previous and current hub. Two code paths
selected `glob(...)[0]`:

- `workflow/scripts/sasa_pocket.py`;
- `workflow/builders/html_builder.py`.

The arbitrary first result caused nine F families (`F0`, `F2`, `F3`, `F4`,
`F6`, `F12`, `F13`, `F19`, and `F26`) to display or export a previous hub's
P2Rank prediction. fpocket used an accession-specific directory and was not
affected. A raw-output audit found exact fpocket agreement for all 354 targets.

### v4.1 behavior

For every P2Rank target, the workflow now:

1. removes that target's previous `out/` directory;
2. runs P2Rank in a fresh directory;
3. requires exactly one prediction CSV whose filename ends with the current
   reference accession;
4. fails the rule if the expected current-reference output is missing or
   ambiguous.

The HTML builder uses the same reference-specific selection rule when enriching
legacy results and when exporting the detector-native P2Rank table. A lone
legacy prediction file remains readable, but a directory with multiple files
and no current-reference match is not guessed.

### AlphaFold profile

The new configuration field is:

```yaml
pocket:
  p2rank_profile: alphafold
```

`alphafold` is the default because the atlas input contract expects predicted
structures whose PDB B-factor column contains pLDDT. The workflow passes
`-c alphafold` to P2Rank. Users processing experimental structures can select
`default` in the portal advanced options or configuration.

The profile is recorded in `used_config.yaml`, `pockets.json` for new runs, and
the family workbook pocket summary. P2Rank `score` and calibrated
`probability` remain separate columns. Changing the profile changes the
scientific prediction, so v4.0 and v4.1 P2Rank scores must not be treated as a
software-only presentation difference.

## Full-length family RNA-seq

### Root cause

The builder loaded the valid run-level `results/rnaseq_expression.csv`, but
full-length families still looked only for obsolete files named
`results/families/<F>/<F>_expression.csv`. The renderer then created
`<img src="undefined">`.

### v4.1 behavior

The builder first supports a legacy per-family table when one exists, then
subsets the run-level expression table by the current family's members. The
result is used for:

- the RNA-seq heatmap in the full-length family workbench;
- the `RNAseq` sheet in the family workbook.

If no family expression exists, the renderer shows an explicit unavailable
message and never creates an image with an undefined source.

## Selectable domain-family superposition

The D-family workbench now defaults to the focused segment plus the D-family
hub. The user can select domains with:

- one checkbox per segment;
- `All`;
- `None`;
- `Neighbors`, which selects the focused segment and its retained local
  Foldseek neighbors.

`Superpose selected (N)` is enabled only when at least two segments are
selected. The selected cropped domains are transformed to the deterministic
D-family hub using the existing US-align transforms. This changes the
interaction and display only; D-family membership and transforms are not
recomputed by the browser.

Each selected segment receives a stable color and a legend. The hub is always
gold. The focused full-length parent can be added as a translucent gray
`Parent context`; this option is off by default so a large parent does not
shrink the aligned-domain view.

Every member also has a linear parent-protein map and an explicit label such as:

```text
TDZ22296.1: domain 376-592 / 592 aa (36.7%, C-terminal)
```

This distinguishes an embedded local domain from a whole-protein match. The
selected-domain superposition can be downloaded as a multi-model PDB.

## Portal and configuration

The build portal adds an advanced P2Rank input-structure selector:

- `AlphaFold / predicted structures` maps to `alphafold` and is the default;
- `Experimental PDB structures` maps to `default`.

The choice is included in the run manifest display and applied to the generated
effective configuration. Existing analysis toggles and clustering thresholds
are unchanged.

## Data-contract changes

New configuration:

```text
pocket.p2rank_profile
```

New or expanded pocket workbook fields:

```text
pocket_summary.profile
pocket_summary.top_probability
pocket_predictions.probability
```

The domain viewer adds browser state for selected segment IDs and parent
context. `domain_workbench.json` and `domain_members.csv` do not change schema.

## Acceptance checklist for Claude

1. Run the full Python test suite and compile `portal/suss_portal.py`,
   `html_builder.py`, and `sasa_pocket.py`.
2. Run a JavaScript syntax check on `renderer.js`.
3. Confirm `config/config.yaml.template` defaults
   `pocket.p2rank_profile` to `alphafold`.
4. Confirm a new P2Rank command includes `-c alphafold`.
5. Place stale and current-reference P2Rank CSV files in one test directory and
   verify only the current-reference table is parsed and exported.
6. Audit every `pockets.json` key against its `ref` and confirm the raw P2Rank
   output filename contains that accession.
7. Open F4 and confirm its reference is `TDZ22880.1`.
8. Compare the F4 P2Rank residues with
   `results/p2rank/F4/out/*TDZ22880.1*_predictions.csv`.
9. Compare the F4 fpocket residues with
   `results/fpocket/F4/TDZ22880.1_out/pockets/pocket1_atm.pdb`.
10. Open a full-length family with RNA-seq and confirm the image has nonzero
    natural width and height and no `src="undefined"`.
11. Download that family's workbook and confirm the `RNAseq` sheet contains
    only family members.
12. Open D35 with `TDZ22296.1:376-592` focused.
13. Confirm the parent map marks residues 376-592 of 592 and labels the segment
    C-terminal.
14. Confirm the initial superposition selection is the focused segment plus the
    gold hub.
15. Confirm `All`, `None`, and `Neighbors` update the checkbox count and
    superposition button.
16. Superpose two or more selected segments and confirm only those colored
    segments are displayed.
17. Enable `Parent context` and confirm only the focused full-length parent is
    added in translucent gray.
18. Download the selected superposition and confirm it is a valid multi-model
    PDB.
19. Submit a small portal build and confirm the selected P2Rank profile appears
    in the run parameters and `used_config.yaml`.

## Validated 4070 reference

Validation completed on 2026-07-28 against:

```text
http://100.80.77.29:8600/atlas?id=20260728-v4.1-reference
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
| Pocket targets | 354 |
| P2Rank current-reference files | 354 / 354 |
| P2Rank parsed values matching raw tables | 354 / 354 |
| P2Rank profile | `alphafold`, 354 / 354 |
| fpocket parsed values matching raw outputs | 354 / 354 |

F4 uses `TDZ22880.1`. Its v4.1 P2Rank top pocket has score `3.12`,
probability `0.084`, and residues:

```text
79, 81, 111, 116, 118, 121, 200, 201, 205
```

Its fpocket top pocket has score `0.314` and residues:

```text
81, 111, 116, 118, 200, 201, 202, 203, 204, 205
```

Seven residues are shared by the two independently generated top predictions.
This confirms reference mapping and detector agreement for this acceptance
case; it is not experimental proof that every predicted cavity binds a ligand.

The F4 workbook has 23 sheets. Its `RNAseq` sheet contains 43 unique F4
members and the five conditions `Co_C`, `Co_VH`, `Co_1DPI`, `Co_3DPI`, and
`Co_7DPI`. The pocket summary records profile `alphafold`, P2Rank score `3.12`,
and probability `0.084`.

Browser validation confirmed:

- the F4 RNA-seq SVG has natural size `470 x 804` and a valid data URI;
- F4 shows the current fpocket and P2Rank values above;
- D35 opens `TDZ22296.1:376-592` with the correct `592 aa`, `36.7%`, and
  `C-terminal` label;
- D35 initially selects the focused segment and gold hub (`2` members);
- `All`, `None`, and `Neighbors` produced selected counts `7`, `0`, and `7`;
- the superposition button is disabled at zero selected members;
- selected gold and magenta domains are visibly aligned without loading the
  other five segments;
- `Parent context` adds the focused full-length parent and is off by default;
- the portal exposes both AlphaFold/predicted and experimental-PDB P2Rank
  choices.

Generated atlas:

```text
results/cor_suss_atlas.html
size:   106,912,326 bytes
sha256: 6f72dcbd919e6695b0b459474a2fe3362fc30503151a77aa76bd7d068f1a10fa
```

`results/downloads/` contains 459 files. Portal streaming checks returned HTTP
`200` for the atlas, F4 workbook (`170,151` bytes), and F4 member-structure ZIP
(`1,334,764` bytes).

Automated validation:

```text
pytest:                       44 passed
extracted engine pytest:      44 passed
Python compile check:         passed
renderer.js syntax check:     passed
354-target pocket audit:      passed with zero errors
Snakemake final dry run:      nothing to be done
Portal and artifact routes:   HTTP 200
```

Known non-blocking warnings:

- Python 3.11 reports that the standard-library `cgi` module is deprecated for
  Python 3.13.
- Four baseline outputs created before complete Snakemake provenance tracking
  (`foldseek`, `qc`, `seqs`, and `validate`) report missing historical
  metadata; the final dry run reports no pending jobs.

Rollback copies on 4070:

```text
/home/claude/suss_atlas_staging/20260728-v3/backups/v4.1-pre-code-20260728.tar.gz
/home/claude/suss_atlas_staging/20260728-v3/backups/pockets.v4.1-pre-profile.json
/home/claude/suss_portal/suss_engine.tar.gz.pre-v4.1.0
/home/claude/suss_portal/suss_portal.py.pre-v4.1.0
```
