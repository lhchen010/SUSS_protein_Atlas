# SUSS Protein Atlas v5.0.2 handoff

This maintenance release corrects molecular-viewer hierarchy, maps Domain
structural conservation to every aligned segment, and reorganizes the atlas
navigation around the two independent clustering axes.

## Domain viewer hierarchy

The previous Domain viewer rendered the complete parent structure at about 22%
opacity. On a white background, non-domain and non-pocket regions were nearly
invisible.

v5.0.2 uses three explicit visual layers:

- complete parent context: solid light grey;
- matched D-family interval: the segment/family color;
- selected pocket-lining residues: red over a visible grey scaffold.

The hierarchy applies to cartoon, surface, stick, sphere, and line
representations on white or black backgrounds. Pocket mode deliberately greys
all non-pocket residues, matching the Full-family viewer semantics.

## Domain structural conservation

The Domain workbench already contained the FoldMason AA structural MSA and
official per-column lDDT values, but the viewer projected these values only onto
the D-family hub. Selecting a non-hub member therefore produced a grey domain
while still showing a red/blue legend.

v5.0.2 maps each scored FoldMason alignment column through every member's aligned
sequence:

1. Start at the retained segment's parent-protein residue coordinate.
2. Advance the parent coordinate for every non-gap alignment symbol.
3. Assign the shared column lDDT to that member residue.
4. Leave unscored alignment columns and the parent region outside the matched
   segment grey.

The calculation is derived from the existing structural MSA and does not rerun
FoldMason or alter D-family membership. In the reference D1 family, all 45
segments receive mapped FoldMason lDDT values.

## Full-family initial state

Full-family structures previously opened directly in Structural conservation
mode. Clicking the already-active button produced no visible state transition
and could be interpreted as a failed control.

Full families now open in an explicit **pLDDT** mode. Structural conservation is
applied only after the corresponding button is selected. The underlying
conservation files, score direction, and fixed 0-1 lDDT scale are unchanged.

All 105 Full families in the reference contain non-empty, non-constant
FoldMason structural-conservation scores.

## Two-axis navigation

The former flat navigation placed these four buttons at the same level:
Full-length families, Domain families, No D-family match, and Protein
singletons. That implied four mutually exclusive categories, which is
biologically incorrect.

v5.0.2 uses two levels:

- **Full-length analysis**
  - Families
  - Singletons
- **Domain analysis**
  - Families
  - No retained match

The last selected child view is remembered independently for each axis.

Reference cross-tabulation:

| Full-length state | Has D family | No D-family match |
|---|---:|---:|
| F-family member | 837 | 58 |
| Full-length singleton | 158 | 91 |

Totals:

- 895 proteins belong to 105 F families.
- 249 proteins are full-length singletons.
- 995 proteins belong to 172 D families.
- 149 proteins have no retained D-family match.

The 158 full-length singletons with a D-family assignment demonstrate why the
two axes must not be presented as one four-way partition.

## Reference validation

Reference run: `20260730-v5.0.2-reference`

- Existing F/D/S assignments and all scientific artifacts were reused without
  recomputation.
- D1 Domain mode displays the matched segment over visible complete-parent
  context.
- D1 fpocket mode displays red lining residues over a visible grey scaffold.
- D1 non-hub `TDZ19833.1` displays FoldMason structural-conservation colors.
- D1 reports 45/45 segments with mapped FoldMason lDDT.
- F43 opens in pLDDT mode and changes state when Structural conservation is
  selected.
- Full and Domain axis tabs expose only their corresponding child views.
- Switching between axes preserves the previous child selection.
- Python tests: 49 passed.
- Renderer JavaScript syntax check: passed.

## Claude acceptance checklist

1. Open `20260730-v5.0.2-reference` and confirm the first-level
   **Full-length analysis** and **Domain analysis** controls.
2. Under Full-length analysis, switch between **Families (105)** and
   **Singletons (249)**.
3. Under Domain analysis, switch between **Families (172)** and
   **No retained match (149)**.
4. Leave Full-length analysis on Singletons, switch to Domain, then return to
   Full-length analysis and confirm Singletons remains selected.
5. Open D1 and select `TDZ19833.1`.
6. In Domain mode, confirm the retained segment is colored and the remaining
   parent protein is visible in light grey.
7. Select Pocket and confirm non-pocket structure remains visible while pocket
   residues are red.
8. Select Structural conservation and confirm the non-hub segment changes to
   the red/blue FoldMason scale.
9. Open Conservation + Pockets and confirm D1 reports 45 segments with
   FoldMason lDDT.
10. Open F43 and confirm pLDDT is the initial mode.
11. Select Structural conservation and confirm the button and molecular
    coloring change to the FoldMason conservation state.
12. Repeat the viewer checks with black background and surface representation.

## Intentional semantics

- Full-length singleton means no retained Full-family assignment.
- No retained D-family match means no retained local D-family assignment.
- These states belong to different axes and may overlap.
- Domain structural conservation describes shared alignment columns and is
  mapped only within each retained segment; residues outside the segment remain
  parent context.
