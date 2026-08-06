# SUSS Protein Atlas v5.2.0 UI handoff

## Decision

Release candidate status: **UI acceptance pending; automated checks pass**.

v5.2.0 modernizes the portal and atlas presentation without changing the
scientific workflow, clustering definitions, thresholds, output schemas, or
analysis data. The visual direction is a compact scientific console: dark
application chrome, serif hierarchy for project and family titles, monospace
labels for controls and measurements, cyan full-length accents, amber domain
accents, and light analysis canvases where structures, trees, heatmaps, and
network labels require high contrast.

## Changed surfaces

| Surface | Change |
|---|---|
| Portal | Branded top bar, clearer four-step build form, compact analysis controls, persistent run-state sidebar, and consistent status/download actions |
| Full-length atlas | Dark analysis chrome and detail panel, light network/viewer canvases, compact tabs, and cyan active states |
| Domain atlas | Same workbench hierarchy as full-length analysis, with amber mode identity and shared control styling |
| Singleton view | Retains a light, scan-friendly data table inside the same application shell |
| Responsive layout | Portal stacks at 940 px; atlas converts to a vertical main/detail layout at 700 px |

## Deliberate non-changes

- No Foldseek, BLAST, FoldMason, FoldTree, pocket, ESM, Rate4Site, or RNA-seq
  computation changed.
- No family or singleton membership changed.
- No download endpoint, artifact format, or atlas JSON schema changed.
- Structure viewers and quantitative figures remain light-backed for legibility.

## v5.2.1 follow-up

The full-length member-download panel and the equivalent singleton download
panel now inherit the dark detail-workbench theme. Their legacy inline light
background was removed; download behavior and artifact contents are unchanged.

## Verification

```text
pytest -q
66 passed

python3 -m py_compile portal/suss_portal.py workflow/builders/html_builder.py
PASS

git diff --check
PASS
```

The F23 audit atlas was rebuilt on the 4070 host from the existing validated
result bundle using only the `assemble` rule. This isolates the visual change
from scientific recomputation:

```text
/home/claude/data/suss_ui_v52_preview_20260807/results/f23_audit_atlas.html
```

## Claude acceptance checklist

1. Open the portal at `http://100.80.77.29:8600/` and confirm the build form,
   recent runs, help, status, and history pages use one coherent visual system.
2. Open a full-length family and confirm the network remains light, detail data
   remain readable, tabs scroll instead of compressing, and viewer background
   controls still work.
3. Open a domain family and confirm amber mode cues distinguish it from the
   cyan full-length context without changing the available analyses.
4. Open singletons and confirm the table remains dense, sortable, and readable.
5. Check a narrow browser window: the portal sidebar must stack and the atlas
   main/detail views must become vertical without overlapping controls.
6. Exercise search, downloads, structure style/background buttons, and
   full-length/domain navigation to confirm all existing actions remain live.

## Acceptance boundary

Visual preference remains a human-review item. Automated tests protect the
shared tokens and layout contract, but approval should be based on the actual
portal and atlas in a browser. Scientific regression is not expected because
this release does not alter analysis code or data contracts.
