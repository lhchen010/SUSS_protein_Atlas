# Example — expected output

The bundled example is **100 secreted proteins of *Colletotrichum orbiculare*** (validated
"Set A" subset of the published atlas). It is a **core-only** run: clustering + independent-TM
cross-check + conservation + pocket + FoldTree + cards + atlas. It needs only the conda
environment plus US-align and FoldTree — **no external databases, no InterProScan, no GPU** — so
it doubles as a "does my install work?" self-check.

## Run it

```bash
conda activate suss-atlas
# from the repo root:
snakemake --configfile examples/config.example.yaml --cores 8 results/example_suss_atlas.html
```

(Edit `tools.fold_tree` and `tools.usalign` in `examples/config.example.yaml` to your install
paths first. To skip those two, set `steps.foldtree: false` and `steps.tm_refine: false` — the
rest still runs on the conda env alone.)

## Expected result (engine v1.0.1-rc1)

| Metric | Expected value |
|---|---|
| Total members | 100 |
| Clustered families | 6 |
| Singletons | 12 |
| Within-family + cross edges (classification rows) | 536 |
| Atlas HTML size | ~13.5 MB |
| US-align cross-check matrices | 6 (one per family) |
| FoldTree trees | 12 files (3 metrics for each family with n >= 3) |

Top families (renumbered by size each run; **F0 = largest**):

| Family | n | mean Foldseek TM | core-SUSS % |
|---|---|---|---|
| F0 | 28 | 0.594 | 92.8 |
| F1 | 26 | 0.669 | 13.7 |
| F2 | 23 | 0.603 | 79.0 |
| F3 | 7  | 0.605 | 54.5 |
| F4 | 2  | 0.642 | 100.0 |
| F5 | 2  | 0.907 | 0.0 |

Small numeric differences (±1 in family sizes, a few % in TM/SUSS) can occur across tool
versions and platforms — the family **count** and the F0/F1/F2 sizes (28/26/23) are the stable
signal that your install is working. If you get 6 families with F0≈28, you're good.

Open `results/example_suss_atlas.html` in a browser (Chrome/Firefox): click any node to see the
integrated card — structure viewer, dual TM matrices (Foldseek + US-align with the consistency r),
sequence-identity matrix, and the FoldTree.

F4 and F5 each contain only two structures. FoldTree is intentionally skipped for these families:
its sub-workflow cannot infer a meaningful tree from two members, and an absent tree is no longer
represented by a misleading zero-byte Newick file.
