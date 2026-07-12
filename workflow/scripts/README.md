# workflow/scripts/ — Rule Script Inventory

Each rule corresponds to one clean script (consumes snakemake.input/output/params, zero absolute paths).
Source = analysis code recovered from artifact lineage (see recovery observations R1–R6 in ../../pipeline_io_contract.md).

| Script | Rule | Recovery Source | Status |
|---|---|---|---|
| qc.py | 1 qc | recovered/qc.py (removed Drive downloads, kept parsing + filtering) | Pending Refactor |
| cluster.py | 3 cluster | recovered/network_families.py | Pending Refactor |
| classify.py | 4 classify | recovered/classification.py | Pending Refactor |
| sasa_pocket.py | 7 sasa_pocket | recovered/sasa_all.py + fpocket_all.py + sasa_pockets.py | Pending Refactor |
| esm_scan.py | 8 esm | recovered/esm_driver.py | Pending Refactor |
| signature.py | signature | recovered/signature.py | Pending Refactor |
| cards.py | cards | recovered/master.py (card rendering logic) | Pending Refactor |
| annotate.py | 10 annotate | recovered/annotation.py | Pending Refactor |
| superpose.py | assemble-used | recovered/superpose_bb.py (backbone superposition, not CA-only) | Pending Refactor |

## Rules that directly invoke external tools (no Python script, shell call only)
- rule 2 foldseek → `foldseek easy-search` (Snakefile shell)
- rule 5 msa → `foldmason` (shell)
- rule 6 conservation → `rate4site` (shell)
- rule 9 foldtree → official FoldTree snakemake sub-pipeline (see design document Step5 recipe: targets before --config, base env end-to-end)

## 4070 measured tool paths (Literals, validated at dry-run)
foldseek=/home/user/miniconda3/bin/foldseek · foldmason=/home/user/miniconda3/bin/foldmason
rate4site=/home/user/.local/bin/rate4site · p2rank=/home/user/p2rank_app/p2rank/prank (requires java17 conda env)
esm_scan=/home/user/ESM-Scan/esmscan.py · fold_tree=/home/user/fold_tree
interproscan=/opt/interproscan/interproscan.sh · effectorp=/home/user/.local/share/effectorp/EffectorP.py
deeptmhmm=/home/user/.local/share/deeptmhmm/…/predict.py · fpocket (local struct-tools) · freesasa (py)
