"""Run FoldTree with validated, per-metric rooting recovery."""
import json
import os, glob, re
from pathlib import Path
from foldtree_runner import run_foldtree_family
from runtime_utils import resolve_executable
fam      = snakemake.wildcards.fam
famfile  = snakemake.input.famfile
pdb_dir  = snakemake.input.pdb_dir
outputs  = list(snakemake.output.trees)
status_output = str(snakemake.output.status)
ftdir    = snakemake.params.ftdir
smk      = resolve_executable(snakemake.params.snakemake, "FoldTree Snakemake")
foldseek = resolve_executable(snakemake.params.foldseek, "FoldTree Foldseek")
extra_path = str(snakemake.params.get("extra_path", "") or "")
metrics  = [str(metric) for metric in snakemake.params.metrics]
rooting = dict(snakemake.params.rooting)

ftdir = str(Path(ftdir).expanduser().resolve())
if not os.path.isfile(os.path.join(ftdir, "workflow", "fold_tree")):
    raise FileNotFoundError(f"FoldTree workflow not found under {ftdir}")

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
famroot = os.path.abspath(os.path.dirname(outputs[0]))
structures = {}
for line in open(famfile):
    p = line.strip()
    if not p: continue
    acc = accre.search(os.path.basename(p))
    src = p if os.path.exists(p) else ((glob.glob(os.path.join(pdb_dir, f"*{acc.group(0)}*.pdb"))[:1] or [None])[0] if acc else None)
    if src and acc:
        structures[acc.group(0)] = src
status = run_foldtree_family(
    family=fam,
    structures=structures,
    family_root=famroot,
    output_paths=outputs,
    metrics=metrics,
    foldtree_dir=ftdir,
    snakemake_bin=smk,
    foldseek_bin=foldseek,
    extra_path=extra_path,
    rooting=rooting,
)
Path(status_output).parent.mkdir(parents=True, exist_ok=True)
temporary_status = status_output + ".tmp"
with open(temporary_status, "w", encoding="utf-8") as status_handle:
    json.dump(status, status_handle, indent=2, sort_keys=True)
    status_handle.write("\n")
os.replace(temporary_status, status_output)

print(
    f"{fam} foldtree: metrics {metrics}; "
    f"rc={status['nested_return_code']}; "
    f"small_family={status['small_family_policy']}"
)
