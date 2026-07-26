"""Run FoldTree with validated, per-metric rooting recovery."""
import datetime
import json
import os, glob, re, shutil, subprocess
from pathlib import Path
from foldtree_utils import TreeValidationError, recover_metric
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
fallback = str(rooting.get("fallback", "midpoint"))
small_family_max = int(rooting.get("small_family_max", 3))

ftdir = str(Path(ftdir).expanduser().resolve())
if not os.path.isfile(os.path.join(ftdir, "workflow", "fold_tree")):
    raise FileNotFoundError(f"FoldTree workflow not found under {ftdir}")

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
famroot = os.path.abspath(os.path.dirname(outputs[0]))   # the tool's `folder`
struct  = os.path.join(famroot, "structs"); os.makedirs(struct, exist_ok=True)
accs = []
for line in open(famfile):
    p = line.strip()
    if not p: continue
    acc = accre.search(os.path.basename(p))
    src = p if os.path.exists(p) else ((glob.glob(os.path.join(pdb_dir, f"*{acc.group(0)}*.pdb"))[:1] or [None])[0] if acc else None)
    if src and acc:
        shutil.copy(src, os.path.join(struct, acc.group(0) + ".pdb"))
        accs.append(acc.group(0))
# custom_structs mode still needs {folder}/identifiers.txt to exist to seed the DAG
# (dl_ids_sequences reads it, then just writes an empty seq file in custom mode)
with open(os.path.join(famroot, "identifiers.txt"), "w") as fh:
    fh.write("\n".join(accs) + "\n")
if len(accs) < 3:
    raise ValueError(f"{fam}: FoldTree requires at least 3 structures; found {len(accs)}")

work_root = os.path.abspath(os.path.join(famroot, "..", "..", ".."))  # engine workdir
ftpkg = os.path.join(work_root, ".foldtree_pkg", fam)
small_family = len(accs) <= small_family_max
target_suffix = "_struct_tree.PP.nwk" if small_family else "_struct_tree.PP.nwk.rooted.final"
tree_targets = [f"{famroot}/{m}{target_suffix}" for m in metrics]
cmd = [smk, "-s", "workflow/fold_tree", "--cores", "4", "--keep-going", *tree_targets,
       "--config", f"folder={famroot}", "filter=False", "custom_structs=True",
       f"foldseek_path={foldseek}"]
env = dict(os.environ)
if extra_path:
    env["PATH"] = extra_path + os.pathsep + env.get("PATH", "")
log_path = os.path.join(famroot, "foldtree_subworkflow.log")
if not os.path.isdir(os.path.join(ftpkg, "workflow")):
    def link_or_copy(source, destination):
        try:
            return os.link(source, destination)
        except OSError:
            return shutil.copy2(source, destination)

    # Each nested Snakemake needs private writable metadata. Hard-link static package
    # files where possible so family isolation does not multiply FoldTree's disk usage.
    shutil.copytree(ftdir, ftpkg, symlinks=True, copy_function=link_or_copy)
r = subprocess.run(cmd, cwd=ftpkg, env=env, capture_output=True, text=True,
                   timeout=1800, check=False)
attempt = datetime.datetime.now(datetime.timezone.utc).isoformat()
with open(log_path, "a", encoding="utf-8") as log_handle:
    log_handle.write(
        f"\n{'=' * 72}\nATTEMPT {attempt}\n$ {' '.join(cmd)}\n"
        f"RETURN CODE {r.returncode}\n\nSTDOUT\n{r.stdout}\nSTDERR\n{r.stderr}\n"
    )

# Recover each metric independently. A non-zero nested return code is tolerated
# only when every declared tree can still be validated and normalized.
produced = glob.glob(os.path.join(famroot, "**", "*.nwk*"), recursive=True)
metric_status = {}
failures = []
for o, m in zip(outputs, metrics):
    try:
        metric_status[m] = recover_metric(
            famroot,
            m,
            o,
            accs,
            small_family=small_family,
            fallback=fallback,
        )
    except TreeValidationError as exc:
        failures.append(str(exc))
        metric_status[m] = {
            "status": "failed",
            "rooting_method": None,
            "reason": str(exc),
        }

status = {
    "family": fam,
    "n_members": len(accs),
    "small_family_policy": small_family,
    "nested_return_code": r.returncode,
    "attempted_at": attempt,
    "log": log_path,
    "metrics": metric_status,
}
Path(status_output).parent.mkdir(parents=True, exist_ok=True)
temporary_status = status_output + ".tmp"
with open(temporary_status, "w", encoding="utf-8") as status_handle:
    json.dump(status, status_handle, indent=2, sort_keys=True)
    status_handle.write("\n")
os.replace(temporary_status, status_output)

with open(log_path, "a", encoding="utf-8") as log_handle:
    log_handle.write("\nRECOVERY STATUS\n" + json.dumps(status, indent=2, sort_keys=True) + "\n")

if failures:
    stderr_tail = "\n".join((r.stderr or r.stdout).splitlines()[-30:])
    raise RuntimeError(
        f"{fam}: unrecoverable FoldTree metric(s): {'; '.join(failures)}\n"
        f"Nested return code: {r.returncode}; log: {log_path}\n{stderr_tail}"
    )

print(
    f"{fam} foldtree: metrics {metrics}; produced {len(produced)} nwk; "
    f"rc={r.returncode}; small_family={small_family}"
)
