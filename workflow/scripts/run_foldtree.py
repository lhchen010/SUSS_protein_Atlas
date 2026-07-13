"""Run the official FoldTree sub-workflow and validate every declared tree."""
import os, glob, re, shutil, subprocess
from pathlib import Path
from runtime_utils import resolve_executable
fam      = snakemake.wildcards.fam
famfile  = snakemake.input.famfile
pdb_dir  = snakemake.input.pdb_dir
outputs  = list(snakemake.output)
ftdir    = snakemake.params.ftdir
smk      = resolve_executable(snakemake.params.snakemake, "FoldTree Snakemake")
foldseek = resolve_executable(snakemake.params.foldseek, "FoldTree Foldseek")
extra_path = str(snakemake.params.get("extra_path", "") or "")
metrics  = [os.path.basename(o).rsplit("_", 1)[1][:-4] for o in outputs]  # F5_foldtree.nwk -> foldtree

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
ftpkg = os.path.join(work_root, ".foldtree_pkg")
if not os.path.isdir(os.path.join(ftpkg, "workflow")):
    shutil.copytree(ftdir, ftpkg, dirs_exist_ok=True)

tree_targets = [f"{famroot}/{m}_struct_tree.PP.nwk.rooted.final" for m in metrics]
cmd = [smk, "-s", "workflow/fold_tree", "--cores", "4", *tree_targets,
       "--config", f"folder={famroot}", "filter=False", "custom_structs=True",
       f"foldseek_path={foldseek}"]
env = dict(os.environ)
if extra_path:
    env["PATH"] = extra_path + os.pathsep + env.get("PATH", "")
r = subprocess.run(cmd, cwd=ftpkg, env=env, capture_output=True, text=True,
                   timeout=1800, check=True)
# normalize the produced rooted newick to the declared outputs
produced = glob.glob(os.path.join(famroot, "**", "*.nwk*"), recursive=True)
for o, m in zip(outputs, metrics):
    cand = [p for p in produced if m in os.path.basename(p).lower() and "rooted" in p.lower()] \
           or [p for p in produced if m in os.path.basename(p).lower()]
    os.makedirs(os.path.dirname(o), exist_ok=True)
    if not cand:
        raise RuntimeError(f"{fam}: FoldTree completed but produced no {m} tree")
    source = sorted(cand, key=len)[0]
    if os.path.getsize(source) == 0:
        raise RuntimeError(f"{fam}: FoldTree produced an empty {m} tree")
    os.makedirs(os.path.dirname(o), exist_ok=True)
    shutil.copy(source, o)
print(f"{fam} foldtree: metrics {metrics}; produced {len(produced)} nwk; rc={r.returncode}")
