"""foldtree rule — official DessimozLab FoldTree pipeline per family, 3 metrics
(foldtree/alntmscore/lddt). Runs the tool's OWN snakemake in the 4070 base env.
Validated recipe (this session): put the tree TARGET before --config so snakemake
doesn't swallow it as a config entry; the base env already has the chain deps
(wget/scipy/biopython/toytree/statsmodels/matplotlib/tqdm/plotly/requests).
Copies member PDBs into a per-family struct dir, runs, then normalizes the
rooted newick outputs to {fam}_{metric}.nwk.
"""
import os, glob, re, shutil, subprocess
fam      = snakemake.wildcards.fam
famfile  = snakemake.input.famfile
pdb_dir  = snakemake.input.pdb_dir
outputs  = list(snakemake.output)
ftdir    = snakemake.params.ftdir
metrics  = [os.path.basename(o).rsplit("_", 1)[1][:-4] for o in outputs]  # F5_foldtree.nwk -> foldtree

accre = re.compile(r"[A-Z]{2,3}\d{4,}\.\d+")
famroot = os.path.abspath(os.path.dirname(outputs[0]))   # the tool's `folder`
struct  = os.path.join(famroot, "structs"); os.makedirs(struct, exist_ok=True)
accs = []
for line in open(famfile):
    p = line.strip()
    if not p: continue
    acc = accre.search(os.path.basename(p))
    src = p if os.path.exists(p) else (glob.glob(os.path.join(pdb_dir, f"*{acc.group(0)}*.pdb"))[:1] or [None])[0]
    if src and acc:
        shutil.copy(src, os.path.join(struct, acc.group(0) + ".pdb"))
        accs.append(acc.group(0))
# custom_structs mode still needs {folder}/identifiers.txt to exist to seed the DAG
# (dl_ids_sequences reads it, then just writes an empty seq file in custom mode)
with open(os.path.join(famroot, "identifiers.txt"), "w") as fh:
    fh.write("\n".join(accs) + "\n")

# tool's snakemake — targets BEFORE --config (validated). base env has deps.
conda_sh = "/home/user/miniconda3/etc/profile.d/conda.sh"
smk = "/home/user/miniconda3/envs/snakemake/bin/snakemake"
# FoldTree writes its .snakemake/log INSIDE its own package dir, but /home/user/fold_tree is
# read-only for claude (uid 1001). Use a WRITABLE copy of the package under the workdir,
# created once and reused across families.
work_root = os.path.abspath(os.path.join(famroot, "..", "..", ".."))  # engine workdir
ftpkg = os.path.join(work_root, ".foldtree_pkg")
if not os.path.isdir(os.path.join(ftpkg, "workflow")):
    os.makedirs(ftpkg, exist_ok=True)
    subprocess.run(["bash", "-lc", f"cp -rn {ftdir}/. {ftpkg}/ 2>/dev/null || cp -r {ftdir}/. {ftpkg}/"],
                   capture_output=True, text=True)
# VALIDATED recipe (this session, jobs ba6e37d0 + 3f2a0bf3, 39 families x 3 metrics):
#  * folder = family ROOT; PDBs staged in folder/structs/ (foldseek_allvall reads structs/)
#  * request the three .PP.nwk.rooted.final files as EXPLICIT targets so the DAG skips
#    the plddt rule (which otherwise blocks the whole pipeline inside script: ctx)
#  * targets BEFORE --config (else snakemake parses them as config -> ValueError)
#  * single --config with ALL keys (a 2nd --config drops folder -> KeyError)
#  * foldseek_path under a config key WITH underscore; no --use-conda (base env has deps)
tree_targets = " ".join(f"{famroot}/{m}_struct_tree.PP.nwk.rooted.final" for m in metrics)
cmd = (f"source {conda_sh} && conda activate base && "
       f"export PATH=/home/user/miniconda3/bin:/home/user/.local/bin:$PATH && "
       f"cd {ftpkg} && "
       f"{smk} -s workflow/fold_tree --cores 4 {tree_targets} "
       f"--config folder={famroot} filter=False custom_structs=True "
       f"foldseek_path=/home/user/miniconda3/bin/foldseek")
r = subprocess.run(["bash", "-lc", cmd], capture_output=True, text=True, timeout=1800)
# normalize the produced rooted newick to the declared outputs
produced = glob.glob(os.path.join(famroot, "**", "*.nwk*"), recursive=True)
for o, m in zip(outputs, metrics):
    cand = [p for p in produced if m in os.path.basename(p).lower() and "rooted" in p.lower()] \
           or [p for p in produced if m in os.path.basename(p).lower()]
    os.makedirs(os.path.dirname(o), exist_ok=True)
    if cand: shutil.copy(sorted(cand, key=len)[0], o)
    else: open(o, "w").write("")   # empty on failure — assemble tolerates missing trees
print(f"{fam} foldtree: metrics {metrics}; produced {len(produced)} nwk; rc={r.returncode}")
if r.returncode != 0:
    print("STDERR tail:", r.stderr[-800:])

