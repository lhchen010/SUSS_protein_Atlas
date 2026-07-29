"""Reusable FoldTree runner for full-length and cropped-domain families."""

import datetime
import glob
import json
import os
import shutil
import subprocess
from pathlib import Path

from foldtree_utils import TreeValidationError, recover_metric


def _link_or_copy(source, destination):
    try:
        return os.link(source, destination)
    except OSError:
        return shutil.copy2(source, destination)


def run_foldtree_family(
    *,
    family,
    structures,
    family_root,
    output_paths,
    metrics,
    foldtree_dir,
    snakemake_bin,
    foldseek_bin,
    extra_path="",
    rooting=None,
    timeout=1800,
):
    """Run FoldTree for a mapping of safe leaf label to structure path."""
    rooting = dict(rooting or {})
    fallback = str(rooting.get("fallback", "midpoint"))
    small_family_max = int(rooting.get("small_family_max", 3))
    family_root = Path(family_root).resolve()
    family_root.mkdir(parents=True, exist_ok=True)
    structure_dir = family_root / "structs"
    structure_dir.mkdir(exist_ok=True)

    labels = []
    for label, source in structures.items():
        source = Path(source)
        if not source.is_file():
            continue
        destination = structure_dir / f"{label}.pdb"
        if source.resolve() != destination.resolve():
            shutil.copy2(source, destination)
        labels.append(str(label))
    labels = sorted(set(labels))
    (family_root / "identifiers.txt").write_text(
        "\n".join(labels) + "\n", encoding="utf-8"
    )
    if len(labels) < 3:
        raise ValueError(
            f"{family}: FoldTree requires at least 3 structures; found {len(labels)}"
        )

    foldtree_dir = Path(foldtree_dir).expanduser().resolve()
    if not (foldtree_dir / "workflow" / "fold_tree").is_file():
        raise FileNotFoundError(f"FoldTree workflow not found under {foldtree_dir}")

    work_root = family_root.parent.parent
    package = work_root / ".foldtree_pkg" / family
    if not (package / "workflow").is_dir():
        package.parent.mkdir(parents=True, exist_ok=True)
        shutil.copytree(
            foldtree_dir,
            package,
            symlinks=True,
            copy_function=_link_or_copy,
        )

    small_family = len(labels) <= small_family_max
    suffix = (
        "_struct_tree.PP.nwk"
        if small_family
        else "_struct_tree.PP.nwk.rooted.final"
    )
    targets = [str(family_root / f"{metric}{suffix}") for metric in metrics]
    command = [
        str(snakemake_bin),
        "-s",
        "workflow/fold_tree",
        "--cores",
        "4",
        "--keep-going",
        *targets,
        "--config",
        f"folder={family_root}",
        "filter=False",
        "custom_structs=True",
        f"foldseek_path={foldseek_bin}",
    ]
    environment = dict(os.environ)
    if extra_path:
        environment["PATH"] = str(extra_path) + os.pathsep + environment.get("PATH", "")

    result = subprocess.run(
        command,
        cwd=package,
        env=environment,
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
    )
    attempted_at = datetime.datetime.now(datetime.timezone.utc).isoformat()
    log_path = family_root / "foldtree_subworkflow.log"
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            f"\n{'=' * 72}\nATTEMPT {attempted_at}\n"
            f"$ {' '.join(command)}\nRETURN CODE {result.returncode}\n\n"
            f"STDOUT\n{result.stdout}\nSTDERR\n{result.stderr}\n"
        )

    metric_status = {}
    failures = []
    for metric, output_path in zip(metrics, output_paths):
        try:
            metric_status[metric] = recover_metric(
                str(family_root),
                str(metric),
                str(output_path),
                labels,
                small_family=small_family,
                fallback=fallback,
            )
        except TreeValidationError as exc:
            failures.append(str(exc))
            metric_status[metric] = {
                "status": "failed",
                "rooting_method": None,
                "reason": str(exc),
            }

    status = {
        "family": family,
        "n_members": len(labels),
        "small_family_policy": small_family,
        "nested_return_code": result.returncode,
        "attempted_at": attempted_at,
        "log": str(log_path),
        "metrics": metric_status,
        "produced_files": sorted(
            glob.glob(str(family_root / "**" / "*.nwk*"), recursive=True)
        ),
    }
    with log_path.open("a", encoding="utf-8") as handle:
        handle.write(
            "\nRECOVERY STATUS\n"
            + json.dumps(status, indent=2, sort_keys=True)
            + "\n"
        )
    if failures:
        stderr_tail = "\n".join(
            (result.stderr or result.stdout).splitlines()[-30:]
        )
        raise RuntimeError(
            f"{family}: unrecoverable FoldTree metric(s): {'; '.join(failures)}\n"
            f"Nested return code: {result.returncode}; log: {log_path}\n"
            f"{stderr_tail}"
        )
    return status
