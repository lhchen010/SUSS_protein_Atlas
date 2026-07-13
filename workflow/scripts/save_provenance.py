"""Write the effective Snakemake config plus reproducibility metadata."""

import copy
import datetime as dt
import hashlib
import os
import subprocess
from pathlib import Path

import yaml

from runtime_utils import resolve_executable


def sha256_file(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def directory_digest(path, suffix=None):
    digest = hashlib.sha256()
    files = sorted(
        p for p in Path(path).rglob("*")
        if p.is_file()
        and not p.name.startswith("._")
        and "__MACOSX" not in p.parts
        and (suffix is None or p.suffix == suffix)
    )
    for file_path in files:
        digest.update(str(file_path.relative_to(path)).encode())
        digest.update(sha256_file(file_path).encode())
    return digest.hexdigest(), len(files)


def git_commit():
    try:
        return subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                              timeout=5, check=True).stdout.strip()
    except (FileNotFoundError, subprocess.SubprocessError):
        for candidate in (Path("GIT_COMMIT"), Path(__file__).resolve().parents[2] / "GIT_COMMIT"):
            if candidate.is_file():
                return candidate.read_text(encoding="utf-8").strip() or None
        return None


def engine_version():
    for candidate in (Path("VERSION"), Path(__file__).resolve().parents[2] / "VERSION"):
        if candidate.is_file():
            return candidate.read_text(encoding="utf-8").strip()
    return None


def tool_record(name, value):
    raw = str(value or "").strip()
    if not raw:
        return {"configured": raw, "status": "not_configured"}
    if name in {"p2rank_java_env", "foldtree_extra_path"}:
        return {"configured": raw, "status": "configured"}
    path = Path(raw).expanduser()
    if path.is_dir():
        return {"configured": raw, "resolved": str(path.resolve()), "status": "directory"}
    resolved = resolve_executable(raw, name, required=False)
    if resolved:
        return {"configured": raw, "resolved": resolved, "status": "available"}
    if path.is_file():
        return {"configured": raw, "resolved": str(path.resolve()), "status": "file"}
    return {"configured": raw, "status": "not_found"}


effective = copy.deepcopy(dict(snakemake.params.config))
pdb_digest, pdb_count = directory_digest(snakemake.input.pdb_dir, suffix=".pdb")
effective["provenance"] = {
    "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    "engine_version": engine_version(),
    "git_commit": git_commit(),
    "inputs": {
        "pdb_sha256": pdb_digest,
        "pdb_count": pdb_count,
        "seqs_sha256": sha256_file(snakemake.input.seqs),
    },
    "tools": {name: tool_record(name, value) for name, value in dict(snakemake.params.tools).items()},
}
os.makedirs(os.path.dirname(snakemake.output[0]), exist_ok=True)
with open(snakemake.output[0], "w", encoding="utf-8") as handle:
    yaml.safe_dump(effective, handle, sort_keys=False, allow_unicode=True)
