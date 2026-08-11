"""Create the canonical QC-passing structure directory used by all Foldseek searches."""

from pathlib import Path

from runtime_utils import stage_qc_structures


staged = stage_qc_structures(
    snakemake.input.qc,
    snakemake.input.pdb_dir,
    snakemake.output.pdb_dir,
)
manifest = Path(snakemake.output.manifest)
manifest.parent.mkdir(parents=True, exist_ok=True)
manifest.write_text("filename\n" + "\n".join(staged) + "\n", encoding="utf-8")
print(f"staged {len(staged)} QC-passing structures for Foldseek")
