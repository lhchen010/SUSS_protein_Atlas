"""Run BLASTp all-vs-all on extracted domain segments."""

from pathlib import Path
import shutil
import subprocess
import tempfile

from runtime_utils import resolve_executable


query = Path(snakemake.input.fasta)
output = Path(snakemake.output.tsv)
output.parent.mkdir(parents=True, exist_ok=True)
headers = sum(1 for line in query.open(encoding="utf-8") if line.startswith(">"))
if not bool(snakemake.params.enabled):
    output.write_text("", encoding="utf-8")
    print("domain BLASTp: disabled with sequence MSA/conservation")
elif headers < 2:
    output.write_text("", encoding="utf-8")
    print(f"domain BLASTp: {headers} segment; no pairwise search required")
else:
    blastp = resolve_executable(snakemake.params.blastp, "BLASTp")
    makeblastdb = resolve_executable(
        snakemake.params.makeblastdb, "makeblastdb"
    )
    with tempfile.TemporaryDirectory(
        prefix="domain_blast_", dir=output.parent
    ) as temp_name:
        temp = Path(temp_name)
        database = temp / "domain_segments"
        result = temp / "domain_blastp.tsv"
        subprocess.run(
            [
                makeblastdb,
                "-in",
                str(query),
                "-dbtype",
                "prot",
                "-out",
                str(database),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [
                blastp,
                "-query",
                str(query),
                "-db",
                str(database),
                "-outfmt",
                "6 qseqid sseqid pident length evalue bitscore qlen slen",
                "-evalue",
                str(snakemake.params.evalue),
                "-max_target_seqs",
                "100000",
                "-num_threads",
                str(max(1, int(snakemake.threads))),
                "-out",
                str(result),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        shutil.copyfile(result, output)
    print(f"domain BLASTp: {headers} segments -> {output}")
