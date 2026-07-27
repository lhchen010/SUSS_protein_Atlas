"""Run sequence MSA, sequence tree, and Rate4Site for the hub's homologous subgroup."""

import json
import os
import re
import subprocess
from pathlib import Path

import pandas as pd

from runtime_utils import resolve_executable
from v3_utils import fasta_records, protein_id, write_fasta


family = snakemake.wildcards.fam
members_file = Path(snakemake.input.famfile)
hub = protein_id(members_file.read_text(encoding="utf-8").splitlines()[0])
subgroups = pd.read_csv(snakemake.input.subgroups)
family_subgroups = subgroups[subgroups.family.astype(str) == str(family)]
hub_row = family_subgroups[family_subgroups.acc.astype(str) == hub]
minimum = int(snakemake.params.min_sequences)
enabled = bool(snakemake.params.enabled)
conservation_enabled = bool(snakemake.params.conservation_enabled)

outputs = {
    "msa": Path(snakemake.output.msa),
    "tree": Path(snakemake.output.tree),
    "r4s": Path(snakemake.output.r4s),
    "ref": Path(snakemake.output.ref),
    "status": Path(snakemake.output.status),
}
for path in outputs.values():
    path.parent.mkdir(parents=True, exist_ok=True)

status = {
    "family": family,
    "reference": hub,
    "sequence_subgroup": None,
    "n_sequences": 0,
    "msa_status": "not_run",
    "tree_status": "not_run",
    "rate4site_status": "not_run",
    "reason": None,
}
for key in ("msa", "tree", "r4s"):
    outputs[key].write_text("", encoding="utf-8")
outputs["ref"].write_text(hub + "\n", encoding="utf-8")

if not enabled:
    status["reason"] = "sequence_conservation_disabled"
elif hub_row.empty:
    status["reason"] = "reference_not_in_sequence_subgroups"
else:
    subgroup = str(hub_row.iloc[0].sequence_subgroup)
    accessions = family_subgroups[
        family_subgroups.sequence_subgroup.astype(str) == subgroup
    ].acc.astype(str).tolist()
    status["sequence_subgroup"] = subgroup
    status["n_sequences"] = len(accessions)
    if len(accessions) < minimum:
        status["reason"] = f"requires_at_least_{minimum}_homologous_sequences"
        status["msa_status"] = "not_applicable"
        status["tree_status"] = "not_applicable"
        status["rate4site_status"] = "not_applicable"
    else:
        sequences = fasta_records(snakemake.input.seqs)
        selected = {
            accession: sequences[accession]
            for accession in accessions
            if accession in sequences
        }
        if len(selected) < minimum:
            status["reason"] = "homologous_sequences_missing_from_fasta"
            status["msa_status"] = "failed"
            raise RuntimeError(
                f"{family}: only {len(selected)}/{len(accessions)} subgroup sequences "
                "were present in the canonical FASTA"
            )

        raw = outputs["msa"].with_suffix(".input.fasta")
        write_fasta(selected, raw)
        mafft = resolve_executable(snakemake.params.mafft, "MAFFT")
        command = [
            mafft,
            "--thread",
            str(max(1, int(snakemake.threads))),
            "--localpair",
            "--maxiterate",
            "1000",
            str(raw),
        ]
        result = subprocess.run(
            command, capture_output=True, text=True, timeout=1800, check=False
        )
        if result.returncode != 0 or not result.stdout.lstrip().startswith(">"):
            raise RuntimeError(
                f"{family}: MAFFT failed ({result.returncode}): "
                + "\n".join(result.stderr.splitlines()[-20:])
            )
        outputs["msa"].write_text(result.stdout, encoding="utf-8")
        raw.unlink(missing_ok=True)
        status["msa_status"] = "complete"

        fasttree = resolve_executable(
            snakemake.params.fasttree, "FastTree", required=False
        )
        if fasttree:
            tree_result = subprocess.run(
                [fasttree, "-wag", str(outputs["msa"])],
                capture_output=True,
                text=True,
                timeout=900,
                check=False,
            )
            if tree_result.returncode != 0 or ";" not in tree_result.stdout:
                raise RuntimeError(
                    f"{family}: FastTree failed ({tree_result.returncode}): "
                    + "\n".join(tree_result.stderr.splitlines()[-20:])
                )
            outputs["tree"].write_text(tree_result.stdout.strip() + "\n", encoding="utf-8")
            status["tree_status"] = "complete"
        else:
            status["tree_status"] = "not_run_tool_missing"

        if conservation_enabled:
            rate4site = resolve_executable(snakemake.params.rate4site, "Rate4Site")
            r4s_result = subprocess.run(
                [
                    rate4site,
                    "-s",
                    str(outputs["msa"]),
                    "-a",
                    hub,
                    "-o",
                    str(outputs["r4s"]),
                ],
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            if r4s_result.returncode != 0 or not outputs["r4s"].exists():
                raise RuntimeError(
                    f"{family}: Rate4Site failed ({r4s_result.returncode}): "
                    + "\n".join((r4s_result.stderr + r4s_result.stdout).splitlines()[-20:])
                )
            parsed = sum(
                bool(re.match(r"^\s*\d+\s+\S\s+-?[\d.]+", line))
                for line in outputs["r4s"].read_text(
                    encoding="utf-8", errors="replace"
                ).splitlines()
            )
            if parsed == 0:
                raise RuntimeError(f"{family}: Rate4Site produced no parseable residue scores")
            status["rate4site_status"] = "complete"
        else:
            status["rate4site_status"] = "disabled"
        status["reason"] = None

outputs["status"].write_text(
    json.dumps(status, indent=2, sort_keys=True) + "\n", encoding="utf-8"
)
print(
    f"{family} sequence analysis: subgroup={status['sequence_subgroup']} "
    f"n={status['n_sequences']} msa={status['msa_status']} "
    f"rate4site={status['rate4site_status']}"
)
