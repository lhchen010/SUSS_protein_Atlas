"""Build full-featured workbench assets for local structural-domain families."""

import json
import math
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd

from foldtree_runner import run_foldtree_family
from runtime_utils import resolve_executable
from v3_utils import fasta_records, foldmason_column_scores, protein_id


def segment_pdb(path: Path, start: int, end: int) -> str:
    lines = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        try:
            residue = int(line[22:26])
        except (ValueError, IndexError):
            continue
        if start <= residue <= end:
            lines.append(line)
    return "\n".join(lines) + ("\nTER\nEND\n" if lines else "")


def safe_name(segment_id: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "__", segment_id)


def json_safe(value):
    """Convert numeric missing values to explicit JSON nulls."""
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def read_fasta_text(text: str, reverse_ids: dict[str, str]) -> dict[str, str]:
    records = {}
    current = None
    chunks = []
    for line in text.splitlines():
        if line.startswith(">"):
            if current is not None:
                records[reverse_ids.get(current, current)] = "".join(chunks)
            current = line[1:].split()[0]
            chunks = []
        elif current is not None:
            chunks.append(line.strip())
    if current is not None:
        records[reverse_ids.get(current, current)] = "".join(chunks)
    return records


def write_fasta(path: Path, records: dict[str, str], id_map: dict[str, str]):
    path.write_text(
        "".join(
            f">{id_map.get(record_id, record_id)}\n{sequence}\n"
            for record_id, sequence in records.items()
        ),
        encoding="utf-8",
    )


def aligned_identity_matrix(
    records: dict[str, str], labels: list[str]
) -> list[list[float | None]]:
    """Pairwise AA identity over columns where both domain segments have residues."""
    matrix = []
    for left in labels:
        row = []
        for right in labels:
            if left == right:
                row.append(1.0)
                continue
            pairs = [
                (a, b)
                for a, b in zip(records.get(left, ""), records.get(right, ""))
                if a not in {"-", "."} and b not in {"-", "."}
            ]
            row.append(
                round(sum(a == b for a, b in pairs) / len(pairs), 6)
                if pairs
                else None
            )
        matrix.append(row)
    return matrix


def parse_rate4site(path: Path) -> dict[int, float]:
    """Return positive-is-conserved scores indexed by ungapped reference position."""
    scores = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^\s*(\d+)\s+\S\s+(-?[\d.]+)", line.strip())
        if match:
            scores[int(match.group(1))] = -float(match.group(2))
    return scores


def map_conservation_to_segments(
    msa: dict[str, str],
    reference: str,
    reference_scores: dict[int, float],
    starts: dict[str, int],
) -> dict[str, dict[int, float]]:
    """Project reference-column Rate4Site scores onto every aligned segment."""
    reference_alignment = msa.get(reference, "")
    reference_position = 0
    column_scores = []
    for symbol in reference_alignment:
        if symbol not in {"-", "."}:
            reference_position += 1
            column_scores.append(reference_scores.get(reference_position))
        else:
            column_scores.append(None)

    mapped = {}
    for segment_id, aligned in msa.items():
        residue = int(starts[segment_id]) - 1
        values = {}
        for column, symbol in enumerate(aligned):
            if symbol in {"-", "."}:
                continue
            residue += 1
            score = column_scores[column] if column < len(column_scores) else None
            if score is not None and math.isfinite(score):
                values[residue] = round(float(score), 6)
        mapped[segment_id] = values
    return mapped


def parse_usalign(stdout: str) -> tuple[dict, dict]:
    """Parse US-align outfmt 2 summary and the structure1 -> structure2 matrix."""
    rows = []
    in_matrix = False
    summary = {}
    for line in stdout.splitlines():
        if line.startswith("#PDBchain1"):
            continue
        if "rotation matrix" in line.lower():
            in_matrix = True
            continue
        fields = line.split()
        if in_matrix and len(fields) >= 5 and fields[0] in {"0", "1", "2"}:
            rows.append([float(value) for value in fields[1:5]])
        elif (
            not line.startswith("#")
            and len(fields) >= 5
            and not summary
            and all(re.fullmatch(r"[-+.\deE]+", value) for value in fields[2:5])
        ):
            summary = {
                "tm_mobile": float(fields[2]),
                "tm_reference": float(fields[3]),
                "rmsd": float(fields[4]),
            }
    if len(rows) != 3:
        raise ValueError("US-align output did not contain a 3x4 transform matrix")
    translation = [row[0] for row in rows]
    rotation = np.asarray([row[1:] for row in rows], dtype=float).T
    return {
        "rotation": rotation.round(10).tolist(),
        "translation": np.asarray(translation).round(10).tolist(),
    }, summary


def run_usalign(usalign: str, mobile: Path, reference: Path):
    result = subprocess.run(
        [
            usalign,
            str(mobile),
            str(reference),
            "-outfmt",
            "2",
            "-m",
            "-",
        ],
        capture_output=True,
        text=True,
        timeout=180,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"US-align failed ({result.returncode}): {result.stderr[-500:]}"
        )
    return parse_usalign(result.stdout)


def sequence_components(
    segment_ids: list[str],
    accession_by_segment: dict[str, str],
    blast: pd.DataFrame,
    evalue_threshold: float,
    coverage_threshold: float,
) -> list[list[str]]:
    parent = {segment_id: segment_id for segment_id in segment_ids}

    def find(item):
        while parent[item] != item:
            parent[item] = parent[parent[item]]
            item = parent[item]
        return item

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[max(left_root, right_root)] = min(left_root, right_root)

    for index, left in enumerate(segment_ids):
        for right in segment_ids[index + 1:]:
            if accession_by_segment[left] == accession_by_segment[right]:
                union(left, right)

    accessions = set(accession_by_segment.values())
    selected = blast[
        blast.q.isin(accessions)
        & blast.t.isin(accessions)
        & (blast.evalue <= evalue_threshold)
        & (blast.min_coverage >= coverage_threshold)
    ]
    related = {
        tuple(sorted((str(row.q), str(row.t))))
        for row in selected.itertuples()
    }
    for index, left in enumerate(segment_ids):
        for right in segment_ids[index + 1:]:
            pair = tuple(
                sorted(
                    (
                        accession_by_segment[left],
                        accession_by_segment[right],
                    )
                )
            )
            if pair in related:
                union(left, right)

    groups = {}
    for segment_id in segment_ids:
        groups.setdefault(find(segment_id), []).append(segment_id)
    return sorted(
        (sorted(group) for group in groups.values()),
        key=lambda group: (-len(group), group[0]),
    )


def configured_tool(value, name, required=True):
    raw = str(value or "").strip()
    return resolve_executable(raw, name, required=required and bool(raw))


members = pd.read_csv(snakemake.input.members)
edges = pd.read_csv(snakemake.input.edges)
qc = pd.read_csv(snakemake.input.qc)
sequences = fasta_records(snakemake.input.seqs)
pdb_dir = Path(snakemake.input.pdb_dir)
output = Path(snakemake.output.json)
assets_root = Path(snakemake.output.assets)
output.parent.mkdir(parents=True, exist_ok=True)
if assets_root.exists():
    shutil.rmtree(assets_root)
assets_root.mkdir(parents=True)
if not bool(snakemake.params.enabled):
    output.write_text(
        json.dumps({"schema_version": 3, "families": {}}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(f"domain workbench: disabled -> {output}")
    sys.exit(0)

mafft = configured_tool(snakemake.params.mafft, "MAFFT")
fasttree = configured_tool(snakemake.params.fasttree, "FastTree")
usalign = configured_tool(snakemake.params.usalign, "US-align")
foldmason = configured_tool(snakemake.params.foldmason, "FoldMason")
conservation_enabled = bool(snakemake.params.conservation_enabled)
rate4site = (
    configured_tool(snakemake.params.rate4site, "Rate4Site")
    if conservation_enabled
    else ""
)
minimum_conservation_sequences = int(
    snakemake.params.min_rate4site_sequences
)
foldtree_enabled = bool(snakemake.params.foldtree_enabled)
foldtree_dir = str(snakemake.params.foldtree_dir or "").strip()
foldtree_snakemake = (
    configured_tool(snakemake.params.foldtree_snakemake, "FoldTree Snakemake")
    if foldtree_enabled
    else ""
)
foldtree_foldseek = (
    configured_tool(snakemake.params.foldtree_foldseek, "FoldTree Foldseek")
    if foldtree_enabled
    else ""
)
threads = max(1, int(snakemake.threads))
pair_threshold = float(snakemake.params.pair_threshold)

blast_columns = [
    "query",
    "target",
    "pident",
    "alnlen",
    "evalue",
    "bitscore",
    "qlen",
    "slen",
]
try:
    blast = pd.read_csv(snakemake.input.blastp, sep="\t", names=blast_columns)
except (pd.errors.EmptyDataError, FileNotFoundError):
    blast = pd.DataFrame(columns=blast_columns)
if len(blast):
    blast["q"] = blast["query"].map(protein_id)
    blast["t"] = blast["target"].map(protein_id)
    blast["min_coverage"] = np.minimum(
        pd.to_numeric(blast.alnlen, errors="coerce")
        / pd.to_numeric(blast.qlen, errors="coerce"),
        pd.to_numeric(blast.alnlen, errors="coerce")
        / pd.to_numeric(blast.slen, errors="coerce"),
    )
else:
    blast["q"] = pd.Series(dtype=str)
    blast["t"] = pd.Series(dtype=str)
    blast["min_coverage"] = pd.Series(dtype=float)

pdb_paths = {}
for row in qc.itertuples():
    accession = str(row.acc)
    filename = getattr(row, "fn", "")
    candidates = [
        pdb_dir / str(filename),
        pdb_dir / f"{accession}.pdb",
        *sorted(pdb_dir.glob(f"*{accession}*.pdb")),
    ]
    pdb_paths[accession] = next(
        (candidate for candidate in candidates if candidate.is_file()), None
    )

payload = {"schema_version": 3, "families": {}}
for family, group in members.groupby("domain_family", sort=False):
    family = str(family)
    family_dir = assets_root / family
    structure_dir = family_dir / "structures"
    sequence_dir = family_dir / "sequences"
    structure_dir.mkdir(parents=True)
    sequence_dir.mkdir()
    family_edges = edges[edges.domain_family.astype(str) == family].copy()
    segment_ids = sorted(group.segment_id.astype(str))
    scores = {segment_id: 0.0 for segment_id in segment_ids}
    for edge in family_edges.itertuples():
        weight = float(edge.lddt) if pd.notna(edge.lddt) else 0.0
        scores[str(edge.source)] = scores.get(str(edge.source), 0.0) + weight
        scores[str(edge.target)] = scores.get(str(edge.target), 0.0) + weight
    hub = min(segment_ids, key=lambda segment_id: (-scores[segment_id], segment_id))

    segment_files = {}
    segment_sequences = {}
    parent_sequences = {}
    accession_by_segment = {}
    start_by_segment = {}
    id_map = {}
    for row in group.itertuples():
        segment_id = str(row.segment_id)
        accession = str(row.acc)
        safe_id = safe_name(segment_id)
        id_map[segment_id] = safe_id
        accession_by_segment[segment_id] = accession
        start_by_segment[segment_id] = int(row.start)
        path = pdb_paths.get(accession)
        if path:
            text = segment_pdb(path, int(row.start), int(row.end))
            if text:
                segment_path = structure_dir / f"{safe_id}.pdb"
                segment_path.write_text(text, encoding="utf-8")
                segment_files[segment_id] = segment_path
        sequence = sequences.get(accession, "")
        if sequence:
            parent_sequences[accession] = sequence
            segment_sequences[segment_id] = sequence[
                int(row.start) - 1:int(row.end)
            ]
    reverse_ids = {value: key for key, value in id_map.items()}
    write_fasta(
        sequence_dir / f"{family}_domain_segments.fasta",
        segment_sequences,
        id_map,
    )
    (sequence_dir / f"{family}_parent_proteins.fasta").write_text(
        "".join(
            f">{accession}\n{sequence}\n"
            for accession, sequence in sorted(parent_sequences.items())
        ),
        encoding="utf-8",
    )

    status = {
        "foldmason": "not_applicable",
        "sequence": "not_applicable",
        "sequence_conservation": (
            "not_applicable" if conservation_enabled else "disabled"
        ),
        "foldtree": "not_applicable",
        "superposition": "not_applicable",
    }
    structural_msa = {}
    three_di_msa = {}
    foldmason_guide_tree = ""
    structural_scores = []
    if foldmason and len(segment_files) >= 2:
        prefix = family_dir / f"{family}_foldmason"
        temp_dir = family_dir / "foldmason_tmp"
        command = [
            foldmason,
            "easy-msa",
            *[str(segment_files[segment_id]) for segment_id in segment_ids
              if segment_id in segment_files],
            str(prefix),
            str(temp_dir),
            "--report-mode",
            "1",
        ]
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=1800,
            check=False,
        )
        if result.returncode != 0:
            raise RuntimeError(
                f"{family}: FoldMason failed ({result.returncode}): "
                + "\n".join(result.stderr.splitlines()[-20:])
            )
        aa_path = Path(f"{prefix}_aa.fa")
        di_path = Path(f"{prefix}_3di.fa")
        structural_msa = read_fasta_text(
            aa_path.read_text(encoding="utf-8"), reverse_ids
        )
        three_di_msa = read_fasta_text(
            di_path.read_text(encoding="utf-8"), reverse_ids
        )
        guide_path = Path(f"{prefix}.nw")
        if guide_path.is_file():
            foldmason_guide_tree = guide_path.read_text(
                encoding="utf-8"
            ).strip()
        status["foldmason"] = "complete"

        normalized_msa = family_dir / f"{family}_foldmason_normalized.fasta"
        write_fasta(normalized_msa, structural_msa, id_map)
        database = family_dir / f"{family}_foldmason_db"
        score_json = family_dir / f"{family}_column_lddt.json"
        subprocess.run(
            [
                foldmason,
                "createdb",
                str(structure_dir),
                str(database),
                "--threads",
                "1",
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        subprocess.run(
            [
                foldmason,
                "msa2lddtjson",
                str(database),
                str(normalized_msa),
                str(score_json),
                "--pair-threshold",
                str(pair_threshold),
                "--threads",
                str(threads),
            ],
            check=True,
            capture_output=True,
            text=True,
        )
        alignment_length = len(next(iter(structural_msa.values())))
        structural_scores = foldmason_column_scores(
            json.loads(score_json.read_text(encoding="utf-8")),
            alignment_length,
        )

    sequence_subgroups = []
    sequence_conservation = {}
    components = sequence_components(
        segment_ids,
        accession_by_segment,
        blast,
        float(snakemake.params.blast_evalue),
        float(snakemake.params.blast_coverage),
    )
    for subgroup_index, component in enumerate(components):
        subgroup_id = f"{family}.S{subgroup_index}"
        subgroup = {
            "id": subgroup_id,
            "members": component,
            "status": "sequence_singleton",
            "msa": {},
            "newick": "",
            "sequence_conservation_status": "not_applicable",
            "sequence_conservation_reference": "",
        }
        subgroup_sequences = {
            segment_id: segment_sequences[segment_id]
            for segment_id in component
            if segment_id in segment_sequences
        }
        if mafft and len(subgroup_sequences) >= 2:
            input_fasta = sequence_dir / f"{subgroup_id}.fasta"
            write_fasta(input_fasta, subgroup_sequences, id_map)
            result = subprocess.run(
                [mafft, "--auto", "--thread", str(threads), str(input_fasta)],
                capture_output=True,
                text=True,
                timeout=1800,
                check=False,
            )
            if result.returncode != 0 or not result.stdout.lstrip().startswith(">"):
                raise RuntimeError(
                    f"{subgroup_id}: MAFFT failed ({result.returncode}): "
                    + "\n".join(result.stderr.splitlines()[-20:])
                )
            subgroup["msa"] = read_fasta_text(result.stdout, reverse_ids)
            msa_path = sequence_dir / f"{subgroup_id}_MAFFT.fasta"
            msa_path.write_text(result.stdout, encoding="utf-8")
            subgroup["status"] = "complete"
            if fasttree:
                tree_result = subprocess.run(
                    [fasttree, "-wag", str(msa_path)],
                    capture_output=True,
                    text=True,
                    timeout=900,
                    check=False,
                )
                if tree_result.returncode != 0 or ";" not in tree_result.stdout:
                    raise RuntimeError(
                        f"{subgroup_id}: FastTree failed "
                        f"({tree_result.returncode})"
                    )
                subgroup["newick"] = tree_result.stdout.strip()
                (sequence_dir / f"{subgroup_id}_FastTree.nwk").write_text(
                    subgroup["newick"] + "\n", encoding="utf-8"
                )
            if conservation_enabled:
                if len(subgroup_sequences) < minimum_conservation_sequences:
                    subgroup["sequence_conservation_status"] = (
                        f"requires_at_least_{minimum_conservation_sequences}_sequences"
                    )
                else:
                    reference = hub if hub in component else component[0]
                    rate4site_path = (
                        sequence_dir / f"{subgroup_id}_Rate4Site.res"
                    )
                    rate4site_result = subprocess.run(
                        [
                            rate4site,
                            "-s",
                            str(msa_path),
                            "-a",
                            id_map[reference],
                            "-o",
                            str(rate4site_path),
                        ],
                        capture_output=True,
                        text=True,
                        timeout=1800,
                        check=False,
                    )
                    if (
                        rate4site_result.returncode != 0
                        or not rate4site_path.is_file()
                    ):
                        raise RuntimeError(
                            f"{subgroup_id}: Rate4Site failed "
                            f"({rate4site_result.returncode}): "
                            + "\n".join(
                                (
                                    rate4site_result.stderr
                                    + rate4site_result.stdout
                                ).splitlines()[-20:]
                            )
                        )
                    reference_scores = parse_rate4site(rate4site_path)
                    if not reference_scores:
                        raise RuntimeError(
                            f"{subgroup_id}: Rate4Site produced no scores"
                        )
                    mapped = map_conservation_to_segments(
                        subgroup["msa"],
                        reference,
                        reference_scores,
                        start_by_segment,
                    )
                    sequence_conservation.update(mapped)
                    subgroup["sequence_conservation_status"] = "complete"
                    subgroup["sequence_conservation_reference"] = reference
        sequence_subgroups.append(subgroup)
    if any(group["status"] == "complete" for group in sequence_subgroups):
        status["sequence"] = "complete"
    if any(
        group["sequence_conservation_status"] == "complete"
        for group in sequence_subgroups
    ):
        status["sequence_conservation"] = "complete"

    transforms = {}
    fit_stats = {}
    labels = [segment_id for segment_id in segment_ids if segment_id in segment_files]
    usalign_matrix = np.full((len(labels), len(labels)), np.nan)
    np.fill_diagonal(usalign_matrix, 1.0)
    identity = {
        "rotation": np.eye(3).tolist(),
        "translation": [0.0, 0.0, 0.0],
    }
    if hub in segment_files:
        transforms[hub] = identity
        fit_stats[hub] = {
            "reference": hub,
            "method": "reference",
            "tm_mobile": 1.0,
            "tm_reference": 1.0,
            "rmsd": 0.0,
        }
    if usalign and len(labels) >= 2:
        status["superposition"] = "complete"
        pair_jobs = [
            (left_index, right_index, labels[left_index], labels[right_index])
            for left_index in range(len(labels))
            for right_index in range(left_index + 1, len(labels))
        ]
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = {
                executor.submit(
                    run_usalign,
                    usalign,
                    segment_files[left],
                    segment_files[right],
                ): (left_index, right_index, left, right)
                for left_index, right_index, left, right in pair_jobs
            }
            for future in as_completed(futures):
                left_index, right_index, left, right = futures[future]
                transform, stats = future.result()
                score = min(stats["tm_mobile"], stats["tm_reference"])
                usalign_matrix[left_index, right_index] = score
                usalign_matrix[right_index, left_index] = score
                if right == hub:
                    transforms[left] = transform
                    fit_stats[left] = {
                        "reference": hub,
                        "method": "US-align domain segment",
                        **stats,
                    }
                elif left == hub:
                    reverse_transform, reverse_stats = run_usalign(
                        usalign, segment_files[right], segment_files[left]
                    )
                    transforms[right] = reverse_transform
                    fit_stats[right] = {
                        "reference": hub,
                        "method": "US-align domain segment",
                        **reverse_stats,
                    }
        pd.DataFrame(
            usalign_matrix, index=labels, columns=labels
        ).to_csv(family_dir / f"{family}_USalign_TM.csv")

    foldtree_trees = {}
    foldtree_status = {}
    if foldtree_enabled and len(segment_files) >= 3:
        metrics = [str(metric) for metric in snakemake.params.foldtree_metrics]
        foldtree_root = family_dir / "foldtree"
        output_paths = [
            foldtree_root / f"{family}_{metric}.nwk" for metric in metrics
        ]
        safe_structures = {
            id_map[segment_id]: path
            for segment_id, path in segment_files.items()
        }
        foldtree_status = run_foldtree_family(
            family=family,
            structures=safe_structures,
            family_root=foldtree_root,
            output_paths=output_paths,
            metrics=metrics,
            foldtree_dir=foldtree_dir,
            snakemake_bin=foldtree_snakemake,
            foldseek_bin=foldtree_foldseek,
            extra_path=snakemake.params.foldtree_extra_path,
            rooting=dict(snakemake.params.foldtree_rooting),
        )
        foldtree_trees = {
            metric: path.read_text(encoding="utf-8").strip()
            for metric, path in zip(metrics, output_paths)
            if path.is_file()
        }
        status["foldtree"] = "complete"

    lddt_values = pd.to_numeric(
        family_edges.get("lddt", pd.Series(dtype=float)), errors="coerce"
    ).dropna()
    alntm_values = pd.to_numeric(
        family_edges.get("alntmscore", pd.Series(dtype=float)),
        errors="coerce",
    ).dropna()
    primary_subgroup = next(
        (
            group
            for group in sequence_subgroups
            if hub in group["members"] and group["status"] == "complete"
        ),
        next(
            (
                group
                for group in sequence_subgroups
                if group["status"] == "complete"
            ),
            {},
        ),
    )
    hub_structural_conservation = {}
    if (
        hub in structural_msa
        and structural_scores
        and len(structural_msa[hub]) == len(structural_scores)
    ):
        residue = int(
            group.loc[group.segment_id.astype(str) == hub, "start"].iloc[0]
        )
        for symbol, score in zip(structural_msa[hub], structural_scores):
            if symbol in {"-", "."}:
                continue
            if score is not None and math.isfinite(float(score)):
                hub_structural_conservation[residue] = round(float(score), 6)
            residue += 1
    sequence_identity_labels = [
        segment_id for segment_id in segment_ids if segment_id in structural_msa
    ]
    sequence_identity_matrix = aligned_identity_matrix(
        structural_msa, sequence_identity_labels
    )
    if sequence_identity_labels:
        pd.DataFrame(
            sequence_identity_matrix,
            index=sequence_identity_labels,
            columns=sequence_identity_labels,
        ).to_csv(family_dir / f"{family}_domain_sequence_identity.csv")
    conservation_rows = [
        {
            "subgroup": next(
                (
                    group["id"]
                    for group in sequence_subgroups
                    if segment_id in group["members"]
                ),
                "",
            ),
            "segment_id": segment_id,
            "resi": residue,
            "sequence_conservation": score,
        }
        for segment_id, values in sequence_conservation.items()
        for residue, score in values.items()
    ]
    pd.DataFrame(
        conservation_rows,
        columns=[
            "subgroup",
            "segment_id",
            "resi",
            "sequence_conservation",
        ],
    ).to_csv(
        family_dir / f"{family}_domain_sequence_conservation.csv",
        index=False,
    )
    payload["families"][family] = {
        "hub": hub,
        "members": segment_ids,
        "member_ids": id_map,
        "transforms": transforms,
        "fit_stats": fit_stats,
        "structural_msa": structural_msa,
        "three_di_msa": three_di_msa,
        "foldmason_guide_newick": foldmason_guide_tree,
        "structural_conservation": structural_scores,
        "hub_structural_conservation": hub_structural_conservation,
        "sequence_conservation": sequence_conservation,
        "sequence_identity_labels": sequence_identity_labels,
        "sequence_identity_matrix": sequence_identity_matrix,
        "sequence_identity_method": "FoldMason-aligned domain amino-acid identity",
        "sequence_msa": primary_subgroup.get("msa", {}),
        "sequence_newick": primary_subgroup.get("newick", ""),
        "sequence_subgroups": sequence_subgroups,
        "usalign_labels": labels,
        "usalign_matrix": [
            [None if not math.isfinite(value) else round(float(value), 6)
             for value in row]
            for row in usalign_matrix
        ],
        "foldtree_trees": foldtree_trees,
        "foldtree_status": foldtree_status,
        "tree_label_map": {value: key for key, value in id_map.items()},
        "status": status,
        "mean_lddt": (
            round(float(lddt_values.mean()), 4) if len(lddt_values) else None
        ),
        "mean_alntm": (
            round(float(alntm_values.mean()), 4) if len(alntm_values) else None
        ),
    }

output.write_text(
    json.dumps(json_safe(payload), indent=2, sort_keys=True, allow_nan=False)
    + "\n",
    encoding="utf-8",
)
print(
    f"domain workbench: {len(payload['families'])} families -> {output}"
)
