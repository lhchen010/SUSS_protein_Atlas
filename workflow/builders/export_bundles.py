"""Excel workbook and downloadable package builders for F and D families."""

import base64
import io
import json
import math
import os
import re
import sys
import zipfile

import numpy as np
import openpyxl
from openpyxl.styles import Font, PatternFill
import pandas as pd

_BUILDERS_DIR = os.path.dirname(os.path.abspath(__file__))
if _BUILDERS_DIR not in sys.path:
    sys.path.insert(0, _BUILDERS_DIR)
from pdb_io import _pdb_with_bfactors


def _clean_cell(value):
    text = "" if value is None else str(value)
    return "" if text.lower() in {"nan", "<na>"} else text


def _finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _domain_xlsx_b64(fam, members, edges, workbench):
    """Build the complete D-family workbook using segment-aware evidence."""
    buf = io.BytesIO()
    member_table = pd.DataFrame(members).copy()
    for column in (
        "overlap_annotations",
        "expression",
        "p2rank",
        "fpocket",
        "esm_values",
    ):
        if column in member_table:
            member_table[column] = member_table[column].map(
                lambda value: json.dumps(value, separators=(",", ":"))
            )
    with pd.ExcelWriter(buf, engine="openpyxl") as xl:
        member_table.to_excel(xl, sheet_name="members", index=False)
        pd.DataFrame(edges).to_excel(
            xl, sheet_name="foldseek_local_links", index=False
        )
        labels = workbench.get("usalign_labels", [])
        matrix = workbench.get("usalign_matrix", [])
        if labels and matrix:
            pd.DataFrame(matrix, index=labels, columns=labels).to_excel(
                xl, sheet_name="usalign_TM"
            )
        sequence_identity_labels = workbench.get(
            "sequence_identity_labels", []
        )
        sequence_identity_matrix = workbench.get(
            "sequence_identity_matrix", []
        )
        if sequence_identity_labels and sequence_identity_matrix:
            pd.DataFrame(
                sequence_identity_matrix,
                index=sequence_identity_labels,
                columns=sequence_identity_labels,
            ).to_excel(xl, sheet_name="sequence_identity")
        structural_identity_labels = workbench.get(
            "structural_alignment_identity_labels", []
        )
        structural_identity_matrix = workbench.get(
            "structural_alignment_identity_matrix", []
        )
        if structural_identity_labels and structural_identity_matrix:
            pd.DataFrame(
                structural_identity_matrix,
                index=structural_identity_labels,
                columns=structural_identity_labels,
            ).to_excel(xl, sheet_name="foldmason_AA_identity")
        pd.DataFrame(workbench.get("domain_blast_edges", [])).to_excel(
            xl, sheet_name="domain_blastp_hits", index=False
        )
        for sheet, records in (
            ("sequence_MSA", workbench.get("sequence_msa", {})),
            ("foldmason_AA", workbench.get("structural_msa", {})),
            ("foldmason_3Di", workbench.get("three_di_msa", {})),
        ):
            if records:
                pd.DataFrame(
                    [
                        {
                            "segment_id": member,
                            "aligned_sequence": sequence,
                            "status": "complete",
                            "reason": "",
                        }
                        for member, sequence in records.items()
                    ]
                ).to_excel(xl, sheet_name=sheet, index=False)
            else:
                reason = (
                    "No reciprocal-coverage sequence subgroup qualified for MAFFT."
                    if sheet == "sequence_MSA"
                    else "Alignment was not produced."
                )
                pd.DataFrame(
                    [
                        {
                            "segment_id": "",
                            "aligned_sequence": "",
                            "status": "not_applicable",
                            "reason": reason,
                        }
                    ]
                ).to_excel(xl, sheet_name=sheet, index=False)
        subgroup_rows = []
        for subgroup in workbench.get("sequence_subgroups", []):
            for segment_id in subgroup.get("members", []):
                subgroup_rows.append(
                    {
                        "subgroup": subgroup.get("id"),
                        "segment_id": segment_id,
                        "status": subgroup.get("status"),
                        "newick": subgroup.get("newick", ""),
                        "aligned_sequence": subgroup.get("msa", {}).get(
                            segment_id, ""
                        ),
                    }
                )
        pd.DataFrame(subgroup_rows).to_excel(
            xl, sheet_name="sequence_subgroups", index=False
        )
        tree_rows = []
        for metric, newick in workbench.get("foldtree_trees", {}).items():
            tree_status = (
                workbench.get("foldtree_status", {})
                .get("metrics", {})
                .get(metric, {})
            )
            tree_rows.append(
                {
                    "tree_type": "FoldTree structural",
                    "metric": metric,
                    "status": tree_status.get("status"),
                    "rooting_method": tree_status.get("rooting_method"),
                    "newick": newick,
                }
            )
        if workbench.get("foldmason_guide_newick"):
            tree_rows.append(
                {
                    "tree_type": "FoldMason guide",
                    "metric": "guide",
                    "status": "complete",
                    "rooting_method": "",
                    "newick": workbench["foldmason_guide_newick"],
                }
            )
        for subgroup in workbench.get("sequence_subgroups", []):
            if subgroup.get("newick"):
                tree_rows.append(
                    {
                        "tree_type": "MAFFT + FastTree sequence",
                        "metric": subgroup.get("id"),
                        "status": subgroup.get("status"),
                        "rooting_method": "",
                        "newick": subgroup["newick"],
                    }
                )
        pd.DataFrame(tree_rows).to_excel(xl, sheet_name="trees", index=False)
        fit_rows = []
        for segment_id, stats in workbench.get("fit_stats", {}).items():
            row = {"segment_id": segment_id, **stats}
            for key in ("rotation", "translation"):
                if key in row:
                    row[key] = json.dumps(row[key], separators=(",", ":"))
            fit_rows.append(row)
        pd.DataFrame(fit_rows).to_excel(
            xl, sheet_name="superposition", index=False
        )
        expression_rows = []
        pocket_rows = []
        annotation_rows = []
        seen_expression = set()
        for member in members:
            accession = member["acc"]
            if accession not in seen_expression:
                seen_expression.add(accession)
                expression_rows.append(
                    {"acc": accession, **member.get("expression", {})}
                )
            for method in ("p2rank", "fpocket"):
                result = member.get(method, {}) or {}
                pocket_rows.append(
                    {
                        "segment_id": member["segment_id"],
                        "acc": accession,
                        "method": method,
                        "profile": (
                            member.get("p2rank_profile")
                            if method == "p2rank"
                            else ""
                        ),
                        "top_score": result.get("top_score"),
                        "top_probability": result.get("top_probability"),
                        "n_pockets": result.get("n_pockets"),
                        "domain_lining_residues": " ".join(
                            map(str, result.get("domain_lining_residues", []))
                        ),
                    }
                )
            annotation_rows.append(
                {
                    "segment_id": member["segment_id"],
                    "acc": accession,
                    "overlap_annotations": "; ".join(
                        f"{item.get('source')}: {item.get('label')} "
                        f"({item.get('start')}-{item.get('end')})"
                        for item in member.get("overlap_annotations", [])
                    ),
                    "gene": member.get("gene", ""),
                    "effectorp": member.get("eff", ""),
                    "tmr": member.get("tmr"),
                    "pfam": member.get("pfam", ""),
                    "interpro": member.get("ipr", ""),
                    "pdb_hit": member.get("pdb", ""),
                    "afdb_hit": member.get("afdb", ""),
                }
            )
        pd.DataFrame(expression_rows).to_excel(
            xl, sheet_name="RNAseq_parent_proteins", index=False
        )
        pd.DataFrame(pocket_rows).to_excel(
            xl, sheet_name="pockets_parent_mapped", index=False
        )
        pd.DataFrame(annotation_rows).to_excel(
            xl, sheet_name="annotation", index=False
        )
        scores = workbench.get("structural_conservation", [])
        if scores:
            pd.DataFrame(
                {
                    "alignment_column": range(1, len(scores) + 1),
                    "foldmason_lddt": scores,
                }
            ).to_excel(
                xl, sheet_name="structural_conservation", index=False
            )
        sequence_conservation_rows = [
            {
                "segment_id": segment_id,
                "resi": residue,
                "sequence_conservation": score,
            }
            for segment_id, values in workbench.get(
                "sequence_conservation", {}
            ).items()
            for residue, score in values.items()
        ]
        pd.DataFrame(
            sequence_conservation_rows,
            columns=["segment_id", "resi", "sequence_conservation"],
        ).to_excel(
            xl, sheet_name="sequence_conservation", index=False
        )
        pd.DataFrame(
            [
                ("members", "D-family segment coordinates, parent F family, and evidence."),
                ("foldseek_local_links", "Retained local Foldseek 3Di+AA links defining the D family."),
                ("usalign_TM", "Independent all-pairs US-align TM scores on cropped segments."),
                ("sequence_MSA", "MAFFT MSA for the eligible sequence-homologous subgroup."),
                ("sequence_identity", "Best-HSP BLASTp identity from the independently searched domain-segment FASTA."),
                ("domain_blastp_hits", "Domain-segment BLASTp hits used to define sequence-homologous D subgroups."),
                ("foldmason_AA_identity", "Amino-acid identity over the structural FoldMason alignment; kept separate from sequence search."),
                ("foldmason_AA", "FoldMason structure-guided amino-acid MSA for all segments."),
                ("foldmason_3Di", "FoldMason 3Di structural-alphabet MSA for all segments."),
                ("sequence_conservation", "Rate4Site conservation projected only within eligible sequence-homologous subgroups."),
                ("trees", "FoldTree structural trees, FoldMason guide tree, and sequence trees are separate records."),
                ("RNAseq_parent_proteins", "RNA-seq is protein-level evidence; each parent accession appears once."),
                ("pockets_parent_mapped", "Pockets were predicted on complete parents and mapped to segment coordinates."),
            ],
            columns=["sheet", "contents"],
        ).to_excel(xl, sheet_name="README", index=False)
    return base64.b64encode(buf.getvalue()).decode()


def _domain_package_b64(
    fam,
    members,
    edges,
    workbench,
    parent_structures,
    parent_sequences,
    workbook_b64,
):
    """Build one auditable ZIP containing every D-family sequence and structure asset."""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        root = f"{fam}_domain_family"
        domain_fasta = []
        parent_fasta = []
        written_parents = set()
        written_parent_structures = set()
        for member in members:
            accession = member["acc"]
            sequence = parent_sequences.get(accession, "")
            if sequence:
                domain_fasta.append(
                    f">{member['segment_id']}\n"
                    f"{sequence[int(member['start']) - 1:int(member['end'])]}\n"
                )
                if accession not in written_parents:
                    written_parents.add(accession)
                    parent_fasta.append(f">{accession}\n{sequence}\n")
            parent = parent_structures.get(accession, "")
            if parent:
                safe = re.sub(
                    r"[^A-Za-z0-9_.-]+", "_", member["segment_id"]
                )
                cropped = "\n".join(
                    line for line in parent.splitlines()
                    if line.startswith(("ATOM", "HETATM"))
                    and line[22:26].strip().lstrip("-").isdigit()
                    and int(member["start"]) <= int(line[22:26]) <= int(member["end"])
                )
                archive.writestr(
                    f"{root}/structures/domains/{safe}.pdb",
                    cropped.rstrip() + "\nTER\nEND\n",
                )
                if accession not in written_parent_structures:
                    written_parent_structures.add(accession)
                    archive.writestr(
                        f"{root}/structures/parents/{accession}.pdb",
                        parent.rstrip() + "\n",
                    )
        archive.writestr(
            f"{root}/sequences/domain_segments.fasta", "".join(domain_fasta)
        )
        archive.writestr(
            f"{root}/sequences/parent_proteins.fasta", "".join(parent_fasta)
        )
        for filename, records in (
            ("MAFFT_sequence_MSA.fasta", workbench.get("sequence_msa", {})),
            ("FoldMason_AA_MSA.fasta", workbench.get("structural_msa", {})),
            ("FoldMason_3Di_MSA.fasta", workbench.get("three_di_msa", {})),
        ):
            if records:
                archive.writestr(
                    f"{root}/alignments/{filename}",
                    "".join(
                        f">{member}\n{sequence}\n"
                        for member, sequence in records.items()
                    ),
                )
        for metric, newick in workbench.get("foldtree_trees", {}).items():
            archive.writestr(
                f"{root}/trees/FoldTree_{metric}.nwk",
                newick.rstrip() + "\n",
            )
        for subgroup in workbench.get("sequence_subgroups", []):
            if subgroup.get("msa"):
                archive.writestr(
                    f"{root}/alignments/sequence_subgroups/"
                    f"{subgroup['id']}_MAFFT.fasta",
                    "".join(
                        f">{member}\n{sequence}\n"
                        for member, sequence in subgroup["msa"].items()
                    ),
                )
            if subgroup.get("newick"):
                archive.writestr(
                    f"{root}/trees/{subgroup['id']}_sequence_FastTree.nwk",
                    subgroup["newick"].rstrip() + "\n",
                )
        if workbench.get("foldmason_guide_newick"):
            archive.writestr(
                f"{root}/trees/FoldMason_guide.nwk",
                workbench["foldmason_guide_newick"].rstrip() + "\n",
            )
        archive.writestr(
            f"{root}/tables/members.csv", pd.DataFrame(members).to_csv(index=False)
        )
        archive.writestr(
            f"{root}/tables/foldseek_local_links.csv",
            pd.DataFrame(edges).to_csv(index=False),
        )
        labels = workbench.get("usalign_labels", [])
        matrix = workbench.get("usalign_matrix", [])
        if labels and matrix:
            archive.writestr(
                f"{root}/tables/usalign_TM.csv",
                pd.DataFrame(matrix, index=labels, columns=labels).to_csv(),
            )
        sequence_identity_labels = workbench.get(
            "sequence_identity_labels", []
        )
        sequence_identity_matrix = workbench.get(
            "sequence_identity_matrix", []
        )
        if sequence_identity_labels and sequence_identity_matrix:
            archive.writestr(
                f"{root}/tables/domain_sequence_identity.csv",
                pd.DataFrame(
                    sequence_identity_matrix,
                    index=sequence_identity_labels,
                    columns=sequence_identity_labels,
                ).to_csv(),
            )
        structural_identity_labels = workbench.get(
            "structural_alignment_identity_labels", []
        )
        structural_identity_matrix = workbench.get(
            "structural_alignment_identity_matrix", []
        )
        if structural_identity_labels and structural_identity_matrix:
            archive.writestr(
                f"{root}/tables/foldmason_AA_identity.csv",
                pd.DataFrame(
                    structural_identity_matrix,
                    index=structural_identity_labels,
                    columns=structural_identity_labels,
                ).to_csv(),
            )
        archive.writestr(
            f"{root}/tables/domain_blastp_hits.csv",
            pd.DataFrame(workbench.get("domain_blast_edges", [])).to_csv(
                index=False
            ),
        )
        structural_scores = workbench.get("structural_conservation", [])
        if structural_scores:
            archive.writestr(
                f"{root}/tables/foldmason_structural_conservation.csv",
                pd.DataFrame(
                    {
                        "alignment_column": range(
                            1, len(structural_scores) + 1
                        ),
                        "foldmason_lddt": structural_scores,
                    }
                ).to_csv(index=False),
            )
        sequence_conservation_rows = [
            {
                "segment_id": segment_id,
                "resi": residue,
                "sequence_conservation": score,
            }
            for segment_id, values in workbench.get(
                "sequence_conservation", {}
            ).items()
            for residue, score in values.items()
        ]
        archive.writestr(
            f"{root}/tables/rate4site_sequence_conservation.csv",
            pd.DataFrame(
                sequence_conservation_rows,
                columns=["segment_id", "resi", "sequence_conservation"],
            ).to_csv(index=False),
        )
        archive.writestr(
            f"{root}/superposition/transforms.json",
            json.dumps(
                {
                    "hub": workbench.get("hub"),
                    "transforms": workbench.get("transforms", {}),
                    "fit_stats": workbench.get("fit_stats", {}),
                },
                indent=2,
            )
            + "\n",
        )
        archive.writestr(
            f"{root}/status.json",
            json.dumps(
                {
                    "analyses": workbench.get("status", {}),
                    "foldtree": workbench.get("foldtree_status", {}),
                },
                indent=2,
            )
            + "\n",
        )
        archive.writestr(
            f"{root}/{fam}_data.xlsx", base64.b64decode(workbook_b64)
        )
    return base64.b64encode(buf.getvalue()).decode()


def _xlsx_b64(fam, members, annotation, tm, usm, idm, blast_pairs, sig, exp,
              pocket_entry, pocket_raw, trees, tree_status, fit_stats,
              analysis_kind="family", sequence_msa=None, structural_msa=None,
              three_di_msa=None, sequence_status=None, domains=None,
              subgroups=None, structural_cons=None):
    """Build the complete, auditable per-family analysis workbook."""
    try:
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as xl:
            pd.DataFrame({"family": [fam] * len(members), "member": members}).to_excel(
                xl, sheet_name="members", index=False)
            (annotation if annotation is not None else pd.DataFrame()).to_excel(
                xl, sheet_name="annotation", index=False)
            if tm is not None:
                tm.to_excel(xl, sheet_name="foldseek_TM", index=False)
            if usm is not None:
                usm.to_excel(xl, sheet_name="usalign_TM", index=False)
            if idm is not None:
                idm.to_excel(xl, sheet_name="blast_identity", index=False)
            if blast_pairs is not None and len(blast_pairs):
                blast_pairs.to_excel(xl, sheet_name="blast_pairs", index=False)
            if sig is not None:
                sig.to_excel(xl, sheet_name="per_site", index=False)
            if structural_cons is not None:
                structural_cons.to_excel(
                    xl, sheet_name="structural_conservation", index=False
                )
            for sheet_name, records in (
                ("sequence_MSA", sequence_msa),
                ("structural_MSA_AA", structural_msa),
                ("structural_MSA_3Di", three_di_msa),
            ):
                if records:
                    pd.DataFrame(
                        [{"member": member, "aligned_sequence": sequence}
                         for member, sequence in records.items()]
                    ).to_excel(xl, sheet_name=sheet_name, index=False)
            if sequence_status:
                pd.DataFrame([sequence_status]).to_excel(
                    xl, sheet_name="sequence_analysis", index=False
                )
            if domains is not None and len(domains):
                domains.to_excel(xl, sheet_name="domain_families", index=False)
            if subgroups is not None and len(subgroups):
                subgroups.to_excel(xl, sheet_name="sequence_subgroups", index=False)
            if exp is not None:
                exp.to_excel(xl, sheet_name="RNAseq", index=False)

            pocket_entry = pocket_entry or {}
            summaries = []
            predictions = []
            for method in ("fpocket", "p2rank"):
                result = pocket_entry.get(method, {}) or {}
                summaries.append(dict(
                    family=fam, reference=pocket_entry.get("ref", ""), method=method,
                    status=pocket_entry.get(f"{method}_status", "not_run"),
                    profile=(pocket_entry.get("p2rank_profile") if method == "p2rank" else None),
                    n_pockets=result.get("n_pockets"), top_score=result.get("top_score"),
                    top_probability=result.get("top_probability"),
                    top_lining_residues=" ".join(map(str, result.get("lining_residues", [])))))
                for pred in result.get("pockets", []):
                    predictions.append(dict(
                        family=fam, reference=pocket_entry.get("ref", ""), method=method,
                        pocket_id=pred.get("pocket_id"), score=pred.get("score"),
                        probability=pred.get("probability"),
                        n_residues=len(pred.get("lining_residues", [])),
                        lining_residues=" ".join(map(str, pred.get("lining_residues", [])))))
            pd.DataFrame(summaries).to_excel(xl, sheet_name="pocket_summary", index=False)
            pd.DataFrame(predictions, columns=["family", "reference", "method", "pocket_id", "score",
                                                      "probability", "n_residues", "lining_residues"]).to_excel(
                xl, sheet_name="pocket_predictions", index=False)
            residue_rows = []
            for pred in predictions:
                for residue in str(pred["lining_residues"]).split():
                    residue_rows.append({"family": fam, "reference": pred["reference"],
                                         "method": pred["method"], "pocket_id": pred["pocket_id"],
                                         "residue_number": int(residue)})
            pd.DataFrame(residue_rows, columns=["family", "reference", "method", "pocket_id",
                                                "residue_number"]).to_excel(
                xl, sheet_name="pocket_residues", index=False)
            for sheet, table in (pocket_raw or {}).items():
                table.to_excel(xl, sheet_name=sheet[:31], index=False)

            if analysis_kind == "family":
                metric_status = (tree_status or {}).get("metrics", {})
                tree_rows = []
                for metric, newick in trees.items():
                    status = metric_status.get(metric, {})
                    tree_rows.append({
                        "metric": metric,
                        "status": status.get("status", "status_unavailable"),
                        "rooting_method": status.get("rooting_method"),
                        "source_stage": status.get("source_stage"),
                        "reason": status.get("reason"),
                        "newick": newick,
                    })
                pd.DataFrame(
                    tree_rows,
                    columns=["metric", "status", "rooting_method", "source_stage", "reason", "newick"],
                ).to_excel(xl, sheet_name="foldtree", index=False)
                fit_rows = []
                for member, stats in fit_stats.items():
                    row = {"member": member, **stats}
                    for key in ("rotation", "translation"):
                        if key in row:
                            row[key] = json.dumps(row[key], separators=(",", ":"))
                    fit_rows.append(row)
                pd.DataFrame(fit_rows).to_excel(
                    xl, sheet_name="superposition", index=False)
            readme_rows = [
                ("annotation", "Complete per-member annotation and component statuses from member_annotation.csv."),
                ("pocket_summary", "Detector status and top-pocket summary for fpocket and P2Rank."),
                ("pocket_predictions", "Every pocket reported by each detector, including scores and lining residues."),
                ("pocket_residues", "One row per detector, pocket, and lining residue."),
                ("fpocket_pockets", "All descriptors parsed from the detector-native fpocket info file."),
                ("p2rank_pockets", "Complete detector-native P2Rank predictions table with all original columns."),
                ("RNAseq", "Per-member, replicate-collapsed RNA-seq expression for this family."),
                ("per_site", "Reference-residue conservation, SASA, pocket, and other site-level evidence."),
                ("structural_conservation", "FoldMason-column occupancy, AA/3Di entropy, and official per-column LDDT with pair-support masking."),
                ("sequence_MSA", "MAFFT alignment of the reference protein's sequence-homologous subgroup."),
                ("structural_MSA_AA", "FoldMason structure-guided amino-acid alignment."),
                ("structural_MSA_3Di", "FoldMason alignment represented in the 3Di structural alphabet."),
                ("sequence_analysis", "Applicability and tool status for MAFFT, sequence tree, and Rate4Site."),
                ("domain_families", "Local 3Di+AA domain-family segments overlapping this family."),
                ("sequence_subgroups", "Sequence-homologous subgroups within the structural family."),
            ]
            if analysis_kind == "family":
                readme_rows.extend([
                    ("foldseek_TM", "Within-family symmetric Foldseek TM-score matrix used for structural clustering."),
                    ("usalign_TM", "Independent within-family US-align TM-score matrix."),
                    ("blast_identity", "Best-HSP BLASTp identity matrix."),
                    ("blast_pairs", "Pair-level BLAST and SUSS classification data for structural edges."),
                    ("foldtree", "FoldTree Newick trees with per-metric rooting method and recovery status."),
                    ("superposition", "Hub-referenced rigid-body fit method, paired CA count, and RMSD."),
                ])
            pd.DataFrame(readme_rows, columns=["sheet", "contents"]).to_excel(
                xl, sheet_name="README", index=False)
        return base64.b64encode(buf.getvalue()).decode()
    except Exception as exc:
        raise RuntimeError(f"{fam}: failed to build family workbook") from exc

