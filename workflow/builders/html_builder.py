"""Build the self-contained SUSS Atlas family and singleton workspaces.

The renderer (vis.js family network, singleton evidence table, embedded 3Dmol.js viewer,
matrices, trees, and Blob-based downloads) is stored as four template parts:
    prefix.html    head+body+CSS+vis setup, ends with 'var D='
    databridge.js  maps the generated data payload to renderer variables
    renderer.js    family and singleton interaction functions
    tail.html      '</script>' + closing tags
build_atlas regenerates D = {NET, SINGLETONS, EXTRA, REFPDB, PAY, SUMMARY} and ANN from
the engine's rule outputs, then assembles:  prefix + json(D) + databridge +
'var ANN=' + json(ANN) + renderer + tail.

Fields are populated ONLY when the producing rule ran (step toggled on) and its output
file exists.
  NET.nodes[] = {id, n, tm, id_pct, suss, plddt, len, maxid}         # from master.csv
  NET.edges[] = {from, to, tm, tm_max, n}   # cross-family TM; only if cross_family_edges.csv exists
  SINGLETONS[] = direct per-protein table records, keyed by accession
  PAY[fam]    = {members[], order[], struct{acc:pdbtext}, assets, newick, maxid}
  PAY[acc]    = singleton structure, sequence, RNA-seq image, and direct-evidence workbook
  EXTRA[fam]  = cons_min/cons_max/cons_sasa_r (conservation rule) · pocket_src/pocket_resi/
                pocket_score/n_pocket + p2rank_resi/p2rank_score/p2rank_n/p2rank_prob +
                fpocket_resi/fpocket_score (pocket rule, pockets.json; P2Rank preferred,
                both kept so the viewer can switch) · has_esm/esm_min/esm_max/esm_vs_cons_r/
                esm_vs_sasa_r (esm rule, esm_all.csv) · hub/hub_meanTM (highest mean-TM
                member, gold-starred on the FoldTree) · n_cys (CYS on ref) · ref_used.
                A step toggled off (e.g. pocket:false) → those keys absent; the viewer's
                gating (var pock=EXTRA[curFam]?EXTRA[curFam].pocket_resi:[]) handles it.
  EXTRA[acc]  = singleton pockets, ESM tolerance, pLDDT range, cysteines, and expression
  REFPDB["<fam>_cons"] = conservation-B-factor ref PDB text (conservation rule)
  ANN[fam]    = {label, n, pct_novel, pct_eff, members[]}   # from member_annotation.csv
  ANN[acc]    = direct singleton annotation and component statuses
"""
import os, io, json, glob, base64, math, re, zipfile
from collections import Counter
import numpy as np, pandas as pd

_TPL = os.path.join(os.path.dirname(__file__), "template")


_AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
           'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
           'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def _seq_from_pdb(pdbtext):
    """One-letter sequence from CA atoms of a PDB string (fallback when no fasta)."""
    seq = {}
    for ln in pdbtext.split("\n"):
        if ln.startswith("ATOM") and ln[12:16].strip() == "CA":
            try: resi = int(ln[22:26])
            except ValueError: continue
            seq[resi] = _AA3TO1.get(ln[17:20].strip(), "X")
    return "".join(seq[k] for k in sorted(seq))


def _bfactor_range(pdbtext):
    """Return the finite B-factor range, used as the singleton pLDDT scale."""
    values = []
    for line in pdbtext.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                value = float(line[60:66])
            except ValueError:
                continue
            if math.isfinite(value):
                values.append(value)
    return (min(values), max(values)) if values else (0.0, 100.0)


def _pdb_with_bfactors(pdbtext, values, default=0.0):
    """Return PDB text with residue-indexed values written to the B-factor column."""
    output = []
    for line in pdbtext.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
            try:
                residue = int(line[22:26])
            except ValueError:
                output.append(line)
                continue
            output.append(f"{line[:60]}{float(values.get(residue, default)):6.2f}{line[66:]}")
        else:
            output.append(line)
    return "\n".join(output)


def _clean_cell(value):
    text = "" if value is None else str(value)
    return "" if text.lower() in {"nan", "<na>"} else text


def _finite_float(value):
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _annotation_payload(group):
    """Build the renderer annotation object for a family or one singleton."""
    if group is None or not len(group):
        return None
    n = len(group)

    def top_value(column):
        if column not in group:
            return "—", 0.0
        values = [
            str(token)
            for raw in group[column].dropna()
            for token in re.split(r"\s*\|\s*", str(raw))
            if token and token != "nan"
        ]
        values = [re.sub(r"\s*\(.*", "", value).strip() for value in values]
        if not values:
            return "—", 0.0
        counts = pd.Series(values).value_counts()
        return counts.index[0], round(100 * counts.iloc[0] / n, 1)

    top_pfam, top_pfam_frac = top_value("pfam_domains")
    top_pdb, top_pdb_frac = top_value("pdb_hit")
    top_name, _ = top_value("afdbsp_name")
    label = (
        top_name
        if top_name != "—"
        else top_pfam
        if top_pfam != "—"
        else top_pdb
        if top_pdb != "—"
        else "novel/unknown"
    )
    members = []
    for _, row in group.iterrows():
        novel = None if "novel" not in group or pd.isna(row.novel) else bool(row.novel)
        members.append({
            "acc": str(row.acc),
            "gene": (
                _clean_cell(row.gene_name)
                if "gene_name" in group
                else _clean_cell(row.protein_name)
                if "protein_name" in group
                else ""
            ),
            "novel": novel,
            "tm": int(row.n_TMR) if "n_TMR" in group and pd.notna(row.n_TMR) else 0,
            "eff": _clean_cell(row.effectorp) if "effectorp" in group else "",
            "pfam": _clean_cell(row.pfam_domains) if "pfam_domains" in group else "",
            "ipr": _clean_cell(row.interpro_entries) if "interpro_entries" in group else "",
            "pdb": _clean_cell(row.pdb_hit) if "pdb_hit" in group else "",
            "pdb_tm": float(row.pdb_tm) if "pdb_tm" in group and pd.notna(row.pdb_tm) else None,
            "afdb": (
                _clean_cell(row.afdbsp_name)
                if "afdbsp_name" in group
                else _clean_cell(row.afdbsp_hit)
                if "afdbsp_hit" in group
                else ""
            ),
            "afdb_hit": _clean_cell(row.afdbsp_hit) if "afdbsp_hit" in group else "",
            "afdb_tm": (
                float(row.afdbsp_tm)
                if "afdbsp_tm" in group and pd.notna(row.afdbsp_tm)
                else None
            ),
            "annotation_status": (
                _clean_cell(row.annotation_status) if "annotation_status" in group else ""
            ),
            "foldseek_pdb_status": (
                _clean_cell(row.foldseek_pdb_status)
                if "foldseek_pdb_status" in group
                else ""
            ),
            "foldseek_afdb_status": (
                _clean_cell(row.foldseek_afdb_status)
                if "foldseek_afdb_status" in group
                else ""
            ),
        })

    known_novel = group.novel.dropna() if "novel" in group else pd.Series(dtype=float)
    return {
        "label": label,
        "n": int(n),
        "pct_novel": round(100 * known_novel.mean(), 1) if len(known_novel) else None,
        "pct_eff": round(100 * group.is_effector.mean(), 1) if "is_effector" in group else None,
        "pct_domain": round(
            100 * group.has_any_domain.mean(), 1
        ) if "has_any_domain" in group else None,
        "top_pfam": top_pfam,
        "top_pfam_frac": top_pfam_frac,
        "top_pdb": top_pdb,
        "top_pdb_frac": top_pdb_frac,
        "top_ipr": top_value("interpro_entries")[0],
        "n_multi": int(group.multi_domain.sum()) if "multi_domain" in group else 0,
        "fusion": bool(group.multi_domain.sum()) if "multi_domain" in group else False,
        "members": members,
    }


def _overlapping_annotations(row, start, end):
    """Return coordinate-bearing annotations that overlap one domain segment."""
    annotations = []
    if row is None:
        return annotations
    for column, source in (
        ("pfam_domains", "Pfam"),
        ("interpro_entries", "InterPro"),
    ):
        raw = _clean_cell(row.get(column, ""))
        for token in re.split(r"\s*\|\s*", raw):
            token = token.strip()
            match = re.search(r"\((\d+)\s*-\s*(\d+)\)\s*$", token)
            if not match:
                continue
            ann_start, ann_end = sorted(map(int, match.groups()))
            overlap = max(0, min(int(end), ann_end) - max(int(start), ann_start) + 1)
            if overlap:
                annotations.append({
                    "source": source,
                    "label": token[:match.start()].strip(),
                    "start": ann_start,
                    "end": ann_end,
                    "overlap": overlap,
                })
    return annotations


def _read_fasta_records(path):
    """Read FASTA/MSA records without changing headers or aligned sequences."""
    records = {}
    if not path or not os.path.exists(path):
        return records
    header = None
    chunks = []
    with open(path, encoding="utf-8", errors="replace") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records[header] = "".join(chunks)
                header = line[1:].split()[0]
                chunks = []
            elif header is not None:
                chunks.append(line)
    if header is not None:
        records[header] = "".join(chunks)
    return records


def _records_by_member(records, members):
    """Map FoldMason headers (often strain-prefixed) back to member accessions."""
    mapped = {}
    for member in members:
        exact = records.get(member)
        if exact is not None:
            mapped[member] = exact
            continue
        hits = [seq for header, seq in records.items() if member in header]
        if len(hits) == 1:
            mapped[member] = hits[0]
    return mapped


def _ca_coordinates(pdbtext):
    coords = []
    for line in pdbtext.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                coords.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except (ValueError, IndexError):
                continue
    return np.asarray(coords, dtype=float)


def _aligned_ca_pairs(ref_pdb, mobile_pdb, ref_aln=None, mobile_aln=None):
    """Return corresponding CA coordinates, preferably from the FoldMason MSA."""
    ref_ca = _ca_coordinates(ref_pdb)
    mob_ca = _ca_coordinates(mobile_pdb)
    if ref_aln and mobile_aln and len(ref_aln) == len(mobile_aln):
        pairs = []
        ri = mi = 0
        for ra, ma in zip(ref_aln, mobile_aln):
            rpos = ri if ra not in "-." else None
            mpos = mi if ma not in "-." else None
            if rpos is not None:
                ri += 1
            if mpos is not None:
                mi += 1
            if rpos is not None and mpos is not None and rpos < len(ref_ca) and mpos < len(mob_ca):
                pairs.append((rpos, mpos))
        if len(pairs) >= 3:
            return (ref_ca[[p[0] for p in pairs]], mob_ca[[p[1] for p in pairs]], "foldmason")
    n = min(len(ref_ca), len(mob_ca))
    if n < 3:
        raise ValueError("at least three paired CA atoms are required for superposition")
    return ref_ca[:n], mob_ca[:n], "ca_order"


def _superpose_pdb(mobile_pdb, ref_pdb, mobile_aln=None, ref_aln=None):
    """Rigidly align a PDB to a reference and return transformed text plus fit metadata."""
    ref_xyz, mob_xyz, method = _aligned_ca_pairs(ref_pdb, mobile_pdb, ref_aln, mobile_aln)
    def fit(mask):
        ref_fit = ref_xyz[mask]
        mob_fit = mob_xyz[mask]
        ref_center = ref_fit.mean(axis=0)
        mob_center = mob_fit.mean(axis=0)
        u, _, vt = np.linalg.svd((mob_fit - mob_center).T @ (ref_fit - ref_center))
        rotation = u @ vt
        if np.linalg.det(rotation) < 0:
            u[:, -1] *= -1
            rotation = u @ vt
        return rotation, ref_center - mob_center @ rotation

    # FoldMason columns provide correspondence, while iterative rejection prevents long
    # flexible loops from pulling the conserved core away from the hub. Four angstroms
    # is deliberately permissive; if an initial fit has too few inliers, retain the
    # closest half and iterate rather than failing a divergent but valid family.
    mask = np.ones(len(ref_xyz), dtype=bool)
    for _ in range(8):
        rot, tran = fit(mask)
        distances = np.linalg.norm(mob_xyz @ rot + tran - ref_xyz, axis=1)
        new_mask = distances <= 4.0
        if new_mask.sum() < 3:
            keep = max(3, int(math.ceil(len(ref_xyz) * 0.5)))
            new_mask = np.zeros(len(ref_xyz), dtype=bool)
            new_mask[np.argsort(distances)[:keep]] = True
        if np.array_equal(new_mask, mask):
            break
        mask = new_mask
    rot, tran = fit(mask)
    fitted = mob_xyz @ rot + tran
    squared = np.sum((fitted - ref_xyz) ** 2, axis=1)
    rmsd = float(np.sqrt(np.mean(squared[mask])))
    rmsd_all = float(np.sqrt(np.mean(squared)))
    output = []
    for line in mobile_pdb.splitlines():
        if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
            try:
                xyz = np.asarray([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            except ValueError:
                output.append(line)
                continue
            x, y, z = xyz @ rot + tran
            line = f"{line[:30]}{x:8.3f}{y:8.3f}{z:8.3f}{line[54:]}"
        output.append(line)
    return "\n".join(output) + "\n", {
        "method": method, "n_ca": int(mask.sum()), "n_ca_total": len(ref_xyz),
        "rmsd": round(rmsd, 4), "rmsd_all": round(rmsd_all, 4),
        "rotation": rot.round(10).tolist(), "translation": tran.round(10).tolist(),
    }


def _structures_zip_b64(fam, structures):
    """Create a family ZIP containing one independently usable PDB per member."""
    if not structures:
        return ""
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for member, pdbtext in structures.items():
            archive.writestr(f"{fam}_structures/{member}.pdb", pdbtext.rstrip() + "\n")
        manifest = "family\tmember\tfile\n" + "".join(
            f"{fam}\t{member}\t{member}.pdb\n" for member in structures)
        archive.writestr(f"{fam}_structures/manifest.tsv", manifest)
    return base64.b64encode(buf.getvalue()).decode()


def _domain_structures_zip_b64(fam, members, parent_structures, segments_only):
    """Create a ZIP of cropped domain members or their complete parent proteins."""
    buf = io.BytesIO()
    written = set()
    with zipfile.ZipFile(buf, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        manifest_rows = []
        for member in members:
            accession = str(member["acc"])
            parent = parent_structures.get(accession, "")
            if not parent:
                continue
            if segments_only:
                name = re.sub(r"[^A-Za-z0-9_.-]+", "_", member["segment_id"])
                text = "\n".join(
                    line for line in parent.splitlines()
                    if line.startswith(("ATOM", "HETATM"))
                    and line[22:26].strip().lstrip("-").isdigit()
                    and int(line[22:26]) >= int(member["start"])
                    and int(line[22:26]) <= int(member["end"])
                )
                filename = f"{name}.pdb"
                archive.writestr(
                    f"{fam}_domain_structures/{filename}",
                    text.rstrip() + "\nTER\nEND\n",
                )
                manifest_rows.append(
                    {
                        "segment_id": member["segment_id"],
                        "accession": accession,
                        "start": member["start"],
                        "end": member["end"],
                        "file": filename,
                    }
                )
            elif accession not in written:
                written.add(accession)
                filename = f"{accession}.pdb"
                archive.writestr(
                    f"{fam}_parent_structures/{filename}",
                    parent.rstrip() + "\n",
                )
                manifest_rows.append(
                    {
                        "segment_id": "",
                        "accession": accession,
                        "start": "",
                        "end": "",
                        "file": filename,
                    }
                )
        manifest = pd.DataFrame(
            manifest_rows,
            columns=["segment_id", "accession", "start", "end", "file"],
        ).to_csv(index=False)
        root = (
            f"{fam}_domain_structures"
            if segments_only
            else f"{fam}_parent_structures"
        )
        archive.writestr(f"{root}/manifest.csv", manifest)
    return base64.b64encode(buf.getvalue()).decode()


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


def _p2rank_prediction_csv(results_dir, fam, ref):
    """Return the P2Rank table for the current reference, never a stale hub."""
    candidates = sorted(glob.glob(
        os.path.join(results_dir, "p2rank", fam, "out", "*_predictions.csv")
    ))
    if not candidates:
        return None
    matches = [
        path for path in candidates
        if ref and (
            os.path.basename(path).endswith(f"{ref}.pdb_predictions.csv")
            or os.path.basename(path).endswith(f"{ref}_predictions.csv")
        )
    ]
    if len(matches) == 1:
        return matches[0]
    if len(matches) > 1:
        raise RuntimeError(
            f"{fam}: multiple P2Rank prediction tables match reference {ref}"
        )
    # Older single-file runs did not consistently preserve the source accession
    # in the output name. They are safe only when there is no competing stale file.
    return candidates[0] if len(candidates) == 1 else None


def _family_expression(expression_all, members):
    """Subset the run-level RNA-seq table for one full-length family."""
    if (
        expression_all is None
        or "acc" not in expression_all.columns
        or not members
    ):
        return None
    subset = expression_all[
        expression_all["acc"].astype(str).isin(set(map(str, members)))
    ].copy()
    return subset if len(subset) else None


def _esm_tolerance(esm_all, key):
    """Return residue-indexed mean ESM LLR values for one protein/family key."""
    if esm_all is None or "family" not in esm_all.columns:
        return {}
    rows = esm_all[esm_all.family.astype(str) == str(key)]
    if not len(rows):
        return {}
    aa_columns = [
        column for column in rows.columns
        if re.fullmatch(r"[A-Z]", str(column))
    ]
    if not aa_columns:
        return {}
    position_column = next(
        (
            column for column in rows.columns
            if str(column).lower() in ("", "unnamed: 0", "pos", "site")
        ),
        rows.columns[0],
    )
    positions = pd.to_numeric(
        rows[position_column].astype(str).str.extract(r"(\d+)$")[0],
        errors="coerce",
    )
    means = rows[aa_columns].apply(pd.to_numeric, errors="coerce").mean(axis=1)
    return {
        int(position): float(value)
        for position, value in zip(positions, means)
        if pd.notna(position) and pd.notna(value) and math.isfinite(float(value))
    }


def _enrich_pocket_entry(results_dir, fam, entry):
    """Backfill all pocket predictions from raw outputs for pre-v1.0.2 runs."""
    entry = json.loads(json.dumps(entry or {}))
    ref = entry.get("ref", "")
    storage_key = entry.get("storage_key", fam)

    p2 = entry.get("p2rank", {}) or {}
    if not p2.get("pockets"):
        prediction_csv = _p2rank_prediction_csv(results_dir, storage_key, ref)
        if prediction_csv:
            table = pd.read_csv(prediction_csv)
            table.columns = [str(c).strip() for c in table.columns]
            predictions = []
            for idx, pred in table.iterrows():
                tokens = str(pred.get("residue_ids", "")).split()
                residues = sorted({int(x.split("_")[-1]) for x in tokens
                                   if x.split("_")[-1].isdigit()})
                predictions.append({
                    "pocket_id": int(pred.get("rank", idx + 1)),
                    "score": float(pred.get("score", 0)),
                    "probability": (
                        float(pred.get("probability"))
                        if pd.notna(pred.get("probability"))
                        else None
                    ),
                    "lining_residues": residues,
                })
            if predictions:
                top = max(predictions, key=lambda p: p["score"])
                p2.update(top_score=top["score"],
                          top_probability=top["probability"],
                          n_pockets=len(predictions),
                          lining_residues=top["lining_residues"], pockets=predictions)
                entry["p2rank"] = p2

    fp = entry.get("fpocket", {}) or {}
    if ref and not fp.get("pockets"):
        root = os.path.join(results_dir, "fpocket", storage_key, f"{ref}_out")
        info = os.path.join(root, f"{ref}_info.txt")
        if os.path.exists(info):
            text = open(info, encoding="utf-8", errors="replace").read()
            scores = {int(n): float(score) for n, score in
                      re.findall(r"Pocket\s+(\d+)\s*:\s*\n\s*Score\s*:\s*([-\d.]+)", text)}
            predictions = []
            for pocket_id, score in sorted(scores.items()):
                residues = set()
                atom_file = os.path.join(root, "pockets", f"pocket{pocket_id}_atm.pdb")
                if os.path.exists(atom_file):
                    for line in open(atom_file, encoding="utf-8", errors="replace"):
                        if line.startswith(("ATOM", "HETATM")):
                            try:
                                residues.add(int(line[22:26]))
                            except (ValueError, IndexError):
                                pass
                predictions.append({"pocket_id": pocket_id, "score": score,
                                    "lining_residues": sorted(residues)})
            if predictions:
                top = max(predictions, key=lambda p: p["score"])
                fp.update(top_score=top["score"], n_pockets=len(predictions),
                          lining_residues=top["lining_residues"], pockets=predictions)
                entry["fpocket"] = fp
    return entry


def _pocket_raw_tables(results_dir, fam, entry):
    """Load detector-native pocket tables for lossless workbook export."""
    tables = {}
    ref = (entry or {}).get("ref", "")
    storage_key = (entry or {}).get("storage_key", fam)
    prediction_csv = _p2rank_prediction_csv(results_dir, storage_key, ref)
    if prediction_csv:
        p2_table = pd.read_csv(prediction_csv)
        p2_table.columns = [str(c).strip() for c in p2_table.columns]
        tables["p2rank_pockets"] = p2_table

    info = os.path.join(
        results_dir,
        "fpocket",
        storage_key,
        f"{ref}_out",
        f"{ref}_info.txt",
    ) if ref else ""
    if info and os.path.exists(info):
        text = open(info, encoding="utf-8", errors="replace").read()
        parts = re.split(r"Pocket\s+(\d+)\s*:\s*", text)
        rows = []
        for idx in range(1, len(parts), 2):
            row = {"pocket_id": int(parts[idx])}
            for key, value in re.findall(r"^\s*([^:\n]+?)\s*:\s*([^\n]*)$", parts[idx + 1], re.MULTILINE):
                column = re.sub(r"[^a-z0-9]+", "_", key.strip().lower()).strip("_")
                row[column] = value.strip()
            rows.append(row)
        if rows:
            tables["fpocket_pockets"] = pd.DataFrame(rows)
    return tables


def _read_tpl(name):
    with open(os.path.join(_TPL, name), encoding="utf-8") as fh:
        return fh.read()


def _svg_datauri(svg):
    """The v19 renderer injects every asset into <img src="...">, so an asset MUST be a
    data URI, not raw <svg> markup (raw markup renders blank). Wrap as base64 data URI."""
    if not svg:
        return ""
    if svg.startswith("data:"):
        return svg
    return "data:image/svg+xml;base64," + base64.b64encode(svg.encode("utf-8")).decode()


def _svg_matrix(mat, labels, title, cmap_lo=(240, 245, 250), cmap_hi=(20, 90, 140),
                vmin=0, vmax=1, unit="TM"):
    """Labelled numeric matrix as inline SVG. Every cell shows its VALUE (×100, 0–100 scale
    so 2 chars fit) on a colour ramp, with a colour-bar legend mapping shade→value. Vector,
    so it stays crisp when zoomed."""
    n = len(labels)
    if n == 0:
        return ""
    # cell big enough to hold a 2–3 char number; large families get a smaller but still
    # numbered cell (SVG is vector — zoom keeps it readable)
    cell = max(20, min(40, int(560 / max(n, 1))))
    fs = max(7, min(13, int(cell * 0.42)))
    labfs = max(7, fs - 1)
    # vertical room the -60° rotated column labels actually occupy above the grid, from the
    # LONGEST label — otherwise they shoot up into the title. sin(60)=0.866, ~0.55px per char.
    maxlab = max((len(str(l)) for l in labels), default=8)
    lab_room = int(maxlab * labfs * 0.55 * 0.866) + 10
    title_band = 40                       # title (y16) + subtitle (y30)
    pad_t = title_band + lab_room         # grid top
    pad_l = max(96, int(maxlab * labfs * 0.62) + 12)   # room for row labels on the left
    pad = pad_t                           # kept for the legend math below
    barw = 54
    W = pad_l + n * cell + barw + 34
    H = pad_t + n * cell + 20
    def col(v):
        t = 0.0 if vmax == vmin else max(0.0, min(1.0, (v - vmin) / (vmax - vmin)))
        r = int(cmap_lo[0] + (cmap_hi[0] - cmap_lo[0]) * t)
        g = int(cmap_lo[1] + (cmap_hi[1] - cmap_lo[1]) * t)
        b = int(cmap_lo[2] + (cmap_hi[2] - cmap_lo[2]) * t)
        return f"#{r:02x}{g:02x}{b:02x}"
    def txtcol(v):
        t = 0.0 if vmax == vmin else max(0.0, min(1.0, (v - vmin) / (vmax - vmin)))
        return "#fff" if t > 0.55 else "#1a2b3c"   # contrast against cell fill
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
             f'font-family="sans-serif" font-size="{fs}">',
             f'<text x="6" y="16" font-size="13" font-weight="600">{title}</text>',
             f'<text x="6" y="30" font-size="9" fill="#667">values ×100 (e.g. 59 = {unit} 0.59); diagonal = self</text>']
    for i in range(n):
        for j in range(n):
            v = mat[i][j]
            x = pad_l + j * cell
            y = pad_t + i * cell
            if v is None or (isinstance(v, float) and math.isnan(v)):
                parts.append(f'<rect x="{x}" y="{y}" width="{cell-1}" height="{cell-1}" fill="#eee"/>')
                continue
            parts.append(f'<rect x="{x}" y="{y}" width="{cell-1}" height="{cell-1}" fill="{col(v)}"/>')
            iv = int(round(v * 100))
            parts.append(f'<text x="{x+cell/2-0.5}" y="{y+cell/2+fs*0.35}" text-anchor="middle" '
                         f'fill="{txtcol(v)}" font-size="{fs}">{iv}</text>')
    for i, lab in enumerate(labels):
        y = pad_t + i * cell + cell * 0.62
        parts.append(f'<text x="{pad_l-4}" y="{y}" text-anchor="end" font-size="{labfs}">{lab}</text>')
        x = pad_l + i * cell + cell * 0.5
        # text-anchor="start" so the label BEGINS at the grid edge and extends up-right,
        # sitting ENTIRELY above the grid (middle-anchor dropped half the label into row 1)
        parts.append(f'<text x="{x}" y="{pad_t-4}" text-anchor="start" font-size="{labfs}" '
                     f'transform="rotate(-60 {x} {pad_t-4})">{lab}</text>')
    # colour-bar legend (vertical), mapping shade -> value
    bx = pad_l + n * cell + 18
    by = pad_t
    bh = min(n * cell, 220)
    steps = 40
    for s in range(steps):
        t = s / (steps - 1)
        val = vmin + (vmax - vmin) * (1 - t)   # top = high
        yy = by + bh * t
        parts.append(f'<rect x="{bx}" y="{yy:.1f}" width="16" height="{bh/steps+0.6:.1f}" fill="{col(val)}"/>')
    for frac in (0, 0.25, 0.5, 0.75, 1.0):
        val = vmin + (vmax - vmin) * (1 - frac)
        yy = by + bh * frac
        parts.append(f'<line x1="{bx+16}" y1="{yy:.1f}" x2="{bx+20}" y2="{yy:.1f}" stroke="#333"/>')
        parts.append(f'<text x="{bx+23}" y="{yy+3:.1f}" font-size="9">{int(round(val*100))}</text>')
    parts.append(f'<text x="{bx+2}" y="{by-6}" font-size="9" font-weight="600">{unit}×100</text>')
    parts.append("</svg>")
    return "".join(parts)


def _svg_heat(df, title):
    """RNAseq expression heatmap (members x conditions). Colour = row z-score (so each gene's
    pattern across conditions is visible), but every cell PRINTS its real expression value
    (raw counts) so the user reads actual numbers, plus a z-score colour legend."""
    num = df.select_dtypes("number")
    if num.empty:
        return ""
    raw = num.values.astype(float)
    lv = np.log1p(np.clip(raw, 0, None))
    mu = lv.mean(1, keepdims=True); sd = lv.std(1, keepdims=True); sd[sd == 0] = 1
    z = (lv - mu) / sd
    rows, cols = z.shape
    cw = 62
    ch = max(16, min(26, int(360 / max(rows, 1))))
    fs = max(8, min(11, int(ch * 0.5)))
    labfs = max(8, fs - 1)
    # top padding must clear the -60° rotated condition labels (which sit above the grid):
    # vertical reach of the longest label = chars * ~0.55px/char * sin(60°)
    maxcond = max((len(str(c)) for c in num.columns), default=6)
    pad_l = 130
    pad_t = 40 + int(maxcond * labfs * 0.55 * 0.866) + 10
    barh = 14
    W = pad_l + cols * cw + 30
    H = pad_t + rows * ch + 40
    def col(v):
        t = max(-2.5, min(2.5, v)) / 2.5
        if t >= 0:
            r, g, b = 200 - int(120*t), 60 + int(40*t), 40
        else:
            r, g, b = 40, 80 - int(30*t), 150 + int(80*(-t))
        return f"#{max(0,min(255,r)):02x}{max(0,min(255,g)):02x}{max(0,min(255,b)):02x}"
    def fmt(v):
        if v >= 1000: return f"{v/1000:.1f}k"
        if v >= 10:   return f"{v:.0f}"
        if v >= 1:    return f"{v:.1f}"
        return f"{v:.1f}"
    P = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'font-family="sans-serif" font-size="{fs}">',
         f'<text x="6" y="16" font-size="13" font-weight="600">{title}</text>',
         f'<text x="6" y="30" font-size="9" fill="#667">cell = mean raw count; colour = per-gene z-score across conditions</text>']
    for i in range(rows):
        for j in range(cols):
            x = pad_l + j*cw; y = pad_t + i*ch
            P.append(f'<rect x="{x}" y="{y}" width="{cw-1}" height="{ch-1}" fill="{col(z[i][j])}"/>')
            tc = "#fff" if abs(z[i][j]) > 1.3 else "#1a2b3c"
            P.append(f'<text x="{x+cw/2}" y="{y+ch/2+fs*0.35}" text-anchor="middle" '
                     f'fill="{tc}" font-size="{fs}">{fmt(raw[i][j])}</text>')
    for i, lab in enumerate(df.index.astype(str)):
        P.append(f'<text x="{pad_l-4}" y="{pad_t+i*ch+ch*0.62}" text-anchor="end" font-size="{max(8,fs-1)}">{lab}</text>')
    for j, lab in enumerate(num.columns.astype(str)):
        x = pad_l + j*cw + cw*0.5
        # text-anchor="start" so the rotated label begins at the grid edge and extends
        # up-right, sitting ENTIRELY above the grid (middle-anchor dropped half into row 1)
        P.append(f'<text x="{x}" y="{pad_t-4}" text-anchor="start" font-size="{labfs}" '
                 f'transform="rotate(-60 {x} {pad_t-4})">{lab}</text>')
    # z-score colour legend (horizontal) under the grid
    ly = pad_t + rows * ch + 16
    lx = pad_l
    seg = (cols * cw) if cols else 200
    steps = 40
    for s in range(steps):
        t = s / (steps - 1)
        zz = -2.5 + 5.0 * t
        P.append(f'<rect x="{lx + seg*t:.1f}" y="{ly}" width="{seg/steps+0.6:.1f}" height="{barh}" fill="{col(zz)}"/>')
    for frac, lab in ((0, "-2.5"), (0.5, "0"), (1.0, "+2.5")):
        P.append(f'<text x="{lx + seg*frac:.1f}" y="{ly+barh+11}" text-anchor="middle" font-size="9">{lab}</text>')
    P.append(f'<text x="{lx-4}" y="{ly+barh-2}" text-anchor="end" font-size="9" font-weight="600">z</text>')
    P.append("</svg>")
    return "".join(P)


def _newick_to_svg(nwk, hub=None, title="FoldTree"):
    """Minimal rectangular cladogram with an explicit evidence-source title."""
    try:
        import io as _io
        from Bio import Phylo
        # MAD rooting can emit MULTIPLE candidate-rooted trees (one per line) when the root is
        # ambiguous; Phylo.read expects exactly one, so take the first newick string only.
        first = next((ln for ln in nwk.splitlines() if ln.strip()), nwk).strip()
        tree = Phylo.read(_io.StringIO(first), "newick")
    except Exception:
        return ""
    leaves = tree.get_terminals()
    n = len(leaves)
    if n == 0:
        return ""
    ypos = {lf: i for i, lf in enumerate(leaves)}
    # real branch lengths (structural distance) if present, else topological depth
    has_bl = any(getattr(cl, "branch_length", None) for cl in tree.find_clades())
    depths = tree.depths() if has_bl else tree.depths(unit_branch_lengths=True)
    maxd = max(depths.values()) or 1
    ch = max(10, min(22, int(320 / n))); H = 40 + n*ch
    maxlab = max((len(lf.name or "") for lf in leaves), default=10)
    labw = int(maxlab * 6.0) + 16          # room for the longest accession (+★)
    xs = 20; xe = 20 + 300                 # tree drawing area
    W = xe + 8 + labw                      # viewport includes the label column
    def X(cl): return xs + (xe-xs) * depths.get(cl, 0) / maxd
    def Y(cl):
        if cl.is_terminal(): return 20 + ypos[cl]*ch
        kids = cl.clades
        return (Y(kids[0]) + Y(kids[-1]))/2
    P = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
         f'font-family="sans-serif" font-size="10">',
         f'<text x="6" y="14" font-size="12" font-weight="600">{title}</text>']
    def draw(cl):
        x0 = X(cl); y0 = Y(cl)
        for k in cl.clades:
            x1 = X(k); y1 = Y(k)
            P.append(f'<path d="M{x0},{y0} L{x0},{y1} L{x1},{y1}" fill="none" stroke="#456" stroke-width="1"/>')
            draw(k)
    draw(tree.root)
    for lf in leaves:
        nm = lf.name or ""
        y = 20 + ypos[lf]*ch
        star = ' ★' if hub and nm == hub else ''
        colr = '#c79a00' if star else '#233'
        P.append(f'<text x="{xe+3}" y="{y+3}" font-size="9" fill="{colr}">{nm}{star}</text>')
    if hub:
        P.append(f'<text x="6" y="{H-6}" font-size="9" fill="#c79a00">★ hub = highest mean structural similarity to family members</text>')
    P.append("</svg>")
    return "".join(P)


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


def _hub_from_tm(tm, labels):
    """Hub = member with highest mean off-diagonal TM ('most like everyone')."""
    try:
        M = np.array(tm.set_index(tm.columns[0]).values, dtype=float)
        n = M.shape[0]
        if n < 2:
            return (labels[0] if labels else None), None
        off = (M.sum(1) - np.diag(M)) / (n - 1)
        i = int(np.argmax(off))
        return labels[i], round(float(off[i]), 3)
    except Exception:
        return (labels[0] if labels else None), None


def _matrix_pair_stats(table):
    """Return explicit off-diagonal stats for a labeled symmetric matrix."""
    empty = {
        "mean_all": None,
        "mean_detected": None,
        "maximum": None,
        "n_pairs": 0,
        "n_detected": 0,
    }
    if table is None or len(table) < 2:
        return empty
    try:
        matrix = table.set_index(table.columns[0]).apply(
            pd.to_numeric, errors="coerce"
        ).to_numpy(dtype=float)
        values = matrix[np.triu_indices(len(matrix), k=1)]
        values = values[np.isfinite(values)]
        if not len(values):
            return empty
        detected = values[values > 0]
        return {
            "mean_all": round(float(values.mean()), 3),
            "mean_detected": (
                round(float(detected.mean()), 3) if len(detected) else None
            ),
            "maximum": round(float(values.max()), 3),
            "n_pairs": int(len(values)),
            "n_detected": int(len(detected)),
        }
    except Exception:
        return empty


def build_atlas(master_csv, cards_dir, composition_xlsx, annotation_csv,
                results_dir, out_html, mode="single", atlas_name="atlas", config=None):
    config = config or {}
    famdir = os.path.join(results_dir, "families")
    master = pd.read_csv(master_csv)
    anno = pd.read_csv(annotation_csv) if os.path.exists(annotation_csv) else pd.DataFrame()
    downloads_dir = os.path.join(results_dir, "downloads")
    if mode == "backend":
        os.makedirs(downloads_dir, exist_ok=True)

    def store_download(name, encoded):
        if mode != "backend":
            return encoded
        path = os.path.join(downloads_dir, name)
        with open(path, "wb") as handle:
            handle.write(base64.b64decode(encoded))
        return None

    def load_csv(p):
        return pd.read_csv(p) if os.path.exists(p) else None

    members_all = load_csv(os.path.join(results_dir, "members.csv"))
    expression_all = load_csv(os.path.join(results_dir, "rnaseq_expression.csv"))
    sequence_subgroups_all = load_csv(
        os.path.join(results_dir, "sequence_subgroups.csv")
    )
    domain_members_all = load_csv(os.path.join(results_dir, "domain_members.csv"))
    domain_families_all = load_csv(os.path.join(results_dir, "domain_families.csv"))
    domain_edges_all = load_csv(os.path.join(results_dir, "domain_edges.csv"))
    domain_bridges_all = load_csv(os.path.join(results_dir, "domain_cross_edges.csv"))
    domain_workbench = {"schema_version": 1, "families": {}}
    domain_workbench_path = os.path.join(results_dir, "domain_workbench.json")
    if os.path.exists(domain_workbench_path):
        try:
            domain_workbench = json.load(open(domain_workbench_path))
        except Exception:
            domain_workbench = {"schema_version": 1, "families": {}}
    for workbench in domain_workbench.get("families", {}).values():
        label_map = workbench.get("tree_label_map", {})
        def display_tree(newick):
            for safe_id, segment_id in sorted(
                label_map.items(), key=lambda item: -len(item[0])
            ):
                display_id = str(segment_id).replace(":", "_")
                newick = re.sub(
                    rf"(?<=[(,]){re.escape(str(safe_id))}(?=[:),;])",
                    display_id,
                    newick,
                )
            return newick

        tree = workbench.get("sequence_newick", "")
        if tree:
            workbench["sequence_tree_svg"] = _svg_datauri(
                _newick_to_svg(
                    display_tree(tree),
                    title="Sequence tree · MAFFT + FastTree",
                )
            )
        for subgroup in workbench.get("sequence_subgroups", []):
            if subgroup.get("newick"):
                subgroup["tree_svg"] = _svg_datauri(
                    _newick_to_svg(
                        display_tree(subgroup["newick"]),
                        title=(
                            f"Sequence tree · {subgroup.get('id', '')} · "
                            "MAFFT + FastTree"
                        ),
                    )
                )
        foldmason_tree = workbench.get("foldmason_guide_newick", "")
        if foldmason_tree:
            workbench["foldmason_guide_tree_svg"] = _svg_datauri(
                _newick_to_svg(
                    display_tree(foldmason_tree),
                    title="FoldMason structural guide tree",
                )
            )
        foldtree_svgs = {}
        for metric, newick in workbench.get("foldtree_trees", {}).items():
            foldtree_svgs[metric] = _svg_datauri(
                _newick_to_svg(
                    display_tree(newick),
                    title=f"FoldTree structural tree · {metric}",
                )
            )
        workbench["foldtree_tree_svgs"] = foldtree_svgs
        labels = workbench.get("usalign_labels", [])
        matrix = workbench.get("usalign_matrix", [])
        if labels and matrix:
            workbench["usalign_matrix_svg"] = _svg_datauri(
                _svg_matrix(
                    matrix,
                    labels,
                    "Domain-segment US-align TM-score",
                    unit="TM",
                )
            )
        sequence_identity_labels = workbench.get(
            "sequence_identity_labels", []
        )
        sequence_identity_matrix = workbench.get(
            "sequence_identity_matrix", []
        )
        if sequence_identity_labels and sequence_identity_matrix:
            workbench["sequence_identity_matrix_svg"] = _svg_datauri(
                _svg_matrix(
                    sequence_identity_matrix,
                    sequence_identity_labels,
                    "Domain-segment BLASTp identity",
                    vmin=0,
                    vmax=1.0,
                    unit="%id",
                )
            )
        structural_identity_labels = workbench.get(
            "structural_alignment_identity_labels", []
        )
        structural_identity_matrix = workbench.get(
            "structural_alignment_identity_matrix", []
        )
        if structural_identity_labels and structural_identity_matrix:
            workbench["structural_alignment_identity_matrix_svg"] = (
                _svg_datauri(
                    _svg_matrix(
                        structural_identity_matrix,
                        structural_identity_labels,
                        "FoldMason-aligned amino-acid identity",
                        vmin=0,
                        vmax=1.0,
                        unit="%id",
                    )
                )
            )

    def structure_text(accession, family_dir=None):
        candidates = []
        if family_dir:
            candidates.append(os.path.join(family_dir, f"{accession}.pdb"))
        input_dir = str(config.get("input", {}).get("pdb_dir", "") or "")
        strain = str(config.get("strain", {}).get("code", "") or "")
        if input_dir:
            candidates.extend([
                os.path.join(input_dir, f"{strain}_{accession}.pdb"),
                os.path.join(input_dir, f"{accession}.pdb"),
            ])
        candidates.append(os.path.join(results_dir, "..", "input", "pdb",
                                       f"{strain}_{accession}.pdb"))
        for candidate in candidates:
            if os.path.exists(candidate):
                return open(candidate, encoding="utf-8", errors="replace").read()
        return ""

    # whole-set outputs loaded once (may be absent if that step was toggled off)
    pockets = {}
    pj = os.path.join(results_dir, "pockets.json")
    if os.path.exists(pj):
        try: pockets = json.load(open(pj))
        except Exception: pockets = {}
    esm_all = load_csv(os.path.join(results_dir, "esm_all.csv"))   # long: family,resi,wt,<AA LLRs>
    classification_all = load_csv(os.path.join(results_dir, "classification.csv"))
    # whole-set mature sequences (seqs rule) -> {acc: seq}; used for per-member FASTA download
    seqs_all = {}
    sf = os.path.join(results_dir, "seqs.fasta")
    if os.path.exists(sf):
        acc = None; buf = []
        for ln in open(sf, encoding="utf-8", errors="replace"):
            ln = ln.rstrip("\n")
            if ln.startswith(">"):
                if acc: seqs_all[acc] = "".join(buf)
                acc = ln[1:].split()[0]; buf = []
            else:
                buf.append(ln.strip())
        if acc: seqs_all[acc] = "".join(buf)

    NET_nodes, PAY, EXTRA, REFPDB, REFAVAIL, ANN = [], {}, {}, {}, {}, {}
    esm_by_fam = {}   # fam -> {ref, tol{resi:score}} for building REFPDB["<fam>_esm"]

    def store_reference(key, text):
        if not text:
            return
        if mode == "backend":
            with open(
                os.path.join(downloads_dir, f"{key}.pdb"),
                "w",
                encoding="utf-8",
            ) as handle:
                handle.write(text.rstrip() + "\n")
            REFAVAIL[key] = True
        else:
            REFPDB[key] = text
    fam_accs = {}
    for _, r in master.iterrows():
        fam = r.family
        fd = os.path.join(famdir, fam)
        mem_file = os.path.join(
            results_dir, "family_members", f"{fam}.members.txt"
        )
        if not os.path.exists(mem_file):
            # Backward compatibility for atlases produced before v5.0.1.
            mem_file = os.path.join(
                results_dir, "families", f"{fam}.members.txt"
            )
        members = []
        if os.path.exists(mem_file):
            for line in open(mem_file):
                m = re.search(r"[A-Z]{2,3}\d{4,}\.\d+", line)
                if m: members.append(m.group(0))
        fam_accs[fam] = members
        # assets
        tm = load_csv(os.path.join(fd, f"{fam}_TM.csv"))
        idm = load_csv(os.path.join(fd, f"{fam}_ID.csv"))
        exp = load_csv(os.path.join(fd, f"{fam}_expression.csv"))
        if exp is None:
            exp = _family_expression(expression_all, members)
        tm_stats = _matrix_pair_stats(tm)
        id_stats = _matrix_pair_stats(idm)
        retained_tm = _finite_float(
            r.get("mean_retained_edge_TM", r.get("mean_TM"))
        ) or 0.0
        retained_fident = _finite_float(
            r.get(
                "mean_retained_edge_foldseek_fident",
                r.get("mean_identity"),
            )
        ) or 0.0
        NET_nodes.append(dict(
            id=fam,
            n=int(r.n_members),
            tm=tm_stats["mean_all"] or 0.0,
            retained_tm=retained_tm,
            id_pct=id_stats["mean_all"] or 0.0,
            id_detected=id_stats["mean_detected"],
            maxid=id_stats["maximum"] or 0.0,
            retained_foldseek_fident=retained_fident,
            n_pairs=tm_stats["n_pairs"],
            n_blast_pairs=id_stats["n_detected"],
            suss=float(r.get("suss_pct", 0) or 0),
            plddt=float(r.get("mean_pLDDT", 0) or 0),
            len=float(r.get("mean_len", 0) or 0),
        ))
        assets = {}
        if tm is not None:
            labs = list(tm.columns[1:]) if tm.columns[0].lower() in ("", "unnamed: 0") else list(tm.iloc[:, 0].astype(str))
            M = tm.set_index(tm.columns[0]).values.tolist()
            assets["tm_svg"] = _svg_datauri(_svg_matrix(M, [l[:11] for l in labs],
                                           f"{fam} · Structural similarity (Foldseek TM-score)",
                                           vmin=0.3, vmax=1.0, unit="TM"))
        # US-align TM matrix — algorithm-independent cross-check of the Foldseek TM above
        usm = load_csv(os.path.join(fd, f"{fam}_TM_usalign.csv"))
        tm_us_mean = tm_cons_r = tm_cons_maxdiff = None; tm_disagree = 0
        if usm is not None:
            ulabs = list(usm.iloc[:, 0].astype(str))
            UM = usm.set_index(usm.columns[0]).values.tolist()
            assets["tmus_svg"] = _svg_datauri(_svg_matrix(UM, [l[:11] for l in ulabs],
                                           f"{fam} · Structural similarity (US-align TM, independent)",
                                           vmin=0.3, vmax=1.0, unit="TM"))
            # Foldseek vs US-align consistency on the shared off-diagonal pairs
            if tm is not None:
                fsq = tm.set_index(tm.columns[0]); usq = usm.set_index(usm.columns[0])
                shared = [c for c in fsq.columns if c in usq.columns and c in fsq.index and c in usq.index]
                if len(shared) >= 2:
                    fa = fsq.loc[shared, shared].values; ua = usq.loc[shared, shared].values
                    iu = np.triu_indices(len(shared), k=1)
                    av, bv = fa[iu].astype(float), ua[iu].astype(float)
                    if len(av) >= 2 and np.std(av) > 0 and np.std(bv) > 0:
                        tm_cons_r = round(float(np.corrcoef(av, bv)[0, 1]), 3)
                    tm_us_mean = round(float(bv.mean()), 3)
                    tm_cons_maxdiff = round(float(np.abs(av - bv).max()), 3)
                    tm_disagree = int((np.abs(av - bv) > 0.1).sum())
        if idm is not None:
            labs = list(idm.iloc[:, 0].astype(str))
            M = idm.set_index(idm.columns[0]).values.tolist()
            assets["id_svg"] = _svg_datauri(_svg_matrix(M, [l[:11] for l in labs],
                                           f"{fam} · Sequence identity (BLASTp %)",
                                           cmap_hi=(150, 60, 20), vmin=0, vmax=1.0, unit="%id"))
        if exp is not None:
            assets["rna_svg"] = _svg_datauri(_svg_heat(exp.set_index(exp.columns[0]), f"{fam} · RNAseq"))
        sig = load_csv(os.path.join(fd, f"{fam}_signature.csv"))
        structural_cons = load_csv(
            os.path.join(fd, f"{fam}_structural_conservation.csv")
        )
        sequence_status = {}
        sequence_status_path = os.path.join(
            fd, f"{fam}_sequence_analysis_status.json"
        )
        if os.path.exists(sequence_status_path):
            try:
                sequence_status = json.load(open(sequence_status_path))
            except Exception:
                sequence_status = {}
        # hub = highest mean-TM member (mark on FoldTree); ref_used = first member (analysis ref)
        tm_labels = list(tm.iloc[:, 0].astype(str)) if tm is not None else members
        hub, hub_meanTM = _hub_from_tm(tm, tm_labels) if tm is not None else (members[0] if members else None, None)
        # FoldTree Newick outputs. The configured foldtree metric remains the interactive
        # tree; every available metric is retained in the family workbook.
        newick = ""
        trees = {}
        tree_status = {}
        for metric in config.get("signals", {}).get("foldtree_metrics", ["foldtree", "alntmscore", "lddt"]):
            tree_path = os.path.join(fd, f"{fam}_{metric}.nwk")
            if os.path.exists(tree_path):
                trees[str(metric)] = open(tree_path, encoding="utf-8", errors="replace").read().strip()
        tree_status_path = os.path.join(fd, f"{fam}_foldtree_status.json")
        if os.path.exists(tree_status_path):
            try:
                tree_status = json.load(open(tree_status_path, encoding="utf-8"))
            except Exception:
                tree_status = {}
        nwk_p = os.path.join(fd, f"{fam}_foldtree.nwk")
        if os.path.exists(nwk_p):
            newick = open(nwk_p).read().strip()
            assets["tree_svg"] = _svg_datauri(_newick_to_svg(newick, hub=hub))

        # ---- EXTRA: conservation + pocket (fpocket/P2Rank) + ESM + hub + cysteines ----
        primary_tree_status = tree_status.get("metrics", {}).get("foldtree", {})
        rooting_method = primary_tree_status.get("rooting_method")
        rooting_reason = primary_tree_status.get("reason")
        if rooting_method == "mad":
            rooting_label = "MAD root"
        elif rooting_method == "midpoint" and rooting_reason == "small_family_policy":
            rooting_label = "midpoint root (small-family policy)"
        elif rooting_method == "midpoint":
            rooting_label = "midpoint root (MAD fallback)"
        else:
            rooting_label = "rooting status unavailable"
        ex = dict(
            ref_used=members[0] if members else "",
            hub=hub,
            hub_meanTM=hub_meanTM,
            foldseek_tm_all_pairs=tm_stats,
            blast_identity_all_pairs=id_stats,
            mean_retained_edge_TM=retained_tm,
            mean_retained_edge_foldseek_fident=retained_fident,
            foldtree_status=tree_status,
            foldtree_rooting_label=rooting_label,
            sequence_analysis_status=sequence_status,
        )
        if structural_cons is not None and "structural_lddt" in structural_cons:
            values = structural_cons.structural_lddt.dropna()
            if len(values):
                ex["structural_lddt_min"] = float(values.min())
                ex["structural_lddt_max"] = float(values.max())
                ex["structural_lddt_mean"] = float(values.mean())
                ex["structural_scored_resi"] = [
                    int(value)
                    for value in structural_cons.loc[
                        structural_cons.structural_lddt.notna(), "resi"
                    ]
                ]
                ex["structural_pair_threshold"] = _finite_float(
                    structural_cons.get("pair_threshold", pd.Series([0.5])).iloc[0]
                )
        family_subgroups = (
            sequence_subgroups_all[
                sequence_subgroups_all.family.astype(str) == str(fam)
            ].copy()
            if sequence_subgroups_all is not None
            else pd.DataFrame()
        )
        family_domains = (
            domain_members_all[
                domain_members_all.acc.astype(str).isin(set(members))
            ].copy()
            if domain_members_all is not None
            else pd.DataFrame()
        )
        ex["n_sequence_subgroups"] = (
            int(family_subgroups.sequence_subgroup.nunique())
            if len(family_subgroups)
            else 0
        )
        ex["n_domain_families"] = (
            int(family_domains.domain_family.nunique())
            if len(family_domains)
            else 0
        )
        ex["domain_families"] = (
            sorted(family_domains.domain_family.dropna().astype(str).unique())
            if len(family_domains)
            else []
        )
        # US-align (independent-algorithm) TM cross-check summary for this family
        if tm_us_mean is not None:
            ex["tm_us_mean"] = tm_us_mean
            ex["tm_cons_r"] = tm_cons_r
            ex["tm_cons_maxdiff"] = tm_cons_maxdiff
            ex["tm_disagree"] = tm_disagree
        if sig is not None and "conservation" in sig:
            cons = sig["conservation"].dropna()
            sub = sig.dropna(subset=["rel_sasa", "conservation"]) if {"rel_sasa","conservation"}.issubset(sig.columns) else sig.iloc[0:0]
            if len(cons):
                ex["cons_min"] = float(cons.min())
                ex["cons_max"] = float(cons.max())
            ex["cons_sasa_r"] = float(np.corrcoef(sub.conservation, sub.rel_sasa)[0, 1]) if len(sub) > 2 else None
        # pockets: prefer P2Rank, keep both sources so the viewer can switch (add_p2rank_esmscan).
        # CRITICAL: the renderer's buildStructPane reads ex.fpocket_resi.length /
        # ex.p2rank_resi.length UNCONDITIONALLY — a missing key makes .length throw a
        # TypeError, aborting the whole struct-pane build so initViewer never runs and NO
        # structure ever renders. Initialise every pocket key the renderer touches with
        # safe empties for ALL families so those reads never throw.
        ex.setdefault("p2rank_resi", []); ex.setdefault("p2rank_prob", None)
        ex.setdefault("p2rank_score", None); ex.setdefault("p2rank_n", None)
        ex.setdefault("fpocket_resi", []); ex.setdefault("fpocket_score", None)
        ex.setdefault("pocket_resi", []); ex.setdefault("pocket_src", None)
        ex.setdefault("pocket_score", None); ex.setdefault("n_pocket", None)
        ex.setdefault("n_cys", 0)
        pk = _enrich_pocket_entry(results_dir, fam, pockets.get(fam, {}))
        pk.setdefault(
            "p2rank_profile",
            config.get("pocket", {}).get("p2rank_profile", "unknown"),
        )
        p2 = pk.get("p2rank", {}); fp = pk.get("fpocket", {})
        if p2:
            ex.update(p2rank_resi=p2.get("lining_residues", []), p2rank_score=p2.get("top_score"),
                      p2rank_n=p2.get("n_pockets"),
                      p2rank_prob=p2.get("top_probability"))
        if fp:
            ex.update(fpocket_resi=fp.get("lining_residues", []), fpocket_score=fp.get("top_score"))
        # default pocket shown = P2Rank if present else fpocket
        if p2:
            ex.update(pocket_src="p2rank", pocket_resi=p2.get("lining_residues", []),
                      pocket_score=p2.get("top_score"), n_pocket=p2.get("n_pockets"))
        elif fp:
            ex.update(pocket_src="fpocket", pocket_resi=fp.get("lining_residues", []),
                      pocket_score=fp.get("top_score"), n_pocket=fp.get("n_pockets"))
        # ESM per-site (mean LLR) for the reference; correlate with conservation & SASA
        if esm_all is not None and "family" in esm_all.columns:
            ef = esm_all[esm_all.family == fam]
            aa_cols = [c for c in ef.columns if len(str(c)) == 1 and str(c).isalpha()]
            if len(ef) and aa_cols:
                mean_llr = ef[aa_cols].mean(axis=1).values
                ex["has_esm"] = True
                ex["esm_min"] = float(np.nanmin(mean_llr)); ex["esm_max"] = float(np.nanmax(mean_llr))
                # record ESM ref + per-residue tolerance so REFPDB["<fam>_esm"] can be built
                poscol = next((c for c in ef.columns if str(c).lower() in ("", "unnamed: 0", "pos", "site")), ef.columns[0])
                pos = ef[poscol].astype(str).str.extract(r"(\d+)$")[0].astype(float)
                erf = str(ef["ref"].iloc[0]) if "ref" in ef.columns else None
                if erf:
                    esm_by_fam[fam] = {"ref": erf,
                                       "tol": {int(p): float(v) for p, v in zip(pos, mean_llr) if not np.isnan(p)}}
                if sig is not None and "conservation" in sig:
                    n = min(len(mean_llr), len(sig))
                    if n > 2:
                        cc = sig["conservation"].values[:n]
                        ss = sig["rel_sasa"].values[:n] if "rel_sasa" in sig else None
                        m2 = mean_llr[:n]
                        good = ~np.isnan(cc) & ~np.isnan(m2)
                        if good.sum() > 2:
                            ex["esm_vs_cons_r"] = float(np.corrcoef(cc[good], m2[good])[0, 1])
                        if ss is not None:
                            g2 = ~np.isnan(ss) & ~np.isnan(m2)
                            if g2.sum() > 2:
                                ex["esm_vs_sasa_r"] = float(np.corrcoef(ss[g2], m2[g2])[0, 1])
            else:
                ex["has_esm"] = False
        else:
            ex["has_esm"] = False
        # cysteine count on the reference conservation PDB (disulfide-rich effector signal)
        cons_pdb_p = os.path.join(fd, f"{fam}_conservation.pdb")
        if os.path.exists(cons_pdb_p):
            cys = set()
            for line in open(cons_pdb_p, encoding="utf-8", errors="replace"):
                if line.startswith("ATOM") and line[17:20].strip() == "CYS":
                    try: cys.add(int(line[22:26]))
                    except Exception: pass
            ex["n_cys"] = len(cys)
        EXTRA[fam] = ex
        # conservation-colored ref PDB
        cons_pdb = os.path.join(fd, f"{fam}_conservation.pdb")
        if os.path.exists(cons_pdb):
            store_reference(
                f"{fam}_cons",
                open(cons_pdb, encoding="utf-8", errors="replace").read(),
            )
        structural_pdb = os.path.join(fd, f"{fam}_structural_conservation.pdb")
        if os.path.exists(structural_pdb):
            store_reference(
                f"{fam}_struct",
                open(
                    structural_pdb, encoding="utf-8", errors="replace"
                ).read(),
            )
        # Backend mode keeps source structures for transforms/download generation but
        # does not duplicate them into the HTML payload.
        source_struct = {}
        for a in members:
            for cand in (os.path.join(fd, f"{a}.pdb"),
                         os.path.join(results_dir, "..", "input", "pdb", f"{config.get('strain',{}).get('code','')}_{a}.pdb")):
                if os.path.exists(cand):
                    source_struct[a] = open(
                        cand, encoding="utf-8", errors="replace"
                    ).read()
                    break
        struct = source_struct if mode == "single" else {}
        msa = _records_by_member(
            _read_fasta_records(os.path.join(fd, f"{fam}.aln")), members
        )
        three_di_msa = _records_by_member(
            _read_fasta_records(os.path.join(fd, f"{fam}.fasta")), members
        )
        sequence_msa = _records_by_member(
            _read_fasta_records(os.path.join(fd, f"{fam}_sequence_msa.fasta")),
            members,
        )
        sequence_tree_path = os.path.join(fd, f"{fam}_sequence_tree.nwk")
        sequence_tree = (
            open(sequence_tree_path, encoding="utf-8", errors="replace").read().strip()
            if os.path.exists(sequence_tree_path)
            else ""
        )
        if sequence_tree:
            trees["sequence"] = sequence_tree
            assets["sequence_tree_svg"] = _svg_datauri(
                _newick_to_svg(
                    sequence_tree,
                    hub=hub,
                    title="Sequence tree · MAFFT + FastTree",
                )
            )
        # FoldMason-aware rigid-body alignment to the canonical hub. Compact transforms
        # are embedded once and applied by the viewer and superposed-PDB downloader.
        transforms = {}
        fit_stats = {}
        if source_struct:
            ref_member = hub if hub in source_struct else next(iter(source_struct))
            ref_pdb = source_struct[ref_member]
            identity_rotation = np.eye(3).tolist()
            identity_translation = [0.0, 0.0, 0.0]
            transforms[ref_member] = {"rotation": identity_rotation, "translation": identity_translation}
            ref_n_ca = len(_ca_coordinates(ref_pdb))
            fit_stats[ref_member] = {"reference": ref_member, "method": "reference",
                                     "n_ca": ref_n_ca, "n_ca_total": ref_n_ca,
                                     "rmsd": 0.0, "rmsd_all": 0.0,
                                     "rotation": identity_rotation, "translation": identity_translation}
            for member, pdbtext in source_struct.items():
                if member == ref_member:
                    continue
                _, stats = _superpose_pdb(
                    pdbtext, ref_pdb, mobile_aln=msa.get(member), ref_aln=msa.get(ref_member))
                transforms[member] = {"rotation": stats["rotation"], "translation": stats["translation"]}
                fit_stats[member] = {"reference": ref_member, **stats}

        # Downloads are generated server-side so the self-contained HTML needs no ZIP or
        # spreadsheet runtime. Original structures remain one PDB per ZIP member.
        assets["structures_zip_b64"] = store_download(
            f"{fam}_member_structures.zip",
            _structures_zip_b64(fam, source_struct),
        )
        blast_pairs = None
        if classification_all is not None and {"q", "t"}.issubset(classification_all.columns):
            member_set = set(members)
            blast_pairs = classification_all[
                classification_all.q.astype(str).isin(member_set) &
                classification_all.t.astype(str).isin(member_set)
            ].copy()
        assets["xlsx_b64"] = store_download(
            f"{fam}_data.xlsx",
            _xlsx_b64(
            fam=fam, members=members,
            annotation=(anno[anno["family"].astype(str) == str(fam)].copy()
                        if len(anno) and "family" in anno.columns else None),
            tm=tm, usm=usm, idm=idm, blast_pairs=blast_pairs, sig=sig, exp=exp, pocket_entry=pk,
            pocket_raw=_pocket_raw_tables(results_dir, fam, pk), trees=trees,
            tree_status=tree_status, fit_stats=fit_stats,
            sequence_msa=sequence_msa, structural_msa=msa,
            three_di_msa=three_di_msa, sequence_status=sequence_status,
            domains=family_domains, subgroups=family_subgroups,
            structural_cons=structural_cons),
        )
        # ESM-tolerance-colored ref PDB: the renderer's "ESM" structure mode reads
        # REFPDB["<fam>_esm"]; without it, clicking the ESM button feeds addModel(undefined)
        # and blanks the viewer. Build it from the ESM ref's embedded structure + per-site
        # mean substitution score (tolerance) written into the B-factor column.
        if ex.get("has_esm") and fam in esm_by_fam:
            eref = esm_by_fam[fam]["ref"]
            eref_pdb = source_struct.get(eref) or (open(os.path.join(fd, f"{eref}.pdb"), encoding="utf-8", errors="replace").read()
                                            if os.path.exists(os.path.join(fd, f"{eref}.pdb")) else None)
            if eref_pdb:
                tol = esm_by_fam[fam]["tol"]
                out_l = []
                for line in eref_pdb.split("\n"):
                    if line.startswith(("ATOM", "HETATM")) and len(line) >= 66:
                        try: ri = int(line[22:26])
                        except ValueError: out_l.append(line); continue
                        out_l.append(f"{line[:60]}{tol.get(ri, 0.0):6.2f}{line[66:]}")
                    else:
                        out_l.append(line)
                store_reference(f"{fam}_esm", "\n".join(out_l))
        # per-member mature sequences (from seqs.fasta), fall back to CA-extraction from
        # the embedded structure so a member always has a downloadable sequence.
        seq = {}
        for a in members:
            s = seqs_all.get(a)
            if not s and a in source_struct:
                s = _seq_from_pdb(source_struct[a])
            if s:
                seq[a] = s
        PAY[fam] = dict(members=members, order=members, struct=struct, transforms=transforms,
                        seq=seq, msa=msa, structural_msa=msa,
                        three_di_msa=three_di_msa, sequence_msa=sequence_msa,
                        assets=assets, newick=newick, sequence_newick=sequence_tree,
                        domains=family_domains.to_dict("records") if len(family_domains) else [],
                        subgroups=family_subgroups.to_dict("records") if len(family_subgroups) else [],
                        maxid=float(r.get("max_identity", 0) or 0))
        if len(anno):
            payload = _annotation_payload(anno[anno.family == fam])
            if payload:
                ANN[fam] = payload

    # Singletons are independent proteins, not a synthetic cluster. They receive
    # structure, pocket, ESM, RNA-seq, and database-annotation payloads, while all
    # pairwise/family-only assets remain absent.
    SINGLETONS = []
    if members_all is not None and {"acc", "family"}.issubset(members_all.columns):
        singleton_members = members_all[members_all.family == "singleton"].copy()
        singleton_members = singleton_members.sort_values("acc")
        for _, member_row in singleton_members.iterrows():
            accession = str(member_row.acc)
            annotation_row = (
                anno[anno.acc.astype(str) == accession].copy()
                if len(anno) and "acc" in anno.columns
                else pd.DataFrame()
            )
            annotation_payload = _annotation_payload(annotation_row)
            if annotation_payload:
                ANN[accession] = annotation_payload
                member_annotation = annotation_payload["members"][0]
            else:
                member_annotation = {
                    "acc": accession, "gene": "", "novel": None, "tm": 0, "eff": "",
                    "pfam": "", "ipr": "", "pdb": "", "pdb_tm": None, "afdb": "",
                    "afdb_hit": "", "afdb_tm": None,
                }

            pdbtext = structure_text(accession)
            structures = {accession: pdbtext} if mode == "single" and pdbtext else {}
            sequence = seqs_all.get(accession) or (_seq_from_pdb(pdbtext) if pdbtext else "")
            sequences = {accession: sequence} if sequence else {}
            assets = {}

            expression = None
            expression_values = {}
            if expression_all is not None and "acc" in expression_all.columns:
                expression = expression_all[
                    expression_all.acc.astype(str) == accession
                ].copy()
                if len(expression):
                    numeric = expression.drop(columns=["acc"], errors="ignore").apply(
                        pd.to_numeric, errors="coerce"
                    )
                    for column in numeric.columns:
                        value = _finite_float(numeric.iloc[0][column])
                        if value is not None:
                            expression_values[str(column)] = value
                    if expression_values:
                        assets["rna_svg"] = _svg_datauri(
                            _svg_heat(
                                expression.set_index("acc"),
                                f"{accession} · RNA-seq expression",
                            )
                        )

            pocket_entry = _enrich_pocket_entry(
                results_dir, accession, pockets.get(accession, {})
            )
            pocket_entry.setdefault(
                "p2rank_profile",
                config.get("pocket", {}).get("p2rank_profile", "unknown"),
            )
            p2rank = pocket_entry.get("p2rank", {}) or {}
            fpocket = pocket_entry.get("fpocket", {}) or {}
            extra = {
                "kind": "singleton",
                "ref_used": accession,
                "p2rank_resi": p2rank.get("lining_residues", []),
                "p2rank_score": p2rank.get("top_score"),
                "p2rank_n": p2rank.get("n_pockets"),
                "p2rank_prob": p2rank.get("top_probability"),
                "fpocket_resi": fpocket.get("lining_residues", []),
                "fpocket_score": fpocket.get("top_score"),
                "pocket_resi": (
                    p2rank.get("lining_residues", [])
                    if p2rank
                    else fpocket.get("lining_residues", [])
                ),
                "pocket_src": "p2rank" if p2rank else "fpocket" if fpocket else None,
                "pocket_score": (
                    p2rank.get("top_score")
                    if p2rank
                    else fpocket.get("top_score")
                    if fpocket
                    else None
                ),
                "pocket_metric": (
                    "probability" if p2rank else "score" if fpocket else None
                ),
                "pocket_value": (
                    p2rank.get("top_probability")
                    if p2rank
                    else fpocket.get("top_score")
                    if fpocket
                    else None
                ),
                "n_pocket": (
                    p2rank.get("n_pockets")
                    if p2rank
                    else fpocket.get("n_pockets")
                    if fpocket
                    else None
                ),
                "has_esm": False,
                "n_cys": 0,
                "expression": expression_values,
            }
            if pdbtext:
                plddt_min, plddt_max = _bfactor_range(pdbtext)
                extra.update(plddt_min=plddt_min, plddt_max=plddt_max)
                cysteines = set()
                for line in pdbtext.splitlines():
                    if line.startswith("ATOM") and line[17:20].strip() == "CYS":
                        try:
                            cysteines.add(int(line[22:26]))
                        except ValueError:
                            pass
                extra["n_cys"] = len(cysteines)

            if esm_all is not None and "family" in esm_all.columns:
                esm_rows = esm_all[esm_all.family.astype(str) == accession]
                aa_columns = [
                    column for column in esm_rows.columns
                    if len(str(column)) == 1 and str(column).isalpha()
                ]
                if len(esm_rows) and aa_columns:
                    mean_llr = esm_rows[aa_columns].mean(axis=1).to_numpy()
                    position_column = next(
                        (
                            column for column in esm_rows.columns
                            if str(column).lower() in ("", "unnamed: 0", "pos", "site")
                        ),
                        esm_rows.columns[0],
                    )
                    positions = pd.to_numeric(
                        esm_rows[position_column].astype(str).str.extract(r"(\d+)$")[0],
                        errors="coerce",
                    )
                    tolerance = {
                        int(position): float(value)
                        for position, value in zip(positions, mean_llr)
                        if pd.notna(position) and math.isfinite(float(value))
                    }
                    extra.update(
                        has_esm=True,
                        esm_min=float(np.nanmin(mean_llr)),
                        esm_max=float(np.nanmax(mean_llr)),
                        esm_values=tolerance,
                    )

            peak_condition = None
            peak_expression = None
            if expression_values:
                peak_condition = max(expression_values, key=expression_values.get)
                peak_expression = expression_values[peak_condition]

            assets["xlsx_b64"] = store_download(
                f"{accession}_data.xlsx",
                _xlsx_b64(
                fam=accession,
                members=[accession],
                annotation=annotation_row if len(annotation_row) else None,
                tm=None,
                usm=None,
                idm=None,
                blast_pairs=None,
                sig=None,
                exp=expression if expression is not None and len(expression) else None,
                pocket_entry=pocket_entry,
                pocket_raw=_pocket_raw_tables(results_dir, accession, pocket_entry),
                trees={},
                tree_status={},
                fit_stats={},
                analysis_kind="singleton"),
            )
            PAY[accession] = {
                "kind": "singleton",
                "members": [accession],
                "order": [accession],
                "struct": structures,
                "transforms": {},
                "seq": sequences,
                "msa": {},
                "assets": assets,
                "newick": "",
                "maxid": 0,
            }
            EXTRA[accession] = extra

            plddt = _finite_float(member_row.get("plddt"))
            length = _finite_float(member_row.get("length"))
            SINGLETONS.append({
                "id": accession,
                "acc": accession,
                "gene": member_annotation.get("gene", ""),
                "label": (
                    annotation_payload.get("label", "novel/unknown")
                    if annotation_payload
                    else "annotation unavailable"
                ),
                "eff": member_annotation.get("eff", ""),
                "tmr": member_annotation.get("tm", 0),
                "novel": member_annotation.get("novel"),
                "pfam": member_annotation.get("pfam", ""),
                "ipr": member_annotation.get("ipr", ""),
                "pdb": member_annotation.get("pdb", ""),
                "pdb_tm": member_annotation.get("pdb_tm"),
                "afdb": member_annotation.get("afdb", ""),
                "afdb_hit": member_annotation.get("afdb_hit", ""),
                "afdb_tm": member_annotation.get("afdb_tm"),
                "plddt": plddt,
                "length": int(length) if length is not None else None,
                "pocket": bool(
                    extra.get("p2rank_resi") or extra.get("fpocket_resi")
                ),
                "pocket_method": extra.get("pocket_src"),
                "pocket_score": _finite_float(extra.get("pocket_score")),
                "pocket_metric": extra.get("pocket_metric"),
                "pocket_value": _finite_float(extra.get("pocket_value")),
                "has_esm": bool(extra.get("has_esm")),
                "rna_condition": peak_condition,
                "rna_peak": peak_expression,
                "rna": expression_values,
            })

    # cross-family structural edges: reuse classification/edges if a cross-fam TM file exists
    NET_edges = []
    xfam = os.path.join(results_dir, "cross_family_edges.csv")
    if os.path.exists(xfam):
        xe = pd.read_csv(xfam)
        for _, e in xe.iterrows():
            NET_edges.append(dict(**{"from": e["from"], "to": e["to"]},
                                  tm=float(e.tm), tm_max=float(e.get("tm_max", e.tm)), n=int(e.get("n", 1))))

    # ---- embed the two family-summary CSVs (clustered / singletons) for network-view downloads ----
    import base64 as _b64
    SUMMARY = {}
    for key, fn in [("clustered", "family_summary_clustered.csv"),
                    ("singletons", "family_summary_singletons.csv")]:
        fp = os.path.join(results_dir, fn)
        if os.path.exists(fp):
            with open(fp, "rb") as _fh:
                SUMMARY[key] = _b64.b64encode(_fh.read()).decode()

    member_parent = {}
    if members_all is not None and {"acc", "family"}.issubset(members_all.columns):
        member_parent = {
            str(row.acc): str(row.family)
            for row in members_all.itertuples()
        }
    annotation_by_acc = {}
    if len(anno) and "acc" in anno.columns:
        annotation_by_acc = {
            str(row.acc): row
            for _, row in anno.drop_duplicates("acc").iterrows()
        }

    enriched_domain_members = []
    if domain_members_all is not None:
        for _, segment in domain_members_all.iterrows():
            acc = str(segment.acc)
            parent = member_parent.get(acc, "")
            workspace = acc if parent == "singleton" else parent
            annotation_row = annotation_by_acc.get(acc)
            member_evidence = {}
            annotation_payload = ANN.get(workspace, {})
            for item in annotation_payload.get("members", []):
                if str(item.get("acc")) == acc:
                    member_evidence = item
                    break
            expression_values = {}
            if expression_all is not None and "acc" in expression_all.columns:
                expression_row = expression_all[
                    expression_all.acc.astype(str) == acc
                ]
                if len(expression_row):
                    for column, value in expression_row.iloc[0].items():
                        if column != "acc":
                            number = _finite_float(value)
                            if number is not None:
                                expression_values[str(column)] = number
            start, end = int(segment.start), int(segment.end)
            overlaps = _overlapping_annotations(annotation_row, start, end)
            pocket_source = pockets.get(acc, {})
            legacy_pocket = pockets.get(workspace, {})
            if not pocket_source and str(legacy_pocket.get("ref", "")) == acc:
                pocket_source = legacy_pocket
            protein_pocket = _enrich_pocket_entry(
                results_dir,
                acc,
                pocket_source,
            )
            p2rank = protein_pocket.get("p2rank", {}) or {}
            fpocket = protein_pocket.get("fpocket", {}) or {}
            for result in (p2rank, fpocket):
                result["domain_lining_residues"] = [
                    int(residue)
                    for residue in result.get("lining_residues", [])
                    if start <= int(residue) <= end
                ]
            esm_values = {
                residue: value
                for residue, value in _esm_tolerance(esm_all, acc).items()
                if start <= int(residue) <= end
            }
            parent_sequence = seqs_all.get(acc, "")
            enriched_domain_members.append({
                "domain_family": str(segment.domain_family),
                "segment_id": str(segment.segment_id),
                "acc": acc,
                "start": start,
                "end": end,
                "length": int(segment.length),
                "parent_length": len(parent_sequence) or None,
                "parent_family": parent,
                "workspace": workspace,
                "gene": member_evidence.get("gene", ""),
                "eff": member_evidence.get("eff", ""),
                "tmr": member_evidence.get("tm", 0),
                "novel": member_evidence.get("novel"),
                "pfam": member_evidence.get("pfam", ""),
                "ipr": member_evidence.get("ipr", ""),
                "pdb": member_evidence.get("pdb", ""),
                "pdb_tm": member_evidence.get("pdb_tm"),
                "afdb": member_evidence.get("afdb", ""),
                "afdb_hit": member_evidence.get("afdb_hit", ""),
                "afdb_tm": member_evidence.get("afdb_tm"),
                "overlap_annotations": overlaps,
                "pocket_residues": sorted(set(
                    p2rank.get("domain_lining_residues", [])
                    + fpocket.get("domain_lining_residues", [])
                )),
                "p2rank": p2rank,
                "fpocket": fpocket,
                "p2rank_profile": protein_pocket.get("p2rank_profile", ""),
                "pocket_status": protein_pocket.get("pocket_status", "not_run"),
                "esm_values": esm_values,
                "expression": expression_values,
                "structure_available": bool(
                    (PAY.get(workspace, {}).get("struct", {}) or {}).get(acc)
                ),
                "sequence_available": bool(
                    (PAY.get(workspace, {}).get("seq", {}) or {}).get(acc)
                ),
            })

    annotations_by_domain = {}
    parent_counts_by_domain = {}
    domain_links_by_workspace = {}
    for member in enriched_domain_members:
        labels = [
            annotation["label"]
            for annotation in member["overlap_annotations"]
            if annotation["label"]
        ]
        annotations_by_domain.setdefault(member["domain_family"], []).extend(labels)
        parent_label = (
            member["acc"]
            if member["parent_family"] == "singleton"
            else member["parent_family"]
        )
        parent_counts_by_domain.setdefault(
            member["domain_family"], Counter()
        )[parent_label] += 1
        domain_links_by_workspace.setdefault(member["workspace"], []).append({
            "domain_family": member["domain_family"],
            "segment_id": member["segment_id"],
            "acc": member["acc"],
            "start": member["start"],
            "end": member["end"],
        })
    for workspace, links in domain_links_by_workspace.items():
        if workspace in EXTRA:
            EXTRA[workspace]["domain_links"] = links
            EXTRA[workspace]["domain_families"] = sorted({
                link["domain_family"] for link in links
            })

    domain_family_records = []
    if domain_families_all is not None:
        for _, family in domain_families_all.iterrows():
            family_id = str(family.domain_family)
            labels = Counter(annotations_by_domain.get(family_id, []))
            top_annotation, top_annotation_count = (
                labels.most_common(1)[0] if labels else ("", 0)
            )
            domain_family_records.append({
                "domain_family": family_id,
                "n_segments": int(family.n_segments),
                "n_proteins": int(family.n_proteins),
                "n_edges": int(family.n_edges),
                "mean_probability": _finite_float(family.mean_probability),
                "mean_alntm": _finite_float(family.get("mean_alntm")),
                "mean_lddt": _finite_float(family.mean_lddt),
                "mean_query_coverage": _finite_float(
                    family.get("mean_query_coverage")
                ),
                "mean_target_coverage": _finite_float(
                    family.get("mean_target_coverage")
                ),
                "mean_aligned_residues": _finite_float(
                    family.mean_aligned_residues
                ),
                "top_annotation": top_annotation,
                "top_annotation_count": int(top_annotation_count),
                "n_annotated_segments": sum(
                    bool(member["overlap_annotations"])
                    for member in enriched_domain_members
                    if member["domain_family"] == family_id
                ),
                "parent_family_counts": [
                    {"family": parent, "n_segments": int(count)}
                    for parent, count in parent_counts_by_domain.get(
                        family_id, Counter()
                    ).most_common()
                ],
            })

    domain_edge_records = []
    if domain_edges_all is not None:
        for _, edge in domain_edges_all.iterrows():
            domain_edge_records.append({
                "domain_family": str(edge.domain_family),
                "source": str(edge.source),
                "target": str(edge.target),
                "evalue": _finite_float(edge.evalue),
                "prob": _finite_float(edge.prob),
                "bits": _finite_float(edge.bits),
                "lddt": _finite_float(edge.lddt),
                "alntmscore": _finite_float(getattr(edge, "alntmscore", None)),
                "fident": _finite_float(edge.fident),
                "qcov": _finite_float(getattr(edge, "qcov", None)),
                "tcov": _finite_float(getattr(edge, "tcov", None)),
                "alnlen": int(edge.alnlen),
                "shorter_coverage": _finite_float(edge.shorter_coverage),
            })

    domain_records_by_id = {
        record["domain_family"]: record for record in domain_family_records
    }
    for family_id, workbench in domain_workbench.get("families", {}).items():
        family_members = [
            member for member in enriched_domain_members
            if member["domain_family"] == family_id
        ]
        family_edges = [
            edge for edge in domain_edge_records
            if edge["domain_family"] == family_id
        ]
        parent_structures = {
            member["acc"]: structure_text(member["acc"])
            for member in family_members
        }
        parent_structures = {
            accession: text
            for accession, text in parent_structures.items()
            if text
        }
        parent_sequences = {
            member["acc"]: seqs_all.get(member["acc"], "")
            for member in family_members
            if seqs_all.get(member["acc"], "")
        }
        expression_rows = []
        seen_expression = set()
        for member in family_members:
            if member["acc"] in seen_expression or not member["expression"]:
                continue
            seen_expression.add(member["acc"])
            expression_rows.append(
                {"acc": member["acc"], **member["expression"]}
            )
        expression_table = (
            pd.DataFrame(expression_rows).set_index("acc")
            if expression_rows else None
        )
        workbench["rna_svg"] = (
            _svg_datauri(
                _svg_heat(
                    expression_table,
                    f"{family_id} parent-protein RNA-seq",
                )
            )
            if expression_table is not None else ""
        )
        workbook_b64 = _domain_xlsx_b64(
            family_id, family_members, family_edges, workbench
        )
        domain_structures_b64 = _domain_structures_zip_b64(
            family_id,
            family_members,
            parent_structures,
            segments_only=True,
        )
        parent_structures_b64 = _domain_structures_zip_b64(
            family_id,
            family_members,
            parent_structures,
            segments_only=False,
        )
        package_b64 = _domain_package_b64(
            family_id,
            family_members,
            family_edges,
            workbench,
            parent_structures,
            parent_sequences,
            workbook_b64,
        )
        workbench["assets"] = {
            "xlsx_b64": store_download(
                f"{family_id}_data.xlsx", workbook_b64
            ),
            "domain_structures_zip_b64": store_download(
                f"{family_id}_member_structures.zip",
                domain_structures_b64,
            ),
            "parent_structures_zip_b64": store_download(
                f"{family_id}_parent_structures.zip",
                parent_structures_b64,
            ),
            "package_zip_b64": store_download(
                f"{family_id}_domain_package.zip", package_b64
            ),
        }
        # Detector-native pocket records stay in the workbook/package. The browser
        # needs only the mapped top-pocket summary, so do not duplicate every raw
        # pocket and residue list in the HTML payload.
        for member in family_members:
            for method in ("p2rank", "fpocket"):
                result = member.get(method, {}) or {}
                member[method] = {
                    key: result.get(key)
                    for key in (
                        "top_score",
                        "top_probability",
                        "n_pockets",
                        "lining_residues",
                        "domain_lining_residues",
                    )
                    if result.get(key) is not None
                }
        record = domain_records_by_id.get(family_id)
        if record is not None:
            record["has_foldmason"] = bool(workbench.get("structural_msa"))
            record["has_foldtree"] = bool(workbench.get("foldtree_trees"))
            record["has_sequence_tree"] = bool(
                any(
                    subgroup.get("newick")
                    for subgroup in workbench.get("sequence_subgroups", [])
                )
            )
            record["has_usalign"] = bool(workbench.get("usalign_matrix"))

    domain_network_edges = []
    if domain_bridges_all is not None:
        for _, edge in domain_bridges_all.iterrows():
            domain_network_edges.append({
                "from": str(edge.source_family),
                "to": str(edge.target_family),
                "n": int(edge.n_edges),
                "prob": _finite_float(edge.mean_probability),
                "prob_max": _finite_float(edge.max_probability),
                "lddt": _finite_float(edge.mean_lddt),
                "lddt_max": _finite_float(edge.max_lddt),
                "alnlen": _finite_float(edge.mean_aligned_residues),
            })

    assigned_domain_accessions = {
        member["acc"] for member in enriched_domain_members
    }
    domain_unassigned_records = []
    if members_all is not None:
        for _, row in members_all.drop_duplicates("acc").iterrows():
            accession = str(row.acc)
            if accession in assigned_domain_accessions:
                continue
            annotation_row = annotation_by_acc.get(accession)
            annotation_payload = (
                _annotation_payload(pd.DataFrame([annotation_row]))
                if annotation_row is not None else None
            )
            expression_values = {}
            if expression_all is not None and "acc" in expression_all.columns:
                match = expression_all[
                    expression_all.acc.astype(str) == accession
                ]
                if len(match):
                    expression_values = {
                        str(column): number
                        for column, value in match.iloc[0].items()
                        if column != "acc"
                        and (number := _finite_float(value)) is not None
                    }
            pocket_entry = _enrich_pocket_entry(
                results_dir, accession, pockets.get(accession, {})
            )
            p2rank = pocket_entry.get("p2rank", {}) or {}
            fpocket = pocket_entry.get("fpocket", {}) or {}
            domain_unassigned_records.append({
                "acc": accession,
                "parent_family": str(row.family),
                "label": (
                    annotation_payload.get("label", "annotation unavailable")
                    if annotation_payload else "annotation unavailable"
                ),
                "plddt": _finite_float(row.get("plddt")),
                "length": int(row.length) if pd.notna(row.get("length")) else None,
                "p2rank_score": _finite_float(p2rank.get("top_score")),
                "p2rank_probability": _finite_float(
                    p2rank.get("top_probability")
                ),
                "fpocket_score": _finite_float(fpocket.get("top_score")),
                "expression": expression_values,
            })

    D = dict(
        ANALYSIS_SCOPE=str(
            config.get("output", {}).get("analysis_scope", "both")
        ),
        NET=dict(nodes=NET_nodes, edges=NET_edges),
        SINGLETONS=SINGLETONS,
        DOMAIN_UNASSIGNED=domain_unassigned_records,
        DOMAIN_FAMILIES=domain_family_records,
        DOMAIN_MEMBERS=enriched_domain_members,
        DOMAIN_EDGES=domain_edge_records,
        DOMAIN_WORKBENCH=domain_workbench,
        DNET=dict(nodes=domain_family_records, edges=domain_network_edges),
        EXTRA=EXTRA,
        REFPDB=REFPDB,
        REFAVAIL=REFAVAIL,
        PAY=PAY,
        SUMMARY=SUMMARY,
        BACKEND={"enabled": mode == "backend"},
    )

    # assemble via string splice (never re.sub)
    prefix = _read_tpl("prefix.html")            # ends with 'var D='

    # ---- page title: user-editable config.output.project_title, else auto from species ----
    # counts come from the actual data (not hardcoded) so a subset run reports its own size.
    n_fam = len(NET_nodes)
    n_prot = len(members_all) if members_all is not None else (
        int(pd.to_numeric(master.get("n_members"), errors="coerce").fillna(0).sum())
        if "n_members" in master.columns else 0
    )
    import html as _html
    user_title = str(config.get("output", {}).get("project_title", "") or "").strip()
    if user_title:
        head = _html.escape(user_title)
    else:
        species = str(config.get("strain", {}).get("species", "") or "").strip()
        head = (_html.escape(species) + " secretome") if species else "SUSS structural atlas"
        head += " &middot; SUSS structural atlas"
    counts = (
        f"{n_fam} families, {len(SINGLETONS)} singletons, {n_prot} secreted proteins"
        if n_prot else f"{n_fam} families, {len(SINGLETONS)} singletons"
    )
    title_html = f"{head} &mdash; {counts}"
    prefix = prefix.replace("__ATLAS_TITLE__", title_html)
    databridge = _read_tpl("databridge.js")
    renderer = _read_tpl("renderer.js")
    tail = _read_tpl("tail.html")
    # INLINE the JS libraries (vis-network + 3Dmol) so the atlas is fully self-contained:
    # the template references them by CDN <script src="https://cdnjs...">, which fails
    # offline / on lab networks that block cdnjs → blank network graph AND blank 3D viewer.
    # Replace each CDN script tag with the vendored library inlined verbatim.
    vend = os.path.join(_TPL, "vendor")
    for src_url, libfile in [
        ("https://cdnjs.cloudflare.com/ajax/libs/vis-network/9.1.6/dist/vis-network.min.js", "vis-network.min.js"),
        ("https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js", "3Dmol-min.js"),
    ]:
        libpath = os.path.join(vend, libfile)
        tag = f'<script src="{src_url}"></script>'
        if os.path.exists(libpath) and tag in prefix:
            js = open(libpath, encoding="utf-8").read()
            prefix = prefix.replace(tag, "<script>\n" + js + "\n</script>")
    doc = (prefix + json.dumps(D) + databridge + "var ANN=" + json.dumps(ANN)
           + renderer + tail)
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html, "w", encoding="utf-8") as fh:
        fh.write(doc)
    return dict(families=len(NET_nodes), singletons=len(SINGLETONS), edges=len(NET_edges),
                bytes=len(doc), mode=mode, out=out_html)
