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
import sys

_BUILDERS_DIR = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS_DIR = os.path.abspath(os.path.join(_BUILDERS_DIR, "..", "scripts"))
for _path in (_BUILDERS_DIR, _SCRIPTS_DIR):
    if _path not in sys.path:
        sys.path.insert(0, _path)

from id_utils import extract_accession
from pdb_io import (
    _aligned_ca_pairs,
    _bfactor_range,
    _ca_coordinates,
    _domain_structures_zip_b64,
    _pdb_with_bfactors,
    _seq_from_pdb,
    _structure_path,
    _structures_zip_b64,
    _superpose_pdb,
)
from svg_charts import (
    _newick_to_svg,
    _svg_datauri,
    _svg_heat,
    _svg_matrix,
)
from export_bundles import (
    _clean_cell,
    _domain_package_b64,
    _domain_xlsx_b64,
    _finite_float,
    _xlsx_b64,
)

_TPL = os.path.join(os.path.dirname(__file__), "template")


_AA3TO1 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
           'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
           'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}











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












def _hub_from_tm(tm, labels):
    """Select the best-covered member, then highest mean measured off-diagonal TM."""
    try:
        M = np.array(tm.set_index(tm.columns[0]).values, dtype=float)
        n = M.shape[0]
        if n < 2:
            return (labels[0] if labels else None), None, 0, 0
        np.fill_diagonal(M, np.nan)
        counts = np.isfinite(M).sum(axis=1)
        means = np.asarray([
            float(np.nanmean(row)) if np.isfinite(row).any() else np.nan
            for row in M
        ])
        candidates = sorted(
            range(n),
            key=lambda index: (
                -int(counts[index]),
                -means[index] if math.isfinite(means[index]) else math.inf,
                str(labels[index]),
            ),
        )
        i = candidates[0]
        mean = round(float(means[i]), 3) if math.isfinite(means[i]) else None
        return labels[i], mean, int(counts[i]), n - 1
    except Exception:
        return (labels[0] if labels else None), None, 0, max(0, len(labels) - 1)


def _matrix_pair_stats(table):
    """Return explicit off-diagonal stats for a labeled symmetric matrix."""
    empty = {
        "mean_all": None,
        "mean_detected": None,
        "maximum": None,
        "n_pairs": 0,
        "n_detected": 0,
        "n_possible": 0,
        "n_missing": 0,
    }
    if table is None or len(table) < 2:
        return empty
    try:
        matrix = table.set_index(table.columns[0]).apply(
            pd.to_numeric, errors="coerce"
        ).to_numpy(dtype=float)
        values = matrix[np.triu_indices(len(matrix), k=1)]
        values = values[np.isfinite(values)]
        possible = len(matrix) * (len(matrix) - 1) // 2
        if not len(values):
            return {**empty, "n_possible": possible, "n_missing": possible}
        detected = values[values > 0]
        return {
            "mean_all": round(float(values.mean()), 3),
            "mean_detected": (
                round(float(detected.mean()), 3) if len(detected) else None
            ),
            "maximum": round(float(values.max()), 3),
            "n_pairs": int(len(values)),
            "n_detected": int(len(detected)),
            "n_possible": possible,
            "n_missing": int(possible - len(values)),
        }
    except Exception:
        return empty


def _matrix_agreement_stats(foldseek, usalign, threshold=0.1):
    """Compare only pairs measured by both matrices; keep dense US-align mean."""
    result = {
        "usalign_mean": None,
        "pearson_r": None,
        "max_abs_diff": None,
        "n_disagree": 0,
        "n_compared": 0,
    }
    if foldseek is None or usalign is None:
        return result
    try:
        fs = foldseek.set_index(foldseek.columns[0]).apply(
            pd.to_numeric, errors="coerce"
        )
        us = usalign.set_index(usalign.columns[0]).apply(
            pd.to_numeric, errors="coerce"
        )
        us_values = us.to_numpy(dtype=float)
        us_upper = us_values[np.triu_indices(len(us_values), k=1)]
        us_upper = us_upper[np.isfinite(us_upper)]
        if len(us_upper):
            result["usalign_mean"] = round(float(us_upper.mean()), 3)

        shared = [
            label
            for label in fs.columns
            if label in fs.index and label in us.columns and label in us.index
        ]
        if len(shared) < 2:
            return result
        indices = np.triu_indices(len(shared), k=1)
        fs_values = fs.loc[shared, shared].to_numpy(dtype=float)[indices]
        us_values = us.loc[shared, shared].to_numpy(dtype=float)[indices]
        measured = np.isfinite(fs_values) & np.isfinite(us_values)
        fs_values = fs_values[measured]
        us_values = us_values[measured]
        result["n_compared"] = int(len(fs_values))
        if not len(fs_values):
            return result
        differences = np.abs(fs_values - us_values)
        result["max_abs_diff"] = round(float(differences.max()), 3)
        result["n_disagree"] = int((differences > threshold).sum())
        if (
            len(fs_values) >= 3
            and np.std(fs_values) > 0
            and np.std(us_values) > 0
        ):
            result["pearson_r"] = round(
                float(np.corrcoef(fs_values, us_values)[0, 1]), 3
            )
        return result
    except Exception:
        return result


def _site_correlations(signature, residue_values):
    """Correlate site scores by residue identifier, never by row position."""
    result = {"conservation": None, "sasa": None, "n_conservation": 0, "n_sasa": 0}
    if signature is None or "resi" not in signature or not residue_values:
        return result
    table = signature.copy()
    table["_site_value"] = pd.to_numeric(
        table["resi"], errors="coerce"
    ).map(residue_values)
    for source, output in (("conservation", "conservation"), ("rel_sasa", "sasa")):
        if source not in table:
            continue
        paired = table[[source, "_site_value"]].apply(
            pd.to_numeric, errors="coerce"
        ).dropna()
        result[f"n_{output}"] = int(len(paired))
        if (
            len(paired) >= 3
            and paired[source].std() > 0
            and paired["_site_value"].std() > 0
        ):
            result[output] = float(
                np.corrcoef(paired[source], paired["_site_value"])[0, 1]
            )
    return result


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
    domain_diagnostics_all = load_csv(
        os.path.join(results_dir, "domain_match_diagnostics.csv")
    )
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
        input_dir = str(config.get("input", {}).get("pdb_dir", "") or "")
        strain = str(config.get("strain", {}).get("code", "") or "")
        candidate = _structure_path(
            accession,
            input_dir=input_dir,
            strain=strain,
            family_dir=family_dir,
            fallback_dir=os.path.join(results_dir, "..", "input", "pdb"),
        )
        if candidate:
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
            custom_pat = str(config.get("input", {}).get("accession_pattern", "") or "").strip() or None
            strain_code = str(config.get("strain", {}).get("code", "") or "")
            for line in open(mem_file):
                raw_line = line.strip()
                if not raw_line:
                    continue
                acc_val = extract_accession(
                    os.path.basename(raw_line),
                    strain=strain_code,
                    custom_pattern=custom_pat,
                    strip_strain_prefix=True,
                )
                if acc_val:
                    members.append(acc_val)
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
        tm_us_mean = tm_cons_r = tm_cons_maxdiff = None
        tm_disagree = tm_cons_n_pairs = 0
        if usm is not None:
            ulabs = list(usm.iloc[:, 0].astype(str))
            UM = usm.set_index(usm.columns[0]).values.tolist()
            assets["tmus_svg"] = _svg_datauri(_svg_matrix(UM, [l[:11] for l in ulabs],
                                           f"{fam} · Structural similarity (US-align TM, independent)",
                                           vmin=0.3, vmax=1.0, unit="TM"))
            agreement = _matrix_agreement_stats(tm, usm)
            tm_us_mean = agreement["usalign_mean"]
            tm_cons_r = agreement["pearson_r"]
            tm_cons_maxdiff = agreement["max_abs_diff"]
            tm_disagree = agreement["n_disagree"]
            tm_cons_n_pairs = agreement["n_compared"]
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
        hub, hub_meanTM, hub_pairs_measured, hub_pairs_expected = (
            _hub_from_tm(tm, tm_labels)
            if tm is not None
            else (members[0] if members else None, None, 0, max(0, len(members) - 1))
        )
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
            hub_pairs_measured=hub_pairs_measured,
            hub_pairs_expected=hub_pairs_expected,
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
            ex["tm_cons_n_pairs"] = tm_cons_n_pairs
        if sig is not None and "conservation" in sig:
            cons = sig["conservation"].dropna()
            sub = sig.dropna(subset=["rel_sasa", "conservation"]) if {"rel_sasa","conservation"}.issubset(sig.columns) else sig.iloc[0:0]
            if len(cons):
                ex["cons_min"] = float(cons.min())
                ex["cons_max"] = float(cons.max())
                ex["sequence_scored_resi"] = [
                    int(value)
                    for value in sig.loc[sig.conservation.notna(), "resi"]
                ]
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
                tolerance = _esm_tolerance(esm_all, fam)
                mean_llr = np.asarray(list(tolerance.values()), dtype=float)
                ex["has_esm"] = bool(len(mean_llr))
                if len(mean_llr):
                    ex["esm_min"] = float(np.nanmin(mean_llr))
                    ex["esm_max"] = float(np.nanmax(mean_llr))
                    # Record the same residue-indexed values used for correlations.
                    erf = str(ef["ref"].iloc[0]) if "ref" in ef.columns else None
                    if erf:
                        esm_by_fam[fam] = {"ref": erf, "tol": tolerance}
                    site_stats = _site_correlations(sig, tolerance)
                    ex["esm_vs_cons_r"] = site_stats["conservation"]
                    ex["esm_vs_sasa_r"] = site_stats["sasa"]
                    ex["esm_vs_cons_n"] = site_stats["n_conservation"]
                    ex["esm_vs_sasa_n"] = site_stats["n_sasa"]
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
        if (
            os.path.exists(cons_pdb)
            and os.path.getsize(cons_pdb) > 0
            and ex.get("sequence_scored_resi")
            and ex.get("cons_min") is not None
            and ex.get("cons_max") is not None
        ):
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
        source_struct = {
            accession: pdbtext
            for accession in members
            if (pdbtext := structure_text(accession, family_dir=fd))
        }
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
    domain_diagnostics_by_acc = {}
    if domain_diagnostics_all is not None and "acc" in domain_diagnostics_all:
        domain_diagnostics_by_acc = {
            str(row.acc): row
            for row in domain_diagnostics_all.drop_duplicates("acc").itertuples()
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
            member_annotation = (
                annotation_payload.get("members", [{}])[0]
                if annotation_payload else {}
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
            peak_condition = None
            peak_expression = None
            if expression_values:
                peak_condition, peak_expression = max(
                    expression_values.items(), key=lambda item: item[1]
                )
            pocket_entry = _enrich_pocket_entry(
                results_dir, accession, pockets.get(accession, {})
            )
            p2rank = pocket_entry.get("p2rank", {}) or {}
            fpocket = pocket_entry.get("fpocket", {}) or {}
            diagnostic = domain_diagnostics_by_acc.get(accession)
            domain_unassigned_records.append({
                "acc": accession,
                "parent_family": str(row.family),
                "label": (
                    annotation_payload.get("label", "annotation unavailable")
                    if annotation_payload else "annotation unavailable"
                ),
                "gene": member_annotation.get("gene", ""),
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
                "annotation_status": member_annotation.get(
                    "annotation_status", ""
                ),
                "foldseek_pdb_status": member_annotation.get(
                    "foldseek_pdb_status", ""
                ),
                "foldseek_afdb_status": member_annotation.get(
                    "foldseek_afdb_status", ""
                ),
                "plddt": _finite_float(row.get("plddt")),
                "length": int(row.length) if pd.notna(row.get("length")) else None,
                "p2rank_resi": p2rank.get("lining_residues", []),
                "p2rank_score": _finite_float(p2rank.get("top_score")),
                "p2rank_probability": _finite_float(
                    p2rank.get("top_probability")
                ),
                "fpocket_resi": fpocket.get("lining_residues", []),
                "fpocket_score": _finite_float(fpocket.get("top_score")),
                "expression": expression_values,
                "rna_condition": peak_condition,
                "rna_peak": peak_expression,
                "domain_status": (
                    str(diagnostic.status)
                    if diagnostic is not None else "diagnostic_unavailable"
                ),
                "best_match": (
                    str(diagnostic.best_match)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "best_match", None))
                    else ""
                ),
                "protein_start": (
                    int(diagnostic.protein_start)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "protein_start", None))
                    else None
                ),
                "protein_end": (
                    int(diagnostic.protein_end)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "protein_end", None))
                    else None
                ),
                "match_start": (
                    int(diagnostic.match_start)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "match_start", None))
                    else None
                ),
                "match_end": (
                    int(diagnostic.match_end)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "match_end", None))
                    else None
                ),
                "evalue": _finite_float(
                    getattr(diagnostic, "evalue", None)
                    if diagnostic is not None else None
                ),
                "prob": _finite_float(
                    getattr(diagnostic, "prob", None)
                    if diagnostic is not None else None
                ),
                "alntmscore": _finite_float(
                    getattr(diagnostic, "alntmscore", None)
                    if diagnostic is not None else None
                ),
                "lddt": _finite_float(
                    getattr(diagnostic, "lddt", None)
                    if diagnostic is not None else None
                ),
                "fident": _finite_float(
                    getattr(diagnostic, "fident", None)
                    if diagnostic is not None else None
                ),
                "protein_coverage": _finite_float(
                    getattr(diagnostic, "protein_coverage", None)
                    if diagnostic is not None else None
                ),
                "match_coverage": _finite_float(
                    getattr(diagnostic, "match_coverage", None)
                    if diagnostic is not None else None
                ),
                "shorter_coverage": _finite_float(
                    getattr(diagnostic, "shorter_coverage", None)
                    if diagnostic is not None else None
                ),
                "alnlen": _finite_float(
                    getattr(diagnostic, "alnlen", None)
                    if diagnostic is not None else None
                ),
                "raw_hit_count": (
                    int(diagnostic.raw_hit_count)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "raw_hit_count", None))
                    else 0
                ),
                "failed_filters": (
                    str(diagnostic.failed_filters)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "failed_filters", None))
                    else ""
                ),
                "domain_summary": (
                    str(diagnostic.summary)
                    if diagnostic is not None
                    and pd.notna(getattr(diagnostic, "summary", None))
                    else "Domain diagnostic was not generated for this run."
                ),
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
