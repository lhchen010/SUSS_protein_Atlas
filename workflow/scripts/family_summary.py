"""family_summary rule — one row per cluster/family with everything merged into a single
downloadable Excel: members, consensus annotation (same logic as the atlas card),
pocket lining residues (fpocket + P2Rank), mean structural similarity (TM), mean/max
sequence identity, SUSS %, effector/novel %, conservation-buried r. Singletons appear as
their own one-member rows so nothing is dropped.
"""
import os, re, json
import numpy as np, pandas as pd

master_csv = snakemake.input.master
members_csv = snakemake.input.members
anno_csv = snakemake.input.annotation
pockets_json = snakemake.input.pockets
out_xlsx = snakemake.output.xlsx
out_clustered = snakemake.output.clustered
out_singletons = snakemake.output.singletons

mem = pd.read_csv(members_csv)                       # acc, family, community, deg, plddt, length
master = pd.read_csv(master_csv)                     # family-level metrics
anno = pd.read_csv(anno_csv) if os.path.exists(anno_csv) and os.path.getsize(anno_csv) > 0 else pd.DataFrame()
pockets = json.load(open(pockets_json)) if os.path.exists(pockets_json) else {}

# RNAseq: aggregate expression (acc x conditions), optional — present only when the run had RNAseq
expr_csv = getattr(snakemake.input, "expr", None)
expr = None; cond_cols = []
if expr_csv and os.path.exists(expr_csv) and os.path.getsize(expr_csv) > 0:
    try:
        _e = pd.read_csv(expr_csv)
        if len(_e) and _e.columns[0] == "acc" and _e.shape[1] > 1:
            expr = _e.set_index("acc")
            cond_cols = list(expr.columns)
    except Exception:
        expr = None

master_by_fam = master.set_index("family").to_dict("index")

def _top(g, col):
    """most-common token in a '|'-separated column, stripped of trailing '(...)' — matches card."""
    if col not in g: return "—", 0.0
    vals = [x for s in g[col].dropna() for x in re.split(r"\s*\|\s*", str(s)) if x and x != "nan"]
    vals = [re.sub(r"\s*\(.*", "", v).strip() for v in vals]
    if not vals: return "—", 0.0
    vc = pd.Series(vals).value_counts()
    return vc.index[0], round(100 * vc.iloc[0] / len(g), 1)

def pock_resis(fam, method):
    e = pockets.get(fam, {})
    p = e.get(method, {})
    r = p.get("lining_residues", []) if isinstance(p, dict) else []
    return (
        r,
        p.get("top_score") if isinstance(p, dict) else None,
        p.get("top_probability") if isinstance(p, dict) else None,
    )

def _usalign_stats(fam):
    """US-align mean and Foldseek agreement over mutually measured pairs."""
    fd = os.path.join(os.path.dirname(members_csv), "families", fam)
    up = os.path.join(fd, f"{fam}_TM_usalign.csv")
    fp = os.path.join(fd, f"{fam}_TM.csv")
    if not os.path.exists(up):
        return None, None, 0, 0
    try:
        us = pd.read_csv(up, index_col=0)
        iu = np.triu_indices(len(us), k=1)
        bv = us.values[iu].astype(float)
        bv = bv[np.isfinite(bv)]
        us_mean = round(float(bv.mean()), 3) if len(bv) else None
        r = None; disagree = compared = 0
        if os.path.exists(fp):
            fs = pd.read_csv(fp, index_col=0)
            shared = [c for c in fs.columns if c in us.columns and c in fs.index and c in us.index]
            if len(shared) >= 2:
                fa = fs.loc[shared, shared].values.astype(float)
                ua = us.loc[shared, shared].values.astype(float)
                j = np.triu_indices(len(shared), k=1)
                av, bb = fa[j], ua[j]
                measured = np.isfinite(av) & np.isfinite(bb)
                av, bb = av[measured], bb[measured]
                compared = int(len(av))
                if len(av) >= 3 and np.std(av) > 0 and np.std(bb) > 0:
                    r = round(float(np.corrcoef(av, bb)[0, 1]), 3)
                disagree = int((np.abs(av - bb) > 0.1).sum())
        return us_mean, r, disagree, compared
    except Exception:
        return None, None, 0, 0

def _matrix_stats(fam, suffix):
    """Explicit all-pair and detected-pair means for a family matrix."""
    path = os.path.join(
        os.path.dirname(members_csv), "families", fam, f"{fam}_{suffix}.csv"
    )
    if not os.path.exists(path):
        return None, None, None, 0, 0
    try:
        matrix = pd.read_csv(path, index_col=0).apply(
            pd.to_numeric, errors="coerce"
        ).values.astype(float)
        values = matrix[np.triu_indices(len(matrix), k=1)]
        values = values[np.isfinite(values)]
        detected = values[values > 0]
        return (
            round(float(values.mean()), 3) if len(values) else None,
            round(float(detected.mean()), 3) if len(detected) else None,
            round(float(values.max()), 3) if len(values) else None,
            int(len(values)),
            int(len(detected)),
        )
    except Exception:
        return None, None, None, 0, 0

rows = []
# iterate every family label present in members (includes 'singleton' bucket -> split per acc)
fam_groups = mem.groupby("family")
for fam, gmem in fam_groups:
    accs = list(gmem.acc)
    # singletons: one row each (they are their own trivial cluster)
    subgroups = [([a], a) for a in accs] if fam == "singleton" else [(accs, fam)]
    is_single = (fam == "singleton")
    for acc_list, label_fam in subgroups:
        current_mem = gmem[gmem.acc.isin(acc_list)]
        g_anno = anno[anno.acc.isin(acc_list)] if len(anno) else pd.DataFrame()
        m = master_by_fam.get(label_fam, {})
        top_name, _ = _top(g_anno, "afdbsp_name")
        top_pfam, top_pfam_frac = _top(g_anno, "pfam_domains")
        top_pdb, top_pdb_frac = _top(g_anno, "pdb_hit")
        consensus = (top_name if top_name != "—" else
                     top_pfam if top_pfam != "—" else
                     top_pdb if top_pdb != "—" else "novel/unknown")
        ref = pockets.get(label_fam, {}).get("ref", "")
        fp_res, fp_score, _ = pock_resis(label_fam, "fpocket")
        p2_res, p2_score, p2_probability = pock_resis(label_fam, "p2rank")
        n = len(acc_list)
        effector_complete = ("effectorp_status" not in g_anno or
                             (len(g_anno) and g_anno.effectorp_status.eq("complete").all()))
        pct_eff = (round(100 * g_anno.is_effector.mean(), 1)
                   if effector_complete and "is_effector" in g_anno and len(g_anno) else np.nan)
        known_novel = g_anno.novel.dropna() if "novel" in g_anno else pd.Series(dtype=float)
        pct_novel = round(100 * known_novel.mean(), 1) if len(known_novel) else np.nan
        # US-align independent TM cross-check (real families only; singletons have no pairs)
        _us_mean, _us_r, _us_disagree, _us_compared = (
            _usalign_stats(label_fam)
            if not is_single
            else (None, None, 0, 0)
        )
        _fs_all, _fs_detected, _fs_max, _fs_pairs, _fs_hits = (
            _matrix_stats(label_fam, "TM")
            if not is_single
            else (None, None, None, 0, 0)
        )
        _bl_all, _bl_detected, _bl_max, _bl_pairs, _bl_hits = (
            _matrix_stats(label_fam, "ID")
            if not is_single
            else (None, None, None, 0, 0)
        )
        # RNAseq: mean expression per condition across this family's members (blank if no RNAseq)
        rna = {}
        if expr is not None and cond_cols:
            present = [a for a in acc_list if a in expr.index]
            if present:
                mv = expr.loc[present, cond_cols].mean(axis=0)
                for cc in cond_cols:
                    rna[f"rnaseq_{cc}"] = round(float(mv[cc]), 2)
        rows.append(dict(
            family=label_fam,
            is_singleton=is_single,
            n_members=n,
            members="; ".join(sorted(acc_list)),
            consensus_annotation=consensus,
            top_pfam=top_pfam, top_pfam_pct=top_pfam_frac,
            top_pdb_fold=top_pdb, top_pdb_pct=top_pdb_frac,
            top_protein_name=top_name,
            foldseek_TM_all_pairs=_fs_all,
            foldseek_TM_detected_pairs=_fs_detected,
            n_structural_pairs=_fs_pairs,
            mean_retained_edge_TM=m.get(
                "mean_retained_edge_TM", m.get("mean_TM", np.nan)
            ),
            mean_TM_usalign=_us_mean,
            tm_foldseek_usalign_r=_us_r,
            tm_pairs_compared=_us_compared,
            tm_pairs_disagree=_us_disagree,
            blast_identity_all_pairs=_bl_all,
            blast_identity_detected_pairs=_bl_detected,
            max_blast_identity=_bl_max,
            n_blast_detected_pairs=_bl_hits,
            mean_retained_edge_foldseek_fident=m.get(
                "mean_retained_edge_foldseek_fident",
                m.get("mean_identity", np.nan),
            ),
            suss_pct=m.get("suss_pct", np.nan),
            cons_sasa_r=m.get("cons_sasa_r", np.nan),
            pct_effector=m.get("pct_effector", pct_eff),
            pct_novel=m.get("pct_novel", pct_novel),
            mean_pLDDT=m.get("mean_pLDDT", round(float(current_mem.plddt.mean()), 2) if "plddt" in current_mem else np.nan),
            mean_len=m.get("mean_len", int(current_mem.length.mean()) if "length" in current_mem else np.nan),
            pocket_ref=ref,
            fpocket_score=fp_score,
            fpocket_residues=" ".join(map(str, fp_res)),
            p2rank_score=p2_score,
            p2rank_probability=p2_probability,
            p2rank_residues=" ".join(map(str, p2_res)),
            **rna,
        ))

df = pd.DataFrame(rows)
# order: real clusters (F#) first by size desc, then singletons
df["_fnum"] = df.family.str.extract(r"^F(\d+)$").astype(float)
df = df.sort_values(["is_singleton", "_fnum"], na_position="last").drop(columns=["_fnum"])

os.makedirs(os.path.dirname(out_xlsx), exist_ok=True)
with pd.ExcelWriter(out_xlsx, engine="openpyxl") as xw:
    df.to_excel(xw, sheet_name="family_summary", index=False)
    # a short README sheet documenting each column
    doc = [
        ("family", "cluster ID (F0.. by size) or 'singleton'"),
        ("n_members", "number of proteins in the cluster"),
        ("members", "'; '-separated accession list"),
        ("consensus_annotation", "consensus label: real protein name > top Pfam > top PDB fold > novel/unknown"),
        ("top_pfam / top_pfam_pct", "most common Pfam domain and % of members carrying it"),
        ("top_pdb_fold / top_pdb_pct", "most common Foldseek PDB100 hit and % of members"),
        ("top_protein_name", "most common AFDB-SwissProt protein name (real name, not accession)"),
        ("foldseek_TM_all_pairs", "mean Foldseek TM across measured unique within-family pairs; missing hits remain NA"),
        ("foldseek_TM_detected_pairs", "mean Foldseek TM among pairs with a reported Foldseek hit"),
        ("mean_retained_edge_TM", "mean Foldseek TM among graph edges retained by the clustering thresholds"),
        ("mean_TM_usalign", "mean pairwise US-align TM within the cluster — algorithm-independent cross-check of Foldseek TM"),
        ("tm_foldseek_usalign_r", "Pearson r between Foldseek and US-align TM over mutually measured pairs (high = robust clustering)"),
        ("tm_pairs_compared", "number of mutually measured Foldseek/US-align pairs used for agreement statistics"),
        ("tm_pairs_disagree", "number of within-cluster pairs where |Foldseek TM − US-align TM| > 0.1"),
        ("blast_identity_all_pairs", "mean best-HSP BLASTp identity across all unique protein pairs; undetected pairs are zero"),
        ("blast_identity_detected_pairs / max_blast_identity", "mean / max best-HSP BLASTp identity among detected pairs"),
        ("mean_retained_edge_foldseek_fident", "mean Foldseek alignment identity among retained structural graph edges; this is not BLAST identity"),
        ("suss_pct", "% of within-family structural links that are core_SUSS (BLAST-undetectable)"),
        ("cons_sasa_r", "Pearson r(conservation, relative SASA) on the reference; negative expected"),
        ("pct_effector / pct_novel", "% of members called effector (EffectorP) / novel (no fold AND no domain)"),
        ("mean_pLDDT / mean_len", "mean AF2 pLDDT / mean sequence length"),
        ("pocket_ref", "reference structure (family hub) pockets were detected on"),
        ("fpocket_score / fpocket_residues", "fpocket top-pocket score and space-separated lining residue numbers"),
        ("p2rank_score / p2rank_probability / p2rank_residues",
         "P2Rank top-pocket ranking score, calibrated probability, and lining residue numbers"),
    ]
    if cond_cols:
        doc.append(("rnaseq_<condition>", "mean expression across the cluster's members for each RNAseq "
                    "condition (" + ", ".join(cond_cols) + "); replicate-averaged counts. Blank if no RNAseq."))
    pd.DataFrame(doc, columns=["column", "meaning"]).to_excel(xw, sheet_name="README", index=False)
    # split sheets: clustered families vs singletons (mirrors the two network-view downloads)
    clustered = df[~df.is_singleton]
    singles = df[df.is_singleton]
    clustered.to_excel(xw, sheet_name="clustered", index=False)
    singles.to_excel(xw, sheet_name="singletons", index=False)

# Declared split outputs let Snakemake track the atlas download dependencies.
df[~df.is_singleton].to_csv(out_clustered, index=False)
df[df.is_singleton].to_csv(out_singletons, index=False)

print(f"family_summary: {len(df)} rows ({(~df.is_singleton).sum()} clustered + "
      f"{df.is_singleton.sum()} singletons), rnaseq_cols={len(cond_cols)} -> {out_xlsx}")
