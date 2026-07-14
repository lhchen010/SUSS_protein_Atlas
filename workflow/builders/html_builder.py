"""html_builder.build_atlas — interactive atlas HTML builder (ported from v19).

Ported from the validated co_suss_network.html v19 (artifact c1b281c8). The renderer
(vis.js family network + draggable/hideable panels, embedded 3Dmol.js viewer with
cartoon/surface/stick/sphere/line × B-factor coloring for conservation/ESM/pocket,
FoldTree with hub gold-star, TM/ID/RNAseq matrices, Blob-based structure & Excel
downloads) is stored verbatim as four constant template halves under template/:
    prefix.html    head+body+CSS+vis setup, ends with 'var D='
    databridge.js  ';var NET=D.NET,EXTRA=D.EXTRA,REFPDB=D.REFPDB,PAY=D.PAY;'
    renderer.js    the ~21KB of pure renderer functions (constant)
    tail.html      '</script>' + closing tags
build_atlas regenerates the DATA objects (D = {NET, EXTRA, REFPDB, PAY} and ANN) from
the engine's rule outputs, then assembles:  prefix + json(D) + databridge +
'var ANN=' + json(ANN) + renderer + tail.

Payload schema (v19 target on the left; what build_atlas populates noted inline).
Fields are populated ONLY when the producing rule ran (step toggled on) and its output
file exists; a family missing a given output simply omits those keys.
  NET.nodes[] = {id, n, tm, id_pct, suss, plddt, len, maxid}         # from master.csv
  NET.edges[] = {from, to, tm, tm_max, n}   # cross-family TM; only if cross_family_edges.csv exists
  PAY[fam]    = {members[], order[], struct{acc:pdbtext}, assets, newick, maxid}
                struct: populated in single mode (embedded PDBs); empty in backend mode.
                assets: tm_svg/id_svg (matrix rules), rna_svg (rnaseq), tree_svg (foldtree),
                        xlsx_b64 (per-family workbook TM/seq-id/per-site/RNAseq — always built).
  EXTRA[fam]  = cons_min/cons_max/cons_sasa_r (conservation rule) · pocket_src/pocket_resi/
                pocket_score/n_pocket + p2rank_resi/p2rank_score/p2rank_n/p2rank_prob +
                fpocket_resi/fpocket_score (pocket rule, pockets.json; P2Rank preferred,
                both kept so the viewer can switch) · has_esm/esm_min/esm_max/esm_vs_cons_r/
                esm_vs_sasa_r (esm rule, esm_all.csv) · hub/hub_meanTM (highest mean-TM
                member, gold-starred on the FoldTree) · n_cys (CYS on ref) · ref_used.
                A step toggled off (e.g. pocket:false) → those keys absent; the viewer's
                gating (var pock=EXTRA[curFam]?EXTRA[curFam].pocket_resi:[]) handles it.
  REFPDB["<fam>_cons"] = conservation-B-factor ref PDB text (conservation rule)
  ANN[fam]    = {label, n, pct_novel, pct_eff, members[]}   # from member_annotation.csv

Critical implementation notes (kept — hard-won this session):
- String-splice assembly, NOT re.sub (backslashes un-escape JSON control chars).
- Downloads via Blob+URL.createObjectURL (<a download href="data:"> blocked in artifact iframe).
- Superposed structures embed backbone N/CA/C/O (not CA-only) so 3Dmol draws cartoon/stick/line.
- pocket gating generalized: var pock=EXTRA[curFam]?EXTRA[curFam].pocket_resi:[]
- single mode: structures embedded (foldcomp-compressible); backend mode: metadata only,
  structures lazy-load from the 4070 portal.
"""
import os, io, json, glob, base64, math, re, zipfile
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


def _enrich_pocket_entry(results_dir, fam, entry):
    """Backfill all pocket predictions from raw outputs for pre-v1.0.2 runs."""
    entry = json.loads(json.dumps(entry or {}))
    ref = entry.get("ref", "")

    p2 = entry.get("p2rank", {}) or {}
    if not p2.get("pockets"):
        candidates = glob.glob(os.path.join(results_dir, "p2rank", fam, "out", "*_predictions.csv"))
        if candidates:
            table = pd.read_csv(candidates[0])
            table.columns = [str(c).strip() for c in table.columns]
            predictions = []
            for idx, pred in table.iterrows():
                tokens = str(pred.get("residue_ids", "")).split()
                residues = sorted({int(x.split("_")[-1]) for x in tokens
                                   if x.split("_")[-1].isdigit()})
                predictions.append({
                    "pocket_id": int(pred.get("rank", idx + 1)),
                    "score": float(pred.get("score", 0)),
                    "lining_residues": residues,
                })
            if predictions:
                top = max(predictions, key=lambda p: p["score"])
                p2.update(top_score=top["score"], n_pockets=len(predictions),
                          lining_residues=top["lining_residues"], pockets=predictions)
                entry["p2rank"] = p2

    fp = entry.get("fpocket", {}) or {}
    if ref and not fp.get("pockets"):
        root = os.path.join(results_dir, "fpocket", fam, f"{ref}_out")
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
    p2_candidates = glob.glob(os.path.join(results_dir, "p2rank", fam, "out", "*_predictions.csv"))
    if p2_candidates:
        p2_table = pd.read_csv(p2_candidates[0])
        p2_table.columns = [str(c).strip() for c in p2_table.columns]
        tables["p2rank_pockets"] = p2_table

    ref = (entry or {}).get("ref", "")
    info = os.path.join(results_dir, "fpocket", fam, f"{ref}_out", f"{ref}_info.txt") if ref else ""
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


def _newick_to_svg(nwk, hub=None):
    """Minimal rectangular cladogram from a Newick string; gold-star the hub leaf."""
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
         '<text x="6" y="14" font-size="12" font-weight="600">FoldTree</text>']
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


def _xlsx_b64(fam, members, tm, usm, idm, blast_pairs, sig, exp, pocket_entry,
              pocket_raw, trees, fit_stats):
    """Build the complete, auditable per-family analysis workbook."""
    try:
        buf = io.BytesIO()
        with pd.ExcelWriter(buf, engine="openpyxl") as xl:
            pd.DataFrame({"family": [fam] * len(members), "member": members}).to_excel(
                xl, sheet_name="members", index=False)
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
                    n_pockets=result.get("n_pockets"), top_score=result.get("top_score"),
                    top_lining_residues=" ".join(map(str, result.get("lining_residues", [])))))
                for pred in result.get("pockets", []):
                    predictions.append(dict(
                        family=fam, reference=pocket_entry.get("ref", ""), method=method,
                        pocket_id=pred.get("pocket_id"), score=pred.get("score"),
                        n_residues=len(pred.get("lining_residues", [])),
                        lining_residues=" ".join(map(str, pred.get("lining_residues", [])))))
            pd.DataFrame(summaries).to_excel(xl, sheet_name="pocket_summary", index=False)
            pd.DataFrame(predictions, columns=["family", "reference", "method", "pocket_id", "score",
                                                      "n_residues", "lining_residues"]).to_excel(
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

            pd.DataFrame([{"metric": metric, "newick": newick} for metric, newick in trees.items()]).to_excel(
                xl, sheet_name="foldtree", index=False)
            fit_rows = []
            for member, stats in fit_stats.items():
                row = {"member": member, **stats}
                for key in ("rotation", "translation"):
                    if key in row:
                        row[key] = json.dumps(row[key], separators=(",", ":"))
                fit_rows.append(row)
            pd.DataFrame(fit_rows).to_excel(
                xl, sheet_name="superposition", index=False)
            pd.DataFrame([
                ("foldseek_TM", "Within-family symmetric Foldseek TM-score matrix used for structural clustering."),
                ("usalign_TM", "Independent within-family US-align TM-score matrix."),
                ("blast_identity", "Best-HSP BLASTp identity matrix."),
                ("blast_pairs", "Pair-level BLAST and SUSS classification data for structural edges."),
                ("pocket_summary", "Detector status and top-pocket summary for fpocket and P2Rank."),
                ("pocket_predictions", "Every pocket reported by each detector, including scores and lining residues."),
                ("pocket_residues", "One row per detector, pocket, and lining residue."),
                ("fpocket_pockets", "All descriptors parsed from the detector-native fpocket info file."),
                ("p2rank_pockets", "Complete detector-native P2Rank predictions table with all original columns."),
                ("foldtree", "All available FoldTree Newick trees, one row per configured metric."),
                ("RNAseq", "Per-member, replicate-collapsed RNA-seq expression for this family."),
                ("per_site", "Reference-residue conservation, SASA, pocket, and other site-level evidence."),
                ("superposition", "Hub-referenced rigid-body fit method, paired CA count, and RMSD."),
            ], columns=["sheet", "contents"]).to_excel(xl, sheet_name="README", index=False)
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


def build_atlas(master_csv, cards_dir, composition_xlsx, annotation_csv,
                results_dir, out_html, mode="single", atlas_name="atlas", config=None):
    config = config or {}
    famdir = os.path.join(results_dir, "families")
    master = pd.read_csv(master_csv)
    anno = pd.read_csv(annotation_csv) if os.path.exists(annotation_csv) else pd.DataFrame()

    def load_csv(p):
        return pd.read_csv(p) if os.path.exists(p) else None

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

    NET_nodes, PAY, EXTRA, REFPDB, ANN = [], {}, {}, {}, {}
    esm_by_fam = {}   # fam -> {ref, tol{resi:score}} for building REFPDB["<fam>_esm"]
    fam_accs = {}
    for _, r in master.iterrows():
        fam = r.family
        fd = os.path.join(famdir, fam)
        mem_file = os.path.join(results_dir, "families", f"{fam}.members.txt")
        members = []
        if os.path.exists(mem_file):
            import re
            for line in open(mem_file):
                m = re.search(r"[A-Z]{2,3}\d{4,}\.\d+", line)
                if m: members.append(m.group(0))
        fam_accs[fam] = members
        NET_nodes.append(dict(id=fam, n=int(r.n_members),
                              tm=float(r.get("mean_TM", 0) or 0), id_pct=float(r.get("mean_identity", 0) or 0),
                              suss=float(r.get("suss_pct", 0) or 0), plddt=float(r.get("mean_pLDDT", 0) or 0),
                              len=float(r.get("mean_len", 0) or 0), maxid=float(r.get("max_identity", 0) or 0)))
        # assets
        tm = load_csv(os.path.join(fd, f"{fam}_TM.csv"))
        idm = load_csv(os.path.join(fd, f"{fam}_ID.csv"))
        exp = load_csv(os.path.join(fd, f"{fam}_expression.csv"))
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
        # hub = highest mean-TM member (mark on FoldTree); ref_used = first member (analysis ref)
        tm_labels = list(tm.iloc[:, 0].astype(str)) if tm is not None else members
        hub, hub_meanTM = _hub_from_tm(tm, tm_labels) if tm is not None else (members[0] if members else None, None)
        # FoldTree Newick outputs. The configured foldtree metric remains the interactive
        # tree; every available metric is retained in the family workbook.
        newick = ""
        trees = {}
        for metric in config.get("signals", {}).get("foldtree_metrics", ["foldtree", "alntmscore", "lddt"]):
            tree_path = os.path.join(fd, f"{fam}_{metric}.nwk")
            if os.path.exists(tree_path):
                trees[str(metric)] = open(tree_path, encoding="utf-8", errors="replace").read().strip()
        nwk_p = os.path.join(fd, f"{fam}_foldtree.nwk")
        if os.path.exists(nwk_p):
            newick = open(nwk_p).read().strip()
            assets["tree_svg"] = _svg_datauri(_newick_to_svg(newick, hub=hub))

        # ---- EXTRA: conservation + pocket (fpocket/P2Rank) + ESM + hub + cysteines ----
        ex = dict(ref_used=members[0] if members else "", hub=hub, hub_meanTM=hub_meanTM)
        # US-align (independent-algorithm) TM cross-check summary for this family
        if tm_us_mean is not None:
            ex["tm_us_mean"] = tm_us_mean
            ex["tm_cons_r"] = tm_cons_r
            ex["tm_cons_maxdiff"] = tm_cons_maxdiff
            ex["tm_disagree"] = tm_disagree
        if sig is not None and "conservation" in sig:
            cons = sig["conservation"].dropna()
            sub = sig.dropna(subset=["rel_sasa", "conservation"]) if {"rel_sasa","conservation"}.issubset(sig.columns) else sig.iloc[0:0]
            ex["cons_min"] = float(cons.min()); ex["cons_max"] = float(cons.max())
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
        p2 = pk.get("p2rank", {}); fp = pk.get("fpocket", {})
        if p2:
            ex.update(p2rank_resi=p2.get("lining_residues", []), p2rank_score=p2.get("top_score"),
                      p2rank_n=p2.get("n_pockets"), p2rank_prob=p2.get("top_score"))
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
            REFPDB[f"{fam}_cons"] = open(cons_pdb, encoding="utf-8", errors="replace").read()
        # structures (single mode embeds; backend mode omits)
        struct = {}
        if mode == "single":
            for a in members:
                for cand in (os.path.join(fd, f"{a}.pdb"),
                             os.path.join(results_dir, "..", "input", "pdb", f"{config.get('strain',{}).get('code','')}_{a}.pdb")):
                    if os.path.exists(cand):
                        struct[a] = open(cand, encoding="utf-8", errors="replace").read(); break
        # FoldMason-aware rigid-body alignment to the canonical hub. Compact transforms
        # are embedded once and applied by the viewer and superposed-PDB downloader.
        transforms = {}
        fit_stats = {}
        if struct:
            ref_member = hub if hub in struct else next(iter(struct))
            msa = _records_by_member(_read_fasta_records(os.path.join(fd, f"{fam}.aln")), members)
            ref_pdb = struct[ref_member]
            identity_rotation = np.eye(3).tolist()
            identity_translation = [0.0, 0.0, 0.0]
            transforms[ref_member] = {"rotation": identity_rotation, "translation": identity_translation}
            ref_n_ca = len(_ca_coordinates(ref_pdb))
            fit_stats[ref_member] = {"reference": ref_member, "method": "reference",
                                     "n_ca": ref_n_ca, "n_ca_total": ref_n_ca,
                                     "rmsd": 0.0, "rmsd_all": 0.0,
                                     "rotation": identity_rotation, "translation": identity_translation}
            for member, pdbtext in struct.items():
                if member == ref_member:
                    continue
                _, stats = _superpose_pdb(
                    pdbtext, ref_pdb, mobile_aln=msa.get(member), ref_aln=msa.get(ref_member))
                transforms[member] = {"rotation": stats["rotation"], "translation": stats["translation"]}
                fit_stats[member] = {"reference": ref_member, **stats}

        # Downloads are generated server-side so the self-contained HTML needs no ZIP or
        # spreadsheet runtime. Original structures remain one PDB per ZIP member.
        assets["structures_zip_b64"] = _structures_zip_b64(fam, struct)
        blast_pairs = None
        if classification_all is not None and {"q", "t"}.issubset(classification_all.columns):
            member_set = set(members)
            blast_pairs = classification_all[
                classification_all.q.astype(str).isin(member_set) &
                classification_all.t.astype(str).isin(member_set)
            ].copy()
        assets["xlsx_b64"] = _xlsx_b64(
            fam=fam, members=members, tm=tm, usm=usm, idm=idm, blast_pairs=blast_pairs,
            sig=sig, exp=exp, pocket_entry=pk,
            pocket_raw=_pocket_raw_tables(results_dir, fam, pk), trees=trees, fit_stats=fit_stats)
        # ESM-tolerance-colored ref PDB: the renderer's "ESM" structure mode reads
        # REFPDB["<fam>_esm"]; without it, clicking the ESM button feeds addModel(undefined)
        # and blanks the viewer. Build it from the ESM ref's embedded structure + per-site
        # mean substitution score (tolerance) written into the B-factor column.
        if ex.get("has_esm") and fam in esm_by_fam:
            eref = esm_by_fam[fam]["ref"]
            eref_pdb = struct.get(eref) or (open(os.path.join(fd, f"{eref}.pdb"), encoding="utf-8", errors="replace").read()
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
                REFPDB[f"{fam}_esm"] = "\n".join(out_l)
        # per-member mature sequences (from seqs.fasta), fall back to CA-extraction from
        # the embedded structure so a member always has a downloadable sequence.
        seq = {}
        for a in members:
            s = seqs_all.get(a)
            if not s and a in struct:
                s = _seq_from_pdb(struct[a])
            if s:
                seq[a] = s
        PAY[fam] = dict(members=members, order=members, struct=struct, transforms=transforms,
                        seq=seq, assets=assets,
                        newick=newick, maxid=float(r.get("max_identity", 0) or 0))
        # ANN — populate every key the v19 renderer reads (label, pct_domain, pct_eff,
        # pct_novel, top_pfam/top_pfam_frac, top_pdb/top_pdb_frac, n_multi, fusion, members)
        if len(anno):
            g = anno[anno.family == fam]
            if len(g):
                n = len(g)
                def _top(col):
                    if col not in g: return "—", 0.0
                    vals = [str(x) for s in g[col].dropna() for x in re.split(r"\s*\|\s*", str(s)) if x and x != "nan"]
                    vals = [re.sub(r"\s*\(.*", "", v).strip() for v in vals]
                    if not vals: return "—", 0.0
                    vc = pd.Series(vals).value_counts()
                    return vc.index[0], round(100 * vc.iloc[0] / n, 1)
                top_pfam, top_pfam_frac = _top("pfam_domains")
                top_pdb, top_pdb_frac = _top("pdb_hit")
                top_name, top_name_frac = _top("afdbsp_name")   # real protein NAME (not accession)
                pct_domain = round(100 * (g.has_any_domain.mean() if "has_any_domain" in g else 0), 1)
                n_multi = int(g.multi_domain.sum()) if "multi_domain" in g else 0
                # consensus label: real protein name > top Pfam > top PDB fold > novel
                label = (top_name if top_name != "—" else
                         top_pfam if top_pfam != "—" else
                         top_pdb if top_pdb != "—" else "novel/unknown")
                # per-member rows: the renderer's annHTML reads member OBJECTS with fields
                # acc/novel/tm/eff/pfam/pdb/afdb — a bare accession string makes every cell
                # undefined. Build the objects here.
                def _cell(v):
                    s = "" if v is None else str(v)
                    return "" if s.lower() == "nan" else s
                mem_objs = []
                for _, mr in g.iterrows():
                    mem_objs.append(dict(
                        acc=str(mr.acc),
                        novel=bool(mr.novel) if "novel" in g else False,
                        tm=int(mr.n_TMR) if "n_TMR" in g and pd.notna(mr.n_TMR) else 0,
                        eff=_cell(mr.effectorp) if "effectorp" in g else "",
                        pfam=_cell(mr.pfam_domains) if "pfam_domains" in g else "",
                        pdb=_cell(mr.pdb_hit) if "pdb_hit" in g else "",
                        afdb=_cell(mr.afdbsp_name) if "afdbsp_name" in g else
                             (_cell(mr.afdbsp_hit) if "afdbsp_hit" in g else "")))
                ANN[fam] = dict(label=label, n=int(n),
                                pct_novel=round(100 * g.novel.mean(), 1) if "novel" in g else 0,
                                pct_eff=round(100 * g.is_effector.mean(), 1) if "is_effector" in g else 0,
                                pct_domain=pct_domain,
                                top_pfam=top_pfam, top_pfam_frac=top_pfam_frac,
                                top_pdb=top_pdb, top_pdb_frac=top_pdb_frac,
                                top_ipr=_top("interpro_entries")[0],
                                n_multi=n_multi, fusion=(n_multi > 0),
                                members=mem_objs)

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

    D = dict(NET=dict(nodes=NET_nodes, edges=NET_edges), EXTRA=EXTRA, REFPDB=REFPDB, PAY=PAY, SUMMARY=SUMMARY)

    # assemble via string splice (never re.sub)
    prefix = _read_tpl("prefix.html")            # ends with 'var D='

    # ---- page title: user-editable config.output.project_title, else auto from species ----
    # counts come from the actual data (not hardcoded) so a subset run reports its own size.
    n_fam = len(NET_nodes)
    n_prot = int(pd.to_numeric(master.get("n_members"), errors="coerce").fillna(0).sum()) \
        if "n_members" in master.columns else 0
    import html as _html
    user_title = str(config.get("output", {}).get("project_title", "") or "").strip()
    if user_title:
        head = _html.escape(user_title)
    else:
        species = str(config.get("strain", {}).get("species", "") or "").strip()
        head = (_html.escape(species) + " secretome") if species else "SUSS structural atlas"
        head += " &middot; SUSS structural atlas"
    counts = f"{n_fam} families, {n_prot} secreted proteins" if n_prot else f"{n_fam} families"
    title_html = f"{head} &mdash; {counts} &middot; <b>click a node</b>"
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
    return dict(families=len(NET_nodes), edges=len(NET_edges),
                bytes=len(doc), mode=mode, out=out_html)
