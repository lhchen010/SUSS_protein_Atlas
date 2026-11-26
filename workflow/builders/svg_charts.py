"""SVG matrix, heatmap, and Newick tree rendering helpers for SUSS Atlas."""

import base64
import html
import math
import re

import numpy as np


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

