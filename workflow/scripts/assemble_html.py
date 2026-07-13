"""assemble rule — build the interactive atlas HTML. Embeds the family network graph,
per-family 6-panel cards, embedded 3Dmol.js viewer (cartoon/surface/stick/line + B-factor
coloring for conservation/ESM/pocket), FoldTree with hub marker, matrices, RNAseq, and
per-family Excel downloads (Blob-based, artifact-sandbox safe).

html_mode:
  single  — self-contained (structures foldcomp-compressed, embedded). For <~2000 members.
  backend — metadata embedded; structures lazy-loaded from the 4070 portal. For merged
            multi-strain atlases (up to ~22,500 members).

The full validated renderer (network layout, 3Dmol viewer, EXTRA/PAY payload schema,
string-splice injection to avoid re.sub escaping, FoldTree hub star, Blob download row)
is PORTED into builders/html_builder.build_atlas(): the renderer JS is kept verbatim as
four constant template halves under builders/template/ (byte-perfect round-trip verified
against co_suss_network.html v19), and build_atlas regenerates the D + ANN data objects
from the rule outputs and splices them into that template.
"""
import os, sys, glob, json
import pandas as pd

master_csv = snakemake.input.master
cards_dir  = snakemake.input.cards
comp_xlsx  = snakemake.input.composition
anno_csv   = snakemake.input.annotation
out_html   = snakemake.output.html
mode       = snakemake.params.mode
allow_fallback = bool(snakemake.params.allow_fallback)
resdir     = os.path.dirname(out_html)
atlas_name = snakemake.config["output"]["atlas_name"]

# import the validated builder (placed under workflow/builders by the engine)
bdir = os.path.join(os.path.dirname(__file__), "..", "builders")
sys.path.insert(0, os.path.abspath(bdir))
try:
    import html_builder
    html_builder.build_atlas(
        master_csv=master_csv, cards_dir=cards_dir, composition_xlsx=comp_xlsx,
        annotation_csv=anno_csv, results_dir=resdir, out_html=out_html,
        mode=mode, atlas_name=atlas_name, config=dict(snakemake.config))
    print(f"atlas HTML built via html_builder ({mode} mode) -> {out_html}")
except Exception as e:
    if not allow_fallback:
        raise
    # Explicit debug-only fallback requested by output.allow_fallback_html.
    master = pd.read_csv(master_csv)
    ncards = len(glob.glob(os.path.join(cards_dir, "*.png")))
    rows = "\n".join(
        f"<tr><td>{r.family}</td><td>{int(r.n_members)}</td><td>{r.get('suss_pct','')}</td>"
        f"<td>{r.get('cons_sasa_r','')}</td><td>{r.get('pct_novel','')}</td></tr>"
        for _, r in master.iterrows())
    html = (f"<!doctype html><meta charset=utf-8><title>{atlas_name}</title>"
            f"<h1>{atlas_name} — SUSS Effector Atlas</h1>"
            f"<p>{len(master)} families · {ncards} cards · mode={mode}</p>"
            f"<p style='color:#a00'>NOTE: full interactive builder (html_builder.py) not found "
            f"({type(e).__name__}: {e}); this is the tabular fallback index.</p>"
            f"<table border=1 cellpadding=4><tr><th>family</th><th>n</th><th>SUSS%</th>"
            f"<th>cons-SASA r</th><th>novel%</th></tr>{rows}</table>")
    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    open(out_html, "w").write(html)
    print(f"atlas HTML fallback index -> {out_html} (builder pending: {e})")
