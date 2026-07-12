"""cards rule — one 6-panel PNG per family (A structure ref | B TM matrix |
C sequence-identity matrix | D per-site conservation+ESM track | E cons-vs-SASA
scatter | F RNAseq heatmap). Small families (n<=3) get a "n small — tree low-power"
caption. Renders from master.csv + per-family signature/matrix/expression files.
This is a functional baseline; the fully-styled 6-panel body validated this session
(co_card_F*.png v-dnds-removed) is the reference layout — see builders/card_layout.md.
"""
import os, glob, re
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

master_csv = snakemake.input.master
outdir     = snakemake.output[0]
resdir     = os.path.dirname(master_csv)
famdir     = os.path.join(resdir, "families")
os.makedirs(outdir, exist_ok=True)

master = pd.read_csv(master_csv)
for _, row in master.iterrows():
    fam = row.family; n = int(row.n_members)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8))
    fig.suptitle(f"{fam}  (n={n})  SUSS {row.get('suss_pct','?')}%  "
                 f"conservation↔surface r={row.get('cons_sasa_r','?')}"
                 + ("   [n small — tree low-power]" if n <= 3 else ""),
                 fontsize=12, x=0.02, ha="left")
    A, B, C, D, E, F = axes.ravel()
    # A — reference structure placeholder note (interactive Mol* lives in HTML)
    A.text(0.5, 0.5, f"ref structure\n(see atlas HTML\nfor 3D viewer)", ha="center", va="center")
    A.set_title("A · structure"); A.axis("off")
    D.set_title("D · per-site")   # default; overwritten below when signature/ESM present
    # E — cons vs SASA scatter
    sig = os.path.join(famdir, fam, f"{fam}_signature.csv")
    if os.path.exists(sig):
        m = pd.read_csv(sig).dropna(subset=["rel_sasa","conservation"])
        E.scatter(m.conservation, m.rel_sasa, s=8, alpha=0.6)
        E.set_xlabel("conservation (−Rate4Site)"); E.set_ylabel("rel. SASA")
        # D — per-site track: Rate4Site conservation + optional ESM variant-tolerance
        resi = list(m.resi) if "resi" in m else list(range(len(m)))
        D.plot(resi, m.conservation, lw=0.8, color="#2a6b8a", label="conservation")
        D.set_xlabel("residue"); D.set_ylabel("conservation", color="#2a6b8a")
        # overlay ESM mean per-residue LLR, aligned to the SAME residue index as conservation
        esm_plotted = False
        esm_all = os.path.join(os.path.dirname(famdir), "esm_all.csv")
        if os.path.exists(esm_all):
            try:
                ea = pd.read_csv(esm_all)
                ef = ea[ea.family == fam] if "family" in ea.columns else ea.iloc[0:0]
                aa_cols = [c for c in ef.columns if len(str(c)) == 1 and str(c).isalpha()]
                if len(ef) and aa_cols:
                    mean_llr = ef[aa_cols].mean(axis=1).values
                    # align to conservation residues: ESM row i -> residue resi[i]
                    n = min(len(mean_llr), len(resi))
                    if n:
                        D2 = D.twinx()
                        D2.plot(resi[:n], mean_llr[:n], lw=0.8, color="#c0562a",
                                alpha=0.7, label="ESM tolerance")
                        D2.set_ylabel("ESM mean LLR", color="#c0562a")
                        esm_plotted = True
            except Exception:
                esm_plotted = False
        D.set_title("D · per-site (conservation + ESM)" if esm_plotted else "D · per-site conservation")
    E.set_title("E · conserved vs surface")
    # B / C — TM & identity matrices
    for ax, key, lab in [(B, f"{fam}_TM.csv", "B · TM"), (C, f"{fam}_ID.csv", "C · seq identity")]:
        mp = os.path.join(famdir, fam, key)
        if os.path.exists(mp):
            mat = pd.read_csv(mp, index_col=0)
            im = ax.imshow(mat.values, cmap="viridis", aspect="auto"); fig.colorbar(im, ax=ax, fraction=0.046)
        ax.set_title(lab)
    # F — RNAseq
    exp = os.path.join(famdir, fam, f"{fam}_expression.csv")
    if os.path.exists(exp):
        e = pd.read_csv(exp, index_col=0)
        num = e.select_dtypes("number")
        if len(num): im = F.imshow(num.values, cmap="magma", aspect="auto"); fig.colorbar(im, ax=F, fraction=0.046)
    F.set_title("F · RNAseq")
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(outdir, f"co_card_{fam}.png"), dpi=140)
    plt.close(fig)
print(f"cards: {len(master)} family cards -> {outdir}")

