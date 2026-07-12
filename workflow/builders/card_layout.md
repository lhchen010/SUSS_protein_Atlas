# 6-panel card layout (validated reference)

This session's finalized six-panel card layout (co_card_F*.png, dN/dS removed version):
- **A Structure** — reference structure (hub member), interactive Mol* in atlas HTML, static render on card
- **B TM Matrix** — pairwise TM within family (symmetric = min(qtm,ttm)), viridis
- **C Sequence Identity Matrix** — BLASTp pident, labeled by gene name not AF ID
- **D Per-site Trajectories** — Rate4Site conservation + ESM mutation tolerance side-by-side (replaces original dN/dS non-syn)
- **E Conserved-vs-Surface Scatter** — conservation (−Rate4Site) × rel_SASA, r value in title
- **F RNAseq Heatmap** — members × conditions (control/invitro/1dpi/3dpi/7dpi), magma

Small families (n≤3) add "n small — tree low-power" to title.
cards.py is the functional baseline; complete styled body reference this session co_card_F* v-dnds-removed.
