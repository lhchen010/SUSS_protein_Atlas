#!/usr/bin/env python3
"""SUSS Atlas upload portal — intranet single-machine test server (stdlib only).

Flow: browser uploads (PDB tarball + optional seqs.fasta + optional rnaseq.xlsx) +
strain metadata + parameters -> the server writes config.yaml, stages inputs into a
fresh per-job engine copy, runs the preflight `validate` rule (blocks on format errors),
then launches the full snakemake pipeline in the background. A status page polls progress
and links to the finished atlas HTML.

Runs in the `suss` conda env on the 4070. No framework — BaseHTTPRequestHandler + cgi
multipart (python 3.11). Bind 127.0.0.1 behind tailscale; NOT hardened for public exposure.

Env knobs (set before launch):
  SUSS_ENGINE_TAR   path to a suss_engine_vX.tar.gz (fresh engine per job)
  SUSS_RUNS_DIR     where per-job dirs live (default ~/suss_portal_runs)
  SUSS_CONDA        conda env prefix (default /home/claude/.conda/envs/suss)
  SUSS_PORT         listen port (default 8600)
  SUSS_BIND         bind address (default 127.0.0.1; set to the tailscale IP or 0.0.0.0
                    to reach it from a browser over tailscale — tailscale is a private
                    encrypted network, so binding its interface is the intended intranet use)
"""
import os, sys, io, cgi, json, time, shutil, tarfile, subprocess, threading, html, re
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

ENGINE_TAR = os.environ.get("SUSS_ENGINE_TAR", os.path.expanduser("~/suss_engine.tar.gz"))
RUNS_DIR   = os.environ.get("SUSS_RUNS_DIR", os.path.expanduser("~/suss_portal_runs"))
CONDA      = os.environ.get("SUSS_CONDA", "/home/claude/.conda/envs/suss")
PORT       = int(os.environ.get("SUSS_PORT", "8600"))
BIND       = os.environ.get("SUSS_BIND", "127.0.0.1")
CORES      = os.environ.get("SUSS_CORES", "4")
os.makedirs(RUNS_DIR, exist_ok=True)

JOBS = {}          # job_id -> dict(state, dir, log, atlas, msg, started, families, meta, ...)
JOBS_LOCK = threading.Lock()

# ------------------------------------------------------------------ persistence
# Each job dir holds: manifest.json (metadata+state), run.log (full log), and the
# original uploads under inputs/ — so history + inputs + params survive a restart.
def _manifest_path(jdir): return os.path.join(jdir, "manifest.json")

def _save_manifest(j):
    m = {k: j.get(k) for k in ("job_id", "state", "started", "finished", "families",
                               "npdb", "nseq", "atlas_rel", "msg", "meta")}
    try:
        with open(_manifest_path(j["dir"]), "w") as fh:
            json.dump(m, fh, indent=2)
    except Exception:
        pass

def _load_history():
    """Rebuild JOBS from manifests on disk (called at startup)."""
    if not os.path.isdir(RUNS_DIR): return
    for jid in sorted(os.listdir(RUNS_DIR)):
        jdir = os.path.join(RUNS_DIR, jid)
        if not os.path.isdir(jdir): continue
        mp = _manifest_path(jdir)
        if not os.path.isfile(mp):
            # manifest-less dir (pre-persistence run or partial) — register a minimal
            # entry so it shows up in history and can be deleted to reclaim disk
            JOBS[jid] = dict(state="unknown", dir=jdir, log=[], atlas=None, msg="no manifest",
                             started=os.path.getmtime(jdir), finished=None, families=None,
                             npdb=None, nseq=None, atlas_rel=None, meta={}, job_id=jid)
            continue
        try:
            m = json.load(open(mp))
        except Exception:
            continue
        logp = os.path.join(jdir, "run.log")
        log = open(logp, encoding="utf-8", errors="replace").read().splitlines() if os.path.exists(logp) else []
        atlas = os.path.join(jdir, "engine", m["atlas_rel"]) if m.get("atlas_rel") else None
        # a job left mid-run by a restart is no longer running: mark interrupted
        state = m.get("state", "error")
        if state in ("queued", "staging", "validating", "running"):
            state = "error"; m["msg"] = (m.get("msg") or "") + " (interrupted by server restart)"
        JOBS[jid] = dict(state=state, dir=jdir, log=log, atlas=atlas, msg=m.get("msg", ""),
                         started=m.get("started", 0), finished=m.get("finished"),
                         families=m.get("families"), npdb=m.get("npdb"), nseq=m.get("nseq"),
                         atlas_rel=m.get("atlas_rel"), meta=m.get("meta", {}), job_id=jid)

# ------------------------------------------------------------------ HTML pages
PAGE_CSS = """
body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;max-width:760px;margin:24px auto;padding:0 16px;color:#233}
h1{font-size:20px;margin:0 0 2px}.sub{color:#789;font-size:13px;margin-bottom:18px}
fieldset{border:1px solid #dce3ea;border-radius:7px;margin:0 0 14px;padding:12px 14px}
legend{font-weight:600;font-size:13px;color:#456;padding:0 6px}
label{display:block;font-size:13px;margin:7px 0 2px;color:#456}
input[type=text],select{width:100%;padding:6px 8px;border:1px solid #cdd6df;border-radius:5px;font-size:13px;box-sizing:border-box}
input[type=file]{font-size:13px;margin-top:3px}
.row{display:flex;gap:14px}.row>div{flex:1}
.toggles{display:flex;flex-wrap:wrap;gap:4px 16px}.toggles label{display:flex;align-items:center;gap:5px;margin:3px 0;font-size:13px}
.rng{display:flex;align-items:center;gap:8px}.rng input{flex:1}.rng b{min-width:44px;text-align:right;font-variant-numeric:tabular-nums}
button.go{background:#1a6db0;color:#fff;border:0;border-radius:6px;padding:10px 18px;font-size:14px;cursor:pointer}
button.go:hover{background:#155a92}.hint{color:#789;font-size:12px}
.err{background:#fdecea;border:1px solid #f5c2bd;color:#a3271b;padding:10px 12px;border-radius:6px;font-size:13px;white-space:pre-wrap}
.ok{background:#eaf6ec;border:1px solid #bfe2c6;color:#1d6b2c;padding:10px 12px;border-radius:6px;font-size:13px}
code{background:#f2f5f8;padding:1px 4px;border-radius:3px}
"""

def page(body, title="SUSS Atlas portal"):
    return (f"<!doctype html><html><head><meta charset=utf-8><title>{title}</title>"
            f"<meta name=viewport content='width=device-width,initial-scale=1'>"
            f"<style>{PAGE_CSS}</style></head><body>{body}</body></html>").encode()

def form_page(msg_html=""):
    recent = history_rows(limit=5)
    hist = (f"<fieldset><legend>Recent runs</legend><table style='font-size:12.5px;width:100%'>"
            f"<tr style='text-align:left;color:#678'><th>started</th><th>job</th><th>strain</th><th>state</th><th>detail</th><th>links</th></tr>"
            f"{recent}</table><div style='margin-top:6px'><a href=/history>full history →</a></div></fieldset>"
            ) if recent else ""
    b = f"""
    <h1>SUSS Effector Atlas — build portal</h1>
    <div class=sub>Upload a strain's secreted-protein AF2 structures; the pipeline builds an interactive SUSS atlas on this machine. Intranet test server. &middot; <a href=/history>run history</a> &middot; <a href=/help>help</a></div>
    {msg_html}
    {hist}"""
    b += f"""
    <form method=POST action=/run enctype=multipart/form-data>
      <fieldset><legend>Strain</legend>
        <div class=row>
          <div><label>Strain code (structure filename prefix)</label><input type=text name=code value=cor></div>
          <div><label>Species</label><input type=text name=species value="Colletotrichum orbiculare"></div>
        </div>
        <label>Host / host range</label><input type=text name=host_range value="Nicotiana benthamiana / Cucurbitaceae">
        <label>Colletotrichum?</label><select name=is_colleto><option value=true selected>yes</option><option value=false>no (outgroup)</option></select>
      </fieldset>
      <fieldset><legend>Input files</legend>
        <label>Structures — <b>.tar.gz of PDB files</b> (mature, signal peptide removed) <span class=hint>required</span></label>
        <input type=file name=pdb_tar accept=".gz,.tgz,.tar">
        <label>Sequences — <b>seqs.fasta</b> <span class=hint>optional; derived from structures if omitted</span></label>
        <input type=file name=seqs accept=".fasta,.fa,.faa">
        <label>RNAseq — <b>rnaseq.xlsx</b> (sheets: id_mapping, expression) <span class=hint>optional</span></label>
        <input type=file name=rnaseq accept=".xlsx">
        <div class=hint style="margin-top:3px">Format examples:
          <a href="/example?f=rnaseq_example.xlsx">filled example</a> &middot;
          <a href="/example?f=rnaseq_template_blank.xlsx">blank template</a> &middot;
          <a href="/help#rnaseq">format help</a></div>
      </fieldset>
      <fieldset><legend>Analyses to run</legend>
        <div class=toggles>
          <label><input type=checkbox checked disabled> QC</label>
          <label><input type=checkbox checked disabled> Cluster</label>
          <label><input type=checkbox name=classify checked> Classify (BLAST)</label>
          <label><input type=checkbox name=conservation checked> Conservation</label>
          <label><input type=checkbox name=pocket checked> Pocket (fpocket/P2Rank)</label>
          <label><input type=checkbox name=esm checked> ESM tolerance</label>
          <label><input type=checkbox name=foldtree checked> FoldTree</label>
          <label><input type=checkbox name=annotate checked> Annotate (InterPro/Foldseek/EffectorP)</label>
          <label><input type=checkbox name=cards checked> Cards</label>
          <label><input type=checkbox name=rnaseq_step checked> RNAseq (if xlsx given)</label>
        </div>
      </fieldset>
      <fieldset><legend>Parameters <span class=hint>(lab defaults; usually unchanged)</span></legend>
        <label>Foldseek TM threshold <span class=hint>fold-similarity cutoff</span></label>
        <div class=rng><input type=range name=tm min=0.4 max=0.7 step=0.05 value=0.5 oninput="this.nextElementSibling.textContent=this.value"><b>0.5</b></div>
        <label>Leiden resolution</label>
        <div class=rng><input type=range name=res min=0.6 max=1.6 step=0.1 value=1.0 oninput="this.nextElementSibling.textContent=this.value"><b>1.0</b></div>
        <label>Min family size</label>
        <div class=rng><input type=range name=minfam min=2 max=5 step=1 value=2 oninput="this.nextElementSibling.textContent=this.value"><b>2</b></div>
        <label>BLAST e-value (10^x)</label>
        <div class=rng><input type=range name=eexp min=-6 max=-1 step=1 value=-3 oninput="this.nextElementSibling.textContent='1e'+this.value"><b>1e-3</b></div>
        <label>Project title <span class=hint>(shown as atlas heading; blank = auto)</span></label>
        <input type=text name=project_title value="" placeholder="e.g. My Lab — C. orbiculare SUSS atlas">
      </fieldset>
      <button class=go type=submit>Validate &amp; build atlas</button>
      <div class=hint style="margin-top:8px">Files are validated for format before the pipeline runs. A malformed RNAseq workbook or empty structure set is rejected with a clear error.</div>
    </form>"""
    return page(b)

def _params_block(j):
    """Confirmed parameters + inputs used for this run (shown on the result page)."""
    m = j.get("meta") or {}
    if not m:
        return ""
    steps = m.get("steps", {})
    on = [k for k, v in steps.items() if v and k not in ("qc", "cluster", "atlas")]
    def esc(x): return html.escape(str(x))
    rows = [
        ("Strain code", m.get("code", "")),
        ("Species", m.get("species", "")),
        ("Host / range", m.get("host_range", "")),
        ("Colletotrichum", "yes" if m.get("is_colleto") else "no (outgroup)"),
        ("Foldseek TM threshold", m.get("tm")),
        ("Leiden resolution", m.get("res")),
        ("Min family size", m.get("minfam")),
        ("BLAST e-value", m.get("evalue")),
        ("Project title", m.get("project_title") or "(auto)"),
        ("RNAseq", m.get("rnaseq_mode")),
        ("Steps run", ", ".join(["qc", "cluster"] + on + ["atlas"])),
    ]
    if j.get("npdb") is not None: rows.append(("Structures used", j["npdb"]))
    if j.get("nseq"): rows.append(("Sequences provided", j["nseq"]))
    tr = "".join(f"<tr><td style='color:#678;padding:2px 10px 2px 0'>{esc(k)}</td>"
                 f"<td><b>{esc(v)}</b></td></tr>" for k, v in rows)
    # input downloads (original uploads preserved under the job dir)
    updir = os.path.join(j["dir"], "inputs")
    dls = []
    for fn, label in [("structures.tar.gz", "structures (.tar.gz)"),
                      ("seqs.fasta", "sequences (.fasta)"),
                      ("rnaseq.xlsx", "RNAseq (.xlsx)")]:
        if os.path.exists(os.path.join(updir, fn)):
            dls.append(f'<a href="/input?id={j["job_id"]}&f={fn}">{label}</a>')
    dl_html = (" &middot; ".join(dls)) if dls else "(none saved)"
    return (f"<fieldset><legend>Parameters &amp; inputs used</legend>"
            f"<table style='font-size:13px'>{tr}</table>"
            f"<div style='margin-top:8px'><b>Original inputs:</b> {dl_html}</div>"
            f"<div style='margin-top:4px'><a href='/config?id={j['job_id']}'>view full config.yaml</a></div>"
            f"</fieldset>")

def help_page():
    b = """
    <h1>SUSS Atlas portal — help</h1>
    <div class=sub><a href=/>← back to build</a> &middot; <a href=/history>run history</a></div>

    <h3>What this does</h3>
    <p style="font-size:13.5px">You upload one strain's secreted-protein AF2 structures (and optionally
    sequences and RNAseq). The pipeline clusters them into structurally-similar families (SUSS families),
    labels each family's sequence divergence, conservation, pockets, mutation tolerance, structural tree and
    expression, and builds one interactive HTML atlas. SUSS = <b>Sequence-Unrelated Structurally Similar</b>.</p>

    <h3>Inputs</h3>
    <ul style="font-size:13.5px">
      <li><b>Structures (required)</b> — a <code>.tar.gz</code> of mature AF2 PDB files (signal peptide already removed).
          Filenames should be <code>&lt;accession&gt;.pdb</code> or <code>&lt;code&gt;_&lt;accession&gt;.pdb</code>
          (e.g. <code>TDZ15389.1.pdb</code> or <code>cor_TDZ15389.1.pdb</code>). macOS <code>._</code> junk files
          in the tarball are ignored automatically.</li>
      <li><b>Sequences (optional)</b> — <code>seqs.fasta</code>. If omitted, sequences are derived from the structures.</li>
      <li><b>RNAseq (optional)</b> — <code>rnaseq.xlsx</code> (format below). Omit it and the expression panel is simply left blank.</li>
    </ul>

    <h3 id=rnaseq>RNAseq workbook format</h3>
    <p style="font-size:13.5px">Two sheets:</p>
    <ul style="font-size:13.5px">
      <li><b>id_mapping</b> — columns <code>protein_accession</code> | <code>gene_id</code>.
          Maps each protein (matching your structure accession) to its RNAseq gene/locus id.</li>
      <li><b>expression</b> — first column <code>gene_id</code>, then one column per RNAseq sample.
          Values are raw counts (or any unit); the heatmap colours by per-gene z-score across conditions.</li>
    </ul>
    <p style="font-size:13.5px">Replicates named <code>&lt;cond&gt;_1, &lt;cond&gt;_2</code> (or <code>&lt;cond&gt;.1</code>) are
    averaged; conditions are auto-ordered biologically (control/mock → in-vitro → ascending DPI); a <code>.t1</code>
    transcript suffix on gene_id is resolved automatically.<br>
    Downloads: <a href="/example?f=rnaseq_example.xlsx"><b>filled example</b></a> ·
    <a href="/example?f=rnaseq_template_blank.xlsx">blank template</a></p>

    <h3>Parameters (lab defaults)</h3>
    <ul style="font-size:13.5px">
      <li><b>Foldseek TM threshold</b> (0.5) — fold-similarity cutoff for linking two structures. 0.5 is the accepted
          fold-similarity boundary; higher = stricter (fewer, tighter families).</li>
      <li><b>Leiden resolution</b> (1.0) — clustering granularity. Higher = more, smaller families.</li>
      <li><b>Min family size</b> (2) — smallest cluster kept as a family; singletons are set aside.</li>
      <li><b>BLAST e-value</b> (1e-3) — sequence-homology sensitivity used only to <i>label</i> divergence
          (core_SUSS / diverged / moderate / recent). It never splits a structural family.</li>
    </ul>

    <h3>Analyses (toggle per run)</h3>
    <p style="font-size:13.5px">QC and clustering always run. Optional: classify (BLAST divergence), conservation
    (Rate4Site), pocket (fpocket + P2Rank), ESM tolerance (GPU), FoldTree (structural tree), annotate
    (InterPro/Foldseek/EffectorP), cards, RNAseq. Turn off what a strain doesn't need (e.g. no RNAseq, skip pocket).</p>

    <h3 id=structure>The 3D structure viewer</h3>
    <ul style="font-size:13.5px">
      <li><b>Click any residue</b> in the structure to label it with its amino-acid name and position
          (e.g. "CYS 47"); click the label again to remove it. Works in every colouring mode.</li>
      <li><b>Colour modes</b>: Conservation (Rate4Site; blue=variable→red=conserved), ESM tolerance
          (blue=constrained→red=mutation-tolerant), Pocket (grey scaffold, red pocket residues), and
          Superpose selected (overlay picked members from the tree).</li>
      <li><b>Representation</b>: cartoon / surface / stick / sphere / line, independent of colour mode.</li>
    </ul>
    <p style="font-size:13.5px" id=chimerax><b>Opening the pocket PDB in ChimeraX / PyMOL.</b>
    The "Pocket-annotated PDB" download writes the pocket-lining residues into the <b>B-factor column</b>
    (pocket residues = <code>999</code>, everything else = <code>0</code>) — because PDB has no dedicated
    "pocket" field. To see the pocket in ChimeraX:</p>
    <pre style="background:#f2f5f8;padding:9px;border-radius:5px;font-size:12.5px">open F0_fpocket_pocket.pdb
color byattribute bfactor palette lightgrey:red range 0,999
# or select just the pocket residues:
select @@bfactor=999
show sel atoms ; color sel red</pre>
    <p style="font-size:13.5px">In PyMOL: <code>spectrum b, grey_red, F0_fpocket_pocket</code> or
    <code>select pocket, b&gt;500</code>. The residue numbers also come as a plain list — use the
    <b>Pocket residues (CSV)</b> button on the card (method, residue number, amino acid).</p>

    <h3 id=numbers>Reading the numbers on a family card</h3>
    <p style="font-size:13.5px">Every metric is computed on the family's <b>reference structure</b> = the <i>hub</i>
    (the member with the highest mean TM-score to the rest — gold-starred on the FoldTree). Conservation, pockets,
    Cys and ESM all overlay on this one structure so they line up.</p>

    <p style="font-size:13.5px"><b>SUSS spectrum &amp; core_SUSS %</b></p>
    <ul style="font-size:13.5px">
      <li>Each structural link (TM ≥ threshold) between two members is labelled by how detectable their sequence
          relationship is at the chosen BLAST e-value:
          <b>core_SUSS</b> (no BLAST hit — structure conserved, sequence undetectable) →
          <b>diverged_paralog</b> (BLAST hit, &lt;30% identity) →
          <b>moderate_paralog</b> (&lt;50%) →
          <b>recent_duplicate</b> (≥50%).</li>
      <li><b>core_SUSS %</b> — the headline SUSS-ness number. Computed <i>per family</i> over its
          within-family structural links (edges), not over members:
          <div style="background:#f4f6f8;border-radius:6px;padding:9px 12px;margin:6px 0;font-family:monospace;font-size:13px">
          core_SUSS %&nbsp;=&nbsp;100 ×
          (# within-family edges labelled core_SUSS) / (# within-family edges)</div>
          where a <b>within-family edge</b> is a pair of members of that family whose structures are
          similar (Foldseek TM ≥ the clustering threshold), and an edge is <b>core_SUSS</b> when the
          same pair has <i>no</i> BLAST hit at or below the chosen e-value (structure conserved,
          sequence undetectable). Pairs with a BLAST hit fall into diverged/moderate/recent and are
          <i>not</i> counted in the numerator.</li>
      <li style="list-style:none;margin-top:4px"><span class=hint>Worked example: family F0 has 35 members
          forming 245 within-family structural edges; 226 of them have no detectable BLAST relationship
          → core_SUSS % = 100 × 226 / 245 = <b>92.2%</b>. A high % means members share a fold but their
          sequences have diverged past BLAST detection — a "true SUSS" family. A family with only 1
          member (or no structural edges) has no edges to score, so core_SUSS % is blank (—).</span></li>
    </ul>

    <p style="font-size:13.5px"><b>Annotation</b> (Annotation tab; four independent tools, no SignalP — signal
    peptides were already removed before AF2)</p>
    <ul style="font-size:13.5px">
      <li><b>Pfam / InterPro</b> (InterProScan) — sequence domains; <code>% with domain</code> and top domain shown.</li>
      <li><b>Foldseek vs PDB100 / AFDB-SwissProt</b> — closest known structure, giving a real protein <i>name</i>
          (e.g. "Pectinase", "LysM") and a TM/identity to it, rather than just an accession.</li>
      <li><b>EffectorP</b> — effector probability; <code>% effector</code> shown per family.</li>
      <li><b>DeepTMHMM</b> — transmembrane regions (n_TMR).</li>
      <li><b>novel</b> = no known fold (Foldseek) AND no domain (InterPro) — candidate lineage-specific effectors.
          <code>% novel</code> is a key column for finding new SUSS.</li>
      <li><b>fusion / multi-domain</b> — members carrying ≥2 domains (possible domain fusions).</li>
    </ul>

    <p style="font-size:13.5px"><b>"Selection &amp; sites" panel</b></p>
    <ul style="font-size:13.5px">
      <li><b>fpocket — score 0.101, 10 res</b> — fpocket's top pocket: its pocket score (higher = more
          pronounced/druggable-like cavity) and the number of pocket-lining residues on the reference.</li>
      <li><b>P2Rank — prob 0.730, 8 res</b> — P2Rank's top pocket: ligand-binding probability (0–1) and its
          residue count. fpocket and P2Rank are independent detectors; the structure viewer lets you switch
          between them.</li>
      <li><b>Cys anchors — 2</b> — number of cysteines on the reference. Cysteines form disulfide bonds that
          rigidify small secreted effectors; a high count often marks a disulfide-stabilised fold.</li>
      <li><b>Conserved-buried r — −0.20</b> — Pearson r between per-residue <i>conservation</i> and <i>relative
          SASA</i> (surface exposure). <b>Negative is expected</b>: conserved residues tend to be buried in the
          structural core, variable ones exposed on the surface. Values near 0 mean conservation and burial are
          decoupled.</li>
      <li><b>ESM vs conservation r — −0.21</b> — r between ESM mutation <i>tolerance</i> and conservation.
          <b>Negative is expected</b>: conserved positions are the least tolerant of mutation (the two signals agree).</li>
      <li><b>ESM vs SASA r — 0.02</b> — r between ESM tolerance and surface exposure. Positive would mean exposed
          sites tolerate mutation more; near 0 (as here) means exposure alone doesn't predict tolerance.</li>
    </ul>
    <p style="font-size:13.5px" class=hint>All r-values are Pearson correlations over the aligned reference
    residues; n small (2–4 member families) makes them low-power — such cards are flagged "n small".
    "conservation" = −(Rate4Site rate), so higher = more conserved.</p>

    <h3>The atlas HTML</h3>
    <p style="font-size:13.5px">The finished atlas is a <b>single self-contained HTML file</b> — every structure, plot,
    matrix and 3D viewer is embedded, with no external dependencies. Download it and:</p>
    <ul style="font-size:13.5px">
      <li>Open locally by double-clicking (works offline, any modern browser).</li>
      <li>Publish by dropping the one file onto Netlify / a lab website — nothing else needs to ship with it.</li>
    </ul>
    <p style="font-size:13.5px">Inside the atlas, each family also offers per-member sequence (FASTA) and structure (PDB) downloads.</p>

    <h3>History &amp; reproducibility</h3>
    <p style="font-size:13.5px">Every run is saved: the <a href=/history>history page</a> lists all builds (survives server
    restarts). Each result page shows the exact parameters used, links the original uploaded inputs, and the full
    <code>config.yaml</code> — so any atlas can be traced back to its inputs and settings.</p>
    <p><a href=/>← back to build</a></p>
    """
    return page(b, "SUSS portal — help")

def _fmt_time(ts):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(ts)) if ts else "—"

_STATE_BADGE = {"done": "#1d6b2c", "error": "#a3271b", "running": "#1a6db0",
                "validating": "#8a6d1a", "staging": "#8a6d1a", "queued": "#789"}

def history_rows(limit=None):
    items = sorted(JOBS.values(), key=lambda x: x.get("started", 0), reverse=True)
    if limit: items = items[:limit]
    out = []
    for j in items:
        m = j.get("meta") or {}
        st = j["state"]; col = _STATE_BADGE.get(st, "#789")
        fam = j.get("families")
        detail = (f"{fam} families" if st == "done" and fam is not None else
                  (html.escape(j.get("msg", ""))[:60] if st == "error" else st))
        out.append(
            f"<tr>"
            f"<td style='padding:3px 10px 3px 0;font-variant-numeric:tabular-nums'>{_fmt_time(j.get('started'))}</td>"
            f"<td><a href='/status?id={j['job_id']}'>{html.escape(j['job_id'])}</a></td>"
            f"<td>{html.escape(m.get('code',''))}</td>"
            f"<td><span style='color:{col};font-weight:600'>{st}</span></td>"
            f"<td class=hint>{detail}</td>"
            f"<td>" + (f"<a href='/atlas?id={j['job_id']}'>atlas</a> " if st=='done' else "") +
            f"<a href='/log?id={j['job_id']}'>log</a>"
            + ("" if st in ("running", "staging") else
               f" &middot; <form method=post action=/delete style='display:inline' "
               f"onsubmit=\"return confirm('Delete run {html.escape(j['job_id'])} and free its disk space on the server? This cannot be undone.')\">"
               f"<input type=hidden name=id value='{html.escape(j['job_id'])}'>"
               f"<button type=submit style='color:#b00;background:none;border:none;cursor:pointer;padding:0;font-size:13px;text-decoration:underline'>delete</button></form>")
            + f"</td></tr>")
    return "".join(out)

def history_page():
    rows = history_rows()
    body = ("<h1>SUSS atlas — run history</h1>"
            "<div class=sub>All builds on this server (persisted across restarts).</div>"
            "<p><a href=/>← new build</a></p>"
            "<table style='font-size:13px;width:100%'>"
            "<tr style='text-align:left;color:#678'><th>started</th><th>job</th><th>strain</th>"
            "<th>state</th><th>detail</th><th>links</th></tr>"
            + (rows or "<tr><td colspan=6 class=hint>no runs yet</td></tr>") + "</table>")
    return page(body, "SUSS run history")

def status_page(job_id):
    j = JOBS.get(job_id)
    if not j:
        return page("<h1>Unknown job</h1><p><a href=/>back</a></p>")
    st = j["state"]
    head = f"<h1>SUSS atlas build — <code>{job_id}</code></h1>"
    if st == "error":
        return page(head + f"<div class=err>Build failed:\n{html.escape(j['msg'])}</div>"
                    + _params_block(j)
                    + f"<p><a href=/>← new build</a> &middot; <a href=/log?id={job_id}>full log</a> &middot; <a href=/history>all runs</a></p>")
    if st == "done":
        fam = j.get("families", "?")
        return page(head + f"<div class=ok>Done — {fam} families. Atlas ready.</div>"
                    f"<div style='margin:14px 0'>"
                    f"<a href=/atlas?id={job_id} style='background:#1a6db0;color:#fff;padding:9px 16px;border-radius:6px;text-decoration:none;font-size:14px'>▶ Open atlas</a> "
                    f"<a href=/atlas?id={job_id}&dl=1 style='background:#1d6b2c;color:#fff;padding:9px 16px;border-radius:6px;text-decoration:none;font-size:14px;margin-left:6px'>⬇ Download atlas HTML</a>"
                    f"<a href=/summary?id={job_id} style='background:#6a4a9c;color:#fff;padding:9px 16px;border-radius:6px;text-decoration:none;font-size:14px;margin-left:6px'>⬇ Family summary (Excel)</a>"
                    f"</div>"
                    f"<div class=hint>The <b>family summary</b> Excel has one row per cluster: members, consensus annotation, "
                    f"Pfam/PDB fold, pocket residues (fpocket + P2Rank), mean structural similarity (TM), sequence identity, and SUSS %.</div>"
                    f"<div class=hint>The downloaded <code>.html</code> is fully self-contained (structures, plots and viewers all embedded, no external files). "
                    f"Drop it straight onto any web host — Netlify, a lab website, or open it locally by double-clicking. Nothing else needs to ship with it.</div>"
                    f"<p style='margin-top:10px'><a href=/log?id={job_id}>log</a> &middot; <a href=/history>all runs</a> &middot; <a href=/help>help</a></p>"
                    + _params_block(j)
                    + f"<p><a href=/>← new build</a></p>")
    # running / validating -> auto-refresh
    elapsed = int(time.time() - j["started"])
    tail = html.escape("\n".join(j["log"][-12:]))
    return page(f"<meta http-equiv=refresh content=3>" + head +
                f"<div class=sub>State: <b>{st}</b> · {elapsed}s elapsed · this page refreshes every 3s</div>"
                f"<pre style='background:#0d1117;color:#c9d1d9;padding:12px;border-radius:6px;font-size:12px;max-height:340px;overflow:auto'>{tail}</pre>"
                f"<p class=hint>You can leave this open; the build continues on the server.</p>")

# ------------------------------------------------------------------ pipeline
def _write_config(eng, meta):
    """Overlay the user's choices onto the engine's SHIPPED config.yaml (which carries all
    tool paths, the signals block, and other locked defaults) and write it back. Overwriting
    from scratch risks dropping keys the Snakefile reads (tools.p2rank_java_env, signals, …),
    so we patch rather than regenerate."""
    import yaml
    cfgp = os.path.join(eng, "config", "config.yaml")
    cfg = yaml.safe_load(open(cfgp))
    steps_in = meta["steps"]
    cfg.setdefault("strain", {}).update(
        code=meta["code"], species=meta["species"],
        host_range=meta["host_range"], is_colletotrichum=meta["is_colleto"])
    cfg.setdefault("input", {}).update(
        pdb_dir="input/pdb", seqs_fasta="input/seqs.fasta", rnaseq_xlsx=meta["rnaseq_xlsx"])
    cfg.setdefault("clustering", {}).update(
        foldseek_tm=meta["tm"], leiden_resolution=meta["res"], min_family_size=meta["minfam"])
    cfg.setdefault("classification", {}).update(blast_evalue=meta["evalue"])
    st = cfg.setdefault("steps", {})
    st.update(qc=True, cluster=True, atlas=True, rnaseq=meta["rnaseq_mode"],
              **{k: bool(v) for k, v in steps_in.items()})
    cfg.setdefault("output", {}).update(
        html_mode="single", atlas_name=f"{meta['code']}_suss_atlas",
        project_title=meta["project_title"])
    yaml.safe_dump(cfg, open(cfgp, "w"), sort_keys=False, allow_unicode=True)

def _log(j, line):
    j["log"].append(line)
    if len(j["log"]) > 4000: j["log"] = j["log"][-2000:]
    try:
        with open(os.path.join(j["dir"], "run.log"), "a", encoding="utf-8") as fh:
            fh.write(line + "\n")
    except Exception:
        pass

def run_pipeline(job_id, meta, files):
    j = JOBS[job_id]; wd = j["dir"]; eng = os.path.join(wd, "engine")
    try:
        j["state"] = "staging"; _save_manifest(j); _log(j, "extracting engine…")
        os.makedirs(eng, exist_ok=True)
        # preserve the original uploads (downloadable later from the result page)
        updir = os.path.join(wd, "inputs"); os.makedirs(updir, exist_ok=True)
        open(os.path.join(updir, "structures.tar.gz"), "wb").write(files["pdb_tar"])
        if files.get("seqs"):   open(os.path.join(updir, "seqs.fasta"), "wb").write(files["seqs"])
        if files.get("rnaseq"): open(os.path.join(updir, "rnaseq.xlsx"), "wb").write(files["rnaseq"])
        with tarfile.open(ENGINE_TAR) as t: t.extractall(eng)
        # engine tars are packed with -C suss_engine . so contents land directly in eng/
        pdbdir = os.path.join(eng, "input", "pdb"); os.makedirs(pdbdir, exist_ok=True)
        # stage PDB tarball
        _log(j, "extracting structures…")
        ptar = os.path.join(wd, "pdb_upload.tar.gz"); open(ptar, "wb").write(files["pdb_tar"])
        skipped = 0
        with tarfile.open(ptar) as t:
            for m in t.getmembers():
                base = os.path.basename(m.name)
                # skip macOS AppleDouble / metadata junk (._foo.pdb, .DS_Store, __MACOSX/)
                if base.startswith("._") or base == ".DS_Store" or "__MACOSX" in m.name:
                    skipped += 1; continue
                if m.isfile() and base.lower().endswith(".pdb"):
                    acc = base
                    if not acc.startswith(meta["code"] + "_"): acc = f"{meta['code']}_{acc}"
                    with open(os.path.join(pdbdir, acc), "wb") as fh:
                        fh.write(t.extractfile(m).read())
        if skipped: _log(j, f"skipped {skipped} macOS metadata file(s) (._*, .DS_Store)")
        npdb = len([f for f in os.listdir(pdbdir) if f.endswith(".pdb")])
        j["npdb"] = npdb; j["nseq"] = (files.get("seqs") or b"").count(b">") or None
        _log(j, f"staged {npdb} structures")
        if npdb == 0: raise ValueError("no .pdb files found in the uploaded tarball")
        # optional seqs / rnaseq
        if files.get("seqs"):
            open(os.path.join(eng, "input", "seqs.fasta"), "wb").write(files["seqs"])
            _log(j, "staged seqs.fasta")
        if files.get("rnaseq"):
            open(os.path.join(eng, "input", "rnaseq.xlsx"), "wb").write(files["rnaseq"])
            _log(j, "staged rnaseq.xlsx")
        # config (patch the shipped config so no tool/signals keys are dropped)
        _write_config(eng, meta)
        _log(j, "wrote config.yaml")

        env = dict(os.environ, PATH=f"{CONDA}/bin:" + os.environ.get("PATH", ""))
        py = f"{CONDA}/bin/python3"
        smk = f"{CONDA}/bin/snakemake"

        # ---- preflight validate (blocks on format errors) ----
        j["state"] = "validating"; _log(j, "running preflight validation…")
        vcode = ("import sys,yaml;sys.path.insert(0,'workflow/scripts');import validate_inputs as V;"
                 "e,w,i=V.validate(yaml.safe_load(open('config/config.yaml')));"
                 "print('INFO',i);print('WARN',w);\nimport sys\n"
                 "sys.exit(('PREFLIGHT_FAIL: '+' | '.join(e)) if e else 0)")
        pv = subprocess.run([py, "-c", vcode], cwd=eng, env=env,
                            capture_output=True, text=True)
        for ln in (pv.stdout + pv.stderr).splitlines(): _log(j, ln)
        if pv.returncode != 0:
            msg = next((l for l in (pv.stdout+pv.stderr).splitlines() if "PREFLIGHT_FAIL" in l),
                       "preflight validation failed")
            raise ValueError(msg.replace("PREFLIGHT_FAIL: ", ""))
        _log(j, "preflight OK")

        # ---- run pipeline ----
        j["state"] = "running"; _save_manifest(j); _log(j, f"launching snakemake (--cores {CORES})…")
        atlas_rel = f"results/{meta['code']}_suss_atlas.html"
        j["atlas_rel"] = atlas_rel
        # build the atlas AND the per-cluster Excel (family_summary.xlsx) — both explicit
        # targets, since the portal doesn't run `rule all`
        targets = [atlas_rel, "results/family_summary.xlsx"]
        proc = subprocess.Popen([smk, "--configfile", "config/config.yaml", "--cores", CORES,
                                 "-p", "--nocolor", *targets],
                                cwd=eng, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True, bufsize=1)
        for line in proc.stdout:
            _log(j, line.rstrip())
        proc.wait()
        atlas = os.path.join(eng, atlas_rel)
        if proc.returncode != 0 or not os.path.exists(atlas):
            raise RuntimeError(f"snakemake exited {proc.returncode}; atlas not produced (see log)")
        # families count from members.csv
        fam = "?"
        mem = os.path.join(eng, "results", "members.csv")
        if os.path.exists(mem):
            import csv
            fams = set()
            for r in csv.DictReader(open(mem)):
                v = r.get("family") or r.get(list(r.keys())[1])
                if v and v != "singleton": fams.add(v)
            fam = len(fams)
        j["atlas"] = atlas; j["families"] = fam; j["state"] = "done"; j["finished"] = time.time()
        _save_manifest(j)
        _log(j, f"DONE — atlas at {atlas} ({fam} families)")
    except Exception as e:
        j["state"] = "error"; j["msg"] = f"{type(e).__name__}: {e}"; j["finished"] = time.time()
        _save_manifest(j)
        _log(j, "ERROR: " + j["msg"])

# ------------------------------------------------------------------ HTTP
def _field(form, name, default=""):
    return form.getfirst(name, default)

class H(BaseHTTPRequestHandler):
    def log_message(self, *a): pass  # quiet

    def _send(self, body, code=200, ctype="text/html; charset=utf-8", extra=None):
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        for k, v in (extra or {}).items(): self.send_header(k, v)
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):
        from urllib.parse import urlparse, parse_qs
        u = urlparse(self.path); q = parse_qs(u.query)
        if u.path == "/":
            self._send(form_page())
        elif u.path == "/status":
            self._send(status_page(q.get("id", [""])[0]))
        elif u.path == "/history":
            self._send(history_page())
        elif u.path == "/help":
            self._send(help_page())
        elif u.path == "/example":
            fn = os.path.basename(q.get("f", [""])[0])
            p = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets", fn)
            if not fn or not os.path.exists(p):
                self._send(page("<h1>Example not found</h1>"), code=404); return
            self._send(open(p, "rb").read(),
                       ctype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                       extra={"Content-Disposition": f'attachment; filename="{fn}"'})
        elif u.path == "/input":
            j = JOBS.get(q.get("id", [""])[0]); fn = os.path.basename(q.get("f", [""])[0])
            p = os.path.join(j["dir"], "inputs", fn) if j else ""
            if not j or not fn or not os.path.exists(p):
                self._send(page("<h1>Input not found</h1>"), code=404); return
            ct = ("application/gzip" if fn.endswith(".gz") else
                  "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet" if fn.endswith(".xlsx")
                  else "text/plain; charset=utf-8")
            self._send(open(p, "rb").read(), ctype=ct,
                       extra={"Content-Disposition": f'attachment; filename="{j["job_id"]}_{fn}"'})
        elif u.path == "/config":
            j = JOBS.get(q.get("id", [""])[0])
            cp = os.path.join(j["dir"], "engine", "config", "config.yaml") if j else ""
            txt = open(cp, encoding="utf-8", errors="replace").read() if (cp and os.path.exists(cp)) else "config.yaml not found"
            self._send(txt.encode(), ctype="text/plain; charset=utf-8")
        elif u.path == "/log":
            j = JOBS.get(q.get("id", [""])[0])
            txt = "\n".join(j["log"]) if j else "no such job"
            self._send(txt.encode(), ctype="text/plain; charset=utf-8")
        elif u.path == "/summary":
            # downloadable per-cluster Excel (results/family_summary.xlsx) or the
            # cluster_composition.xlsx (f=composition)
            j = JOBS.get(q.get("id", [""])[0])
            which = q.get("f", ["summary"])[0]
            fn = "cluster_composition.xlsx" if which == "composition" else "family_summary.xlsx"
            p = os.path.join(j["dir"], "engine", "results", fn) if j else ""
            if not j or not os.path.exists(p):
                self._send(page("<h1>Summary not found</h1>"), code=404); return
            self._send(open(p, "rb").read(),
                       ctype="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                       extra={"Content-Disposition": f'attachment; filename="{j["job_id"]}_{fn}"'})
        elif u.path == "/atlas":
            j = JOBS.get(q.get("id", [""])[0])
            if not j or not j.get("atlas") or not os.path.exists(j["atlas"]):
                self._send(page("<h1>Atlas not ready</h1>"), code=404); return
            data = open(j["atlas"], "rb").read()
            dispo = "attachment" if q.get("dl") else "inline"
            fn = os.path.basename(j["atlas"])
            self._send(data, ctype="text/html; charset=utf-8",
                       extra={"Content-Disposition": f'{dispo}; filename="{fn}"'})
        else:
            self._send(page("<h1>404</h1>"), code=404)

    def _handle_delete(self):
        ctype = self.headers.get("Content-Type", "")
        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers,
                                environ={"REQUEST_METHOD": "POST", "CONTENT_TYPE": ctype})
        jid = os.path.basename(form.getfirst("id", ""))   # basename guards against path traversal
        j = JOBS.get(jid)
        if j and j.get("state") in ("running", "staging"):
            self._send(page("<h1>Run is still active</h1><p>Wait for it to finish or fail before deleting.</p>"
                            "<p><a href=/history>← history</a></p>"), code=409); return
        # dir from the tracked job, else fall back to an on-disk dir under RUNS_DIR
        # (older runs with no manifest aren't in JOBS but still occupy disk — clean those too)
        d = j.get("dir", "") if j else os.path.join(RUNS_DIR, jid)
        if not jid or not os.path.isdir(d):
            self._send(page("<h1>No such run</h1><p><a href=/history>← history</a></p>"), code=404); return
        removed = ""
        try:
            if d and os.path.isdir(d) and os.path.realpath(d).startswith(os.path.realpath(RUNS_DIR) + os.sep):
                shutil.rmtree(d); removed = d
        except Exception as e:
            self._send(page(f"<h1>Delete failed</h1><pre>{html.escape(str(e))}</pre>"
                            "<p><a href=/history>← history</a></p>"), code=500); return
        JOBS.pop(jid, None)
        self._send(page(f"<div class=ok>Deleted run {html.escape(jid)}"
                        + (f" (freed {html.escape(removed)})" if removed else " (no disk dir found)") + "</div>"
                        "<p><a href=/history>← back to history</a> &middot; <a href=/>new build</a></p>"))

    def do_POST(self):
        if self.path == "/delete":
            self._handle_delete(); return
        if self.path != "/run":
            self._send(page("<h1>404</h1>"), code=404); return
        ctype = self.headers.get("Content-Type", "")
        form = cgi.FieldStorage(fp=self.rfile, headers=self.headers,
                                environ={"REQUEST_METHOD": "POST", "CONTENT_TYPE": ctype})
        # required structure upload
        pdb_item = form["pdb_tar"] if "pdb_tar" in form else None
        pdb_bytes = pdb_item.file.read() if (pdb_item is not None and pdb_item.filename) else b""
        if not pdb_bytes:
            self._send(form_page("<div class=err>No structure tarball uploaded — a .tar.gz of PDB files is required.</div>"))
            return
        def rd(name):
            it = form[name] if name in form else None
            return it.file.read() if (it is not None and getattr(it, "filename", "")) else None
        files = {"pdb_tar": pdb_bytes, "seqs": rd("seqs"), "rnaseq": rd("rnaseq")}
        rnaseq_given = files["rnaseq"] is not None
        def ck(n): return form.getfirst(n) is not None
        meta = dict(
            code=re.sub(r"[^A-Za-z0-9_]", "", _field(form, "code", "strain") or "strain"),
            species=_field(form, "species", ""),
            host_range=_field(form, "host_range", ""),
            is_colleto=(_field(form, "is_colleto", "true") == "true"),
            tm=float(_field(form, "tm", "0.5")),
            res=float(_field(form, "res", "1.0")),
            minfam=int(_field(form, "minfam", "2")),
            evalue=10 ** int(_field(form, "eexp", "-3")),
            project_title=_field(form, "project_title", "").replace('"', "'"),
            rnaseq_xlsx=("input/rnaseq.xlsx" if rnaseq_given else ""),
            rnaseq_mode=(("true" if ck("rnaseq_step") else "false") if rnaseq_given else "false"),
            steps=dict(classify=ck("classify"), conservation=ck("conservation"),
                       pocket=ck("pocket"), esm=ck("esm"), foldtree=ck("foldtree"),
                       annotate=ck("annotate"), cards=ck("cards")),
        )
        job_id = time.strftime("%Y%m%d-%H%M%S") + f"-{meta['code']}"
        jdir = os.path.join(RUNS_DIR, job_id); os.makedirs(jdir, exist_ok=True)
        with JOBS_LOCK:
            JOBS[job_id] = dict(state="queued", dir=jdir, log=[], atlas=None, job_id=job_id,
                                msg="", started=time.time(), finished=None, families=None,
                                npdb=None, nseq=None, atlas_rel=None, meta=meta)
        _save_manifest(JOBS[job_id])
        threading.Thread(target=run_pipeline, args=(job_id, meta, files), daemon=True).start()
        self._send(b"", code=303, extra={"Location": f"/status?id={job_id}"})

if __name__ == "__main__":
    _load_history()
    srv = ThreadingHTTPServer((BIND, PORT), H)
    print(f"SUSS portal on http://{BIND}:{PORT}  (engine={ENGINE_TAR}, runs={RUNS_DIR}, "
          f"{len(JOBS)} prior run(s) loaded)", flush=True)
    srv.serve_forever()

