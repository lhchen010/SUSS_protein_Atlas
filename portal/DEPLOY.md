# SUSS Atlas upload portal — 4070 deploy

Intranet single-machine test server that connects the upload page to the Snakemake engine:
browser upload → format validation → pipeline run → interactive atlas HTML.

## What it does
1. Serves an upload form (strain metadata, PDB tarball + optional seqs.fasta / rnaseq.xlsx,
   step toggles, parameter sliders, editable project title).
2. On submit: stages inputs into a fresh per-job engine copy, patches the engine's shipped
   `config.yaml` with the user's choices (keeps all tool paths / signals defaults), runs the
   preflight `validate` rule (rejects malformed inputs with a clear message), then runs the
   full pipeline in the background.
3. A status page auto-refreshes with the live snakemake log and links to the finished atlas.
4. Every run is persisted: `<run>/manifest.json` (state + confirmed parameters), `<run>/run.log`
   (full log), and `<run>/inputs/` (the original uploaded structures.tar.gz / seqs.fasta /
   rnaseq.xlsx). The result page shows a **Parameters & inputs used** table with download links
   and a link to the exact `config.yaml`. A **/history** page lists all runs and survives restarts
   (the server rebuilds JOBS from the manifests at startup; a run left mid-build by a restart is
   marked errored/interrupted).

## Run on the 4070
```bash
# files live in /home/claude/suss_portal/ : suss_portal.py, suss_engine.tar.gz, launch.sh
bash /home/claude/suss_portal/launch.sh      # PID-scoped restart + local health check
```
Environment knobs (set in launch.sh): `SUSS_ENGINE_TAR`, `SUSS_RUNS_DIR`, `SUSS_CONDA`,
`SUSS_PORT` (8600), `SUSS_CORES` (4), `SUSS_BIND` (0.0.0.0 for tailscale reach).

## Reach it
From a machine on the tailscale network: **http://100.80.77.29:8600**
(4070 = `mpfi-linux1.tailf711c0.ts.net`). Per-job outputs under `SUSS_RUNS_DIR`.

## Scope / caveats (intranet test build)
- Runs one job per submit in a background thread; no auth, no queue, no user isolation.
  Fine for tailscale-internal lab use; NOT hardened for public exposure.
- Binds `0.0.0.0` so the tailscale interface is reachable — safe only because the 4070's
  public interfaces are not exposed and tailscale is a private encrypted network.
- `cgi` module is deprecated (removed in py3.13); the engine's `suss` env is py3.11 so it
  works today. If moving to py3.13, swap the multipart parse for `multipart`/`email` or a
  small framework.
- To take it down, send `TERM` to the PID stored in `/home/claude/suss_portal/suss_portal.pid`.
