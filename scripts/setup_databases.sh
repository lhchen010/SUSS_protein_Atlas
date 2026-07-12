#!/usr/bin/env bash
# ============================================================
# SUSS Atlas — download the Foldseek annotation databases
# ============================================================
# Downloads the two reference databases the `annotate` step needs:
#   pdb100          (known folds; ~ a few GB)
#   afdb_swissprot  (AlphaFold/Swiss-Prot -> real protein NAMES; tens of GB)
#
# Usage:   bash scripts/setup_databases.sh /path/to/foldseek_db
# Then set in config.yaml:   tools.foldseek_db_dir: /path/to/foldseek_db
#
# Requirements: `foldseek` on PATH (conda activate suss-atlas), curl/aria2, and disk:
#   pdb100 ~2-5 GB, afdb_swissprot ~30-40 GB (download + index). Budget ~60 GB free.
# The script is idempotent: an already-built DB is skipped, so it is safe to re-run
# after an interrupted download.
# ============================================================
set -euo pipefail

DB_DIR="${1:-}"
if [[ -z "$DB_DIR" ]]; then
  echo "usage: bash scripts/setup_databases.sh /path/to/foldseek_db" >&2
  exit 1
fi
mkdir -p "$DB_DIR"

if ! command -v foldseek >/dev/null 2>&1; then
  echo "ERROR: 'foldseek' not on PATH. Run: conda activate suss-atlas" >&2
  exit 1
fi

TMP="$DB_DIR/tmp"
mkdir -p "$TMP"

download_db () {   # $1 = foldseek db name, $2 = local basename
  local name="$1" out="$DB_DIR/$2"
  if [[ -f "${out}.dbtype" || -f "${out}" ]]; then
    echo "[skip] $2 already present in $DB_DIR"
    return
  fi
  echo "[download] $name -> $out"
  foldseek databases "$name" "$out" "$TMP"
  echo "[done] $2"
}

echo "=== Foldseek annotation databases -> $DB_DIR ==="
download_db "PDB"            "pdb100"          # foldseek's PDB100 set
download_db "Alphafold/Swiss-Prot" "afdb_swissprot"

echo
echo "All set. Contents of $DB_DIR:"
ls -1 "$DB_DIR" | grep -vE '^tmp$' || true
echo
echo "Now set in config.yaml:  tools.foldseek_db_dir: $DB_DIR"
echo
echo "NOTE: InterProScan (Pfam/CDD/Gene3D) bundles its OWN data with its release download"
echo "      (Linux only) — that is separate from these Foldseek databases."
