#!/usr/bin/env bash
set -euo pipefail

PYTHON_BIN="${SUSS_DEEPTMHMM_PYTHON:-/home/user/miniconda3/bin/python}"
COMPAT_DIR="${SUSS_DEEPTMHMM_COMPAT_DIR:-/home/claude/suss_tools/deeptmhmm-compat}"
USER_SITE="$($PYTHON_BIN -c 'import site; print(site.getusersitepackages())')"

export PYTHONNOUSERSITE=1
export PYTHONPATH="$COMPAT_DIR:$USER_SITE${PYTHONPATH:+:$PYTHONPATH}"
exec "$PYTHON_BIN" "$@"
