#!/usr/bin/env bash
set -euo pipefail

PORTAL_DIR="${SUSS_PORTAL_DIR:-/home/claude/suss_portal}"
PID_FILE="$PORTAL_DIR/suss_portal.pid"
LOG_FILE="$PORTAL_DIR/portal.log"

if [[ -f "$PID_FILE" ]]; then
    old_pid="$(cat "$PID_FILE")"
    if [[ "$old_pid" =~ ^[0-9]+$ ]] && kill -0 "$old_pid" 2>/dev/null; then
        kill "$old_pid"
        for _ in {1..20}; do
            kill -0 "$old_pid" 2>/dev/null || break
            sleep 0.25
        done
    fi
fi

export SUSS_ENGINE_TAR="$PORTAL_DIR/suss_engine.tar.gz"
export SUSS_RUNS_DIR="$PORTAL_DIR/runs"
export SUSS_CONDA="${SUSS_CONDA:-/home/claude/.conda/envs/suss}"
export SUSS_PORT="${SUSS_PORT:-8600}"
export SUSS_CORES="${SUSS_CORES:-4}"
export SUSS_BIND="${SUSS_BIND:-0.0.0.0}"
export SUSS_MAX_UPLOAD_MB="${SUSS_MAX_UPLOAD_MB:-512}"
export SUSS_MAX_EXTRACTED_MB="${SUSS_MAX_EXTRACTED_MB:-2048}"
export SUSS_MAX_ARCHIVE_FILES="${SUSS_MAX_ARCHIVE_FILES:-10000}"
export SUSS_MAX_ACTIVE_JOBS="${SUSS_MAX_ACTIVE_JOBS:-1}"

cd "$PORTAL_DIR"
setsid nohup "$SUSS_CONDA/bin/python" suss_portal.py >"$LOG_FILE" 2>&1 </dev/null &
new_pid=$!
echo "$new_pid" >"$PID_FILE"

for _ in {1..20}; do
    if curl --silent --show-error --fail "http://127.0.0.1:$SUSS_PORT/" >/dev/null; then
        echo "SUSS portal started: pid=$new_pid port=$SUSS_PORT"
        exit 0
    fi
    sleep 0.25
done

echo "SUSS portal failed its startup health check; see $LOG_FILE" >&2
exit 1
