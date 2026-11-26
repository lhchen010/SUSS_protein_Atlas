#!/usr/bin/env bash
# Build a clean, macOS-metadata-free SUSS engine tarball for portal deployment or release.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
VERSION="$(tr -d "[:space:]" < "$REPO_ROOT/VERSION")"
COMMIT="$(git -C "$REPO_ROOT" rev-parse HEAD 2>/dev/null || echo "unknown")"
SHORT="$(echo "$COMMIT" | cut -c1-7)"
OUT="${1:-$REPO_ROOT/suss_engine-v${VERSION}-${SHORT}.tar.gz}"

STAGE="$(mktemp -d)"
trap "rm -rf \"$STAGE\"" EXIT

git -C "$REPO_ROOT" archive HEAD | tar -x -C "$STAGE"
printf "%s\n" "$COMMIT" > "$STAGE/GIT_COMMIT"

# Strip any stray AppleDouble or cache artifacts defensively
find "$STAGE" -name "._*" -o -name ".DS_Store" -o -name "__pycache__" -o -name ".pytest_cache" | xargs -r rm -rf

COPYFILE_DISABLE=1 tar -czf "$OUT" -C "$STAGE" .
echo "Built clean release archive: $OUT ($(du -h "$OUT" | awk "{print \$1}"))"
