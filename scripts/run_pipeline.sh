#!/usr/bin/env bash
set -euo pipefail

# Run from repo root no matter where invoked from
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# Load machine-local config (not committed)
CONFIG_FILE="${CONFIG_FILE:-scripts/config.env}"
if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
else
  echo "ERROR: Missing config file: $CONFIG_FILE"
  echo "Create it from scripts/config.example.env (or your own) and re-run."
  exit 1
fi

# Defaults if not set in config
PLANTCYC_ROOT="${PLANTCYC_ROOT:?PLANTCYC_ROOT must be set in scripts/config.env}"
GPML_SNAPSHOT="${GPML_SNAPSHOT:-gpml2021}"
OUT_BASE="${OUT_BASE:-./output_gpml}"

# Read PlantCyc version from default-version.txt (your PlantCyc layout)
DEFAULT_VER_FILE="$PLANTCYC_ROOT/default-version"
[[ -f "$DEFAULT_VER_FILE" ]] || { echo "ERROR: Not found: $DEFAULT_VER_FILE"; exit 1; }
PLANTCYC_VERSION="$(tr -d '[:space:]' < "$DEFAULT_VER_FILE")"
[[ -n "$PLANTCYC_VERSION" ]] || { echo "ERROR: Empty PlantCyc version in $DEFAULT_VER_FILE"; exit 1; }

PLANTCYC_DATA_DIR="$PLANTCYC_ROOT/$PLANTCYC_VERSION/data"
[[ -d "$PLANTCYC_DATA_DIR" ]] || { echo "ERROR: PlantCyc data dir not found: $PLANTCYC_DATA_DIR"; exit 1; }

# Git metadata for traceability
GIT_COMMIT="$(git rev-parse HEAD)"
GIT_BRANCH="$(git rev-parse --abbrev-ref HEAD)"
GIT_DIRTY="false"
if [[ -n "$(git status --porcelain)" ]]; then
  GIT_DIRTY="true"
fi

# Tag scheme (exactly your style)
TAG_NAME="plantcyc${PLANTCYC_VERSION}-${GPML_SNAPSHOT}"

# Output folder naming:
# output_gpml/plantcyc17.0.0-gpml2021__git<shortsha>__20260128-153012
RUN_TS="$(date +%Y%m%d-%H%M%S)"
SHORT_SHA="$(git rev-parse --short HEAD)"
OUT_DIR="${OUT_BASE%/}/${TAG_NAME}__git${SHORT_SHA}__${RUN_TS}"
mkdir -p "$OUT_DIR"

# Hardcoded: include reactions
CMD=(python "$REPO_ROOT/scripts/build_pathways.py" "$PLANTCYC_DATA_DIR" "$OUT_DIR" --include-reactions --db-version "$PLANTCYC_VERSION")

# Write metadata log (super useful later)
LOG_FILE="$OUT_DIR/run.metadata.txt"
{
  echo "run_timestamp: $(date -Is)"
  echo "repo: $(basename "$REPO_ROOT")"
  echo "git_branch: $GIT_BRANCH"
  echo "git_commit: $GIT_COMMIT"
  echo "git_dirty: $GIT_DIRTY"
  echo "plantcyc_root: $PLANTCYC_ROOT"
  echo "plantcyc_version: $PLANTCYC_VERSION"
  echo "plantcyc_data_dir: $PLANTCYC_DATA_DIR"
  echo "gpml_snapshot: $GPML_SNAPSHOT"
  echo "tag_name: $TAG_NAME"
  echo "output_dir: $OUT_DIR"
  echo "command: ${CMD[*]}"
  echo "python: $(python --version 2>&1)"
} > "$LOG_FILE"

echo "==> PlantCyc version: $PLANTCYC_VERSION"
echo "==> GPML snapshot:   $GPML_SNAPSHOT"
echo "==> Output dir:      $OUT_DIR"
echo "==> Log:             $LOG_FILE"

# Create tag if it doesn't exist (do nothing if it does)
if git rev-parse -q --verify "refs/tags/$TAG_NAME" >/dev/null; then
  echo "==> Tag already exists: $TAG_NAME (leaving as-is)"
else
  echo "==> Creating tag: $TAG_NAME"
  git tag -a "$TAG_NAME" -m "Data build: PlantCyc ${PLANTCYC_VERSION} -> ${GPML_SNAPSHOT}"
  echo "==> Tag created locally. (Push it with: git push origin $TAG_NAME)"
fi

# Run the pipeline
echo "==> Running: ${CMD[*]}"
"${CMD[@]}"
=======

echo "==> Done."
echo "==> Output: $OUT_DIR"