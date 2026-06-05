#!/usr/bin/env bash
# run_pipeline.sh — thin config wrapper around build_pathways.py
#
# All build logic (output naming, validation, metadata, species coverage,
# GPML statistics) lives in build_pathways.py. This script only:
#   1. Reads config.env for machine-local paths
#   2. Resolves the PlantCyc version
#   3. Optionally creates a git tag for the release
#   4. Calls build_pathways.py with the right arguments
#
# Usage:
#   ./scripts/run_pipeline.sh
#   RELEASE_VERSION=v2 ./scripts/run_pipeline.sh   # override release tag
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

# ── Load machine-local config ─────────────────────────────────────────────────
CONFIG_FILE="${CONFIG_FILE:-scripts/config.env}"
if [[ -f "$CONFIG_FILE" ]]; then
  # shellcheck disable=SC1090
  source "$CONFIG_FILE"
else
  echo "ERROR: Missing config file: $CONFIG_FILE"
  echo "Create it from scripts/config.example.env and re-run."
  exit 1
fi

# ── Resolve PlantCyc paths ────────────────────────────────────────────────────
PLANTCYC_ROOT="${PLANTCYC_ROOT:?PLANTCYC_ROOT must be set in scripts/config.env}"
OUT_BASE="${OUT_BASE:-./output_gpml}"
RELEASE_VERSION="${RELEASE_VERSION:-v1}"

DEFAULT_VER_FILE="$PLANTCYC_ROOT/default-version"
[[ -f "$DEFAULT_VER_FILE" ]] || { echo "ERROR: Not found: $DEFAULT_VER_FILE"; exit 1; }
PLANTCYC_VERSION="$(tr -d '[:space:]' < "$DEFAULT_VER_FILE")"
[[ -n "$PLANTCYC_VERSION" ]] || { echo "ERROR: Empty version in $DEFAULT_VER_FILE"; exit 1; }

PLANTCYC_DATA_DIR="$PLANTCYC_ROOT/$PLANTCYC_VERSION/data"
[[ -d "$PLANTCYC_DATA_DIR" ]] || { echo "ERROR: Data dir not found: $PLANTCYC_DATA_DIR"; exit 1; }

# ── Optional release tag ──────────────────────────────────────────────────────
# Tag scheme: plantcyc17.0.0-gpml2021-v1
# Increment RELEASE_VERSION in config.env when re-releasing the same
# PlantCyc+GPML version with updated pipeline code.
TAG_NAME="plantcyc${PLANTCYC_VERSION}-gpml2021-${RELEASE_VERSION}"

if git rev-parse -q --verify "refs/tags/$TAG_NAME" >/dev/null; then
  echo "==> Tag already exists: $TAG_NAME (leaving as-is)"
else
  echo "==> Creating tag: $TAG_NAME"
  git tag -a "$TAG_NAME" -m "Data build: PlantCyc ${PLANTCYC_VERSION} -> gpml2021 (${RELEASE_VERSION})"
  echo "==> Tag created locally. Push with: git push origin $TAG_NAME"
fi

# ── Run the pipeline ──────────────────────────────────────────────────────────
# build_pathways.py handles:
#   - input validation (validate_plantcyc_input.py)
#   - output directory naming: plantcyc{version}-gpml2021__git{sha}__{ts}
#   - run.metadata.txt
#   - GPML_STATISTICS_REPORT.txt
#   - VALIDATION_REPORT.txt + VALIDATION_SUMMARY.tsv
#   - species_coverage.tsv + species_coverage_by_ncbi.tsv
echo "==> PlantCyc version : $PLANTCYC_VERSION"
echo "==> Output base dir  : $OUT_BASE"
echo ""

python "$REPO_ROOT/scripts/build_pathways.py" \
  "$PLANTCYC_DATA_DIR" \
  "$OUT_BASE" \
  --include-reactions \
  --db-version "$PLANTCYC_VERSION"

echo ""
echo "==> Done."
