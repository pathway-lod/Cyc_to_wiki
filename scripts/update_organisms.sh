#!/usr/bin/env bash
# scripts/update_organisms.sh
set -euo pipefail

GPML_DIR="${1:-}"
if [[ -z "$GPML_DIR" ]]; then
  echo "Usage: $0 <gpml_directory>"
  echo "Example: $0 output/plantcyc17.0.0-gpml2021__git020f9429__20260203-115037"
  exit 1
fi

if [[ ! -d "$GPML_DIR" ]]; then
  echo "Error: GPML directory does not exist: $GPML_DIR"
  exit 1
fi

GPML_DIR="${GPML_DIR%/}"                 # strip trailing slash
GPML_DIR_NAME="$(basename "$GPML_DIR")"  # human-friendly id

# Output base (separate from your GPML output pipeline)
OUT_BASE="${OUT_BASE_BRIDGEDB:-./bridgedb_output/bridgedb_organisms_update_from_plantcyc}"
mkdir -p "$OUT_BASE"

# Make a clean run dir name:
# If GPML_DIR_NAME already looks like plantcyc...__git...__YYYY..., reuse it.
# Otherwise fall back to timestamp-only.
TS="$(date +%Y%m%d-%H%M%S)"
if [[ "$GPML_DIR_NAME" == plantcyc*__git*__* ]]; then
  RUN_ID="${GPML_DIR_NAME}__bridgedb_${TS}"
else
  RUN_ID="organisms_update_${TS}"
fi

OUTDIR="${OUT_BASE%/}/${RUN_ID}"
mkdir -p "$OUTDIR"

echo "Output directory: $OUTDIR"

# 1) Download BridgeDb organisms.tsv 
echo "Downloading BridgeDb organisms.tsv (no git clone)..."
BRIDGEDB_TSV_LOCAL="${OUTDIR}/organisms_bridgedb.tsv"
BRIDGEDB_URL="https://raw.githubusercontent.com/bridgedb/datasources/main/organisms.tsv"

if command -v curl >/dev/null 2>&1; then
  curl -L -o "$BRIDGEDB_TSV_LOCAL" "$BRIDGEDB_URL"
elif command -v wget >/dev/null 2>&1; then
  wget -O "$BRIDGEDB_TSV_LOCAL" "$BRIDGEDB_URL"
else
  echo "Error: curl or wget not found."
  exit 1
fi

BRIDGEDB_SHA256="$(sha256sum "$BRIDGEDB_TSV_LOCAL" | awk '{print $1}')"
echo "BridgeDb url   : $BRIDGEDB_URL"
echo "BridgeDb sha256: $BRIDGEDB_SHA256"

# We no longer have a git commit SHA from cloning
BRIDGEDB_COMMIT="unknown"

# 2) Capture generator git commit (your repo) if available
GENERATOR_COMMIT="unknown"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  GENERATOR_COMMIT="$(git rev-parse HEAD)"
fi
echo "Generator commit: $GENERATOR_COMMIT"

# 3) Run python generator (NO gpml manifest hash)
OUT_TSV="${OUTDIR}/organisms.tsv"
OUT_META="${OUTDIR}/organisms.tsv.meta.json"

echo "Generating organisms.tsv and organisms.tsv.meta.json..."
python scripts/utils/generate_organisms_tsv.py "$GPML_DIR" \
  --existing "$BRIDGEDB_TSV_LOCAL" \
  --output "$OUT_TSV" \
  --meta-output "$OUT_META" \
  --bridgedb-commit "$BRIDGEDB_COMMIT" \
  --bridgedb-sha256 "$BRIDGEDB_SHA256" \
  --gpml-input-dir "$GPML_DIR" \
  --gpml-dir-name "$GPML_DIR_NAME" \
  --generator-git-commit "$GENERATOR_COMMIT"

# 4) Extra small run metadata file (optional but useful)
cat > "${OUTDIR}/run.metadata.txt" <<EOF
run_timestamp: $(date -Is)
gpml_input_dir: ${GPML_DIR}
gpml_input_dir_name: ${GPML_DIR_NAME}
bridgedb_git_commit: ${BRIDGEDB_COMMIT}
bridgedb_file_sha256: ${BRIDGEDB_SHA256}
generator_git_commit: ${GENERATOR_COMMIT}
EOF

# 5) Cleanup clone + downloaded TSV
rm -rf "$BRIDGEDB_DIR"
rm -f "$BRIDGEDB_TSV_LOCAL"

echo ""
echo "Done."
echo "  TSV : $OUT_TSV"
echo "  META: $OUT_META"
echo "  RUN : ${OUTDIR}/run.metadata.txt"