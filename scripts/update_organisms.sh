#!/bin/bash
# scripts/update_organisms.sh
set -euo pipefail

GPML_DIR="$1"
GPML_DIR_NAME="$(basename "$GPML_DIR")"

if [[ -z "$GPML_DIR" ]]; then
  echo "Usage: $0 <gpml_directory>"
  echo "Example: $0 output/biocyc_pathways_20260201_103448"
  exit 1
fi

if [[ ! -d "$GPML_DIR" ]]; then
  echo "Error: GPML directory does not exist: $GPML_DIR"
  exit 1
fi

# Timestamped output dir
TS="$(date +%Y%m%d_%H%M%S)"
OUTDIR="output/bridgedb_organisms_update_from_plantcyc/organisms_update_${TS}"
mkdir -p "$OUTDIR"

echo "Output directory: $OUTDIR"

# 1) Get BridgeDb organisms.tsv + commit SHA via shallow clone
echo "Cloning BridgeDb datasources (shallow) to capture commit SHA..."
BRIDGEDB_DIR="${OUTDIR}/bridgedb_datasources"
git clone --depth 1 https://github.com/bridgedb/datasources.git "$BRIDGEDB_DIR" >/dev/null

BRIDGEDB_COMMIT="$(git -C "$BRIDGEDB_DIR" rev-parse HEAD)"
BRIDGEDB_TSV_SRC="${BRIDGEDB_DIR}/organisms.tsv"
BRIDGEDB_TSV_LOCAL="${OUTDIR}/organisms_bridgedb.tsv"

if [[ ! -f "$BRIDGEDB_TSV_SRC" ]]; then
  echo "Error: organisms.tsv not found in BridgeDb repo clone."
  exit 1
fi

cp "$BRIDGEDB_TSV_SRC" "$BRIDGEDB_TSV_LOCAL"
BRIDGEDB_SHA256="$(sha256sum "$BRIDGEDB_TSV_LOCAL" | awk '{print $1}')"

echo "BridgeDb commit : $BRIDGEDB_COMMIT"
echo "BridgeDb sha256 : $BRIDGEDB_SHA256"

# 2) Create GPML manifest (sha256 of each gpml file, stable ordering), then hash that manifest
echo "Creating GPML manifest..."
GPML_MANIFEST="${OUTDIR}/gpml_manifest.sha256"

# sha256 each GPML file, sorted by path for determinism
find "$GPML_DIR" -type f -name "*.gpml" -print0 \
  | sort -z \
  | xargs -0 sha256sum > "$GPML_MANIFEST"

GPML_MANIFEST_SHA256="$(sha256sum "$GPML_MANIFEST" | awk '{print $1}')"
GPML_DIR_NAME="$(basename "$GPML_DIR")"

echo "GPML manifest sha256: $GPML_MANIFEST_SHA256"
echo "GPML dir name       : $GPML_DIR_NAME"

# 3) Capture generator git commit (your repo)
GENERATOR_COMMIT="unknown"
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  GENERATOR_COMMIT="$(git rev-parse HEAD)"
fi
echo "Generator commit    : $GENERATOR_COMMIT"

# 4) Run python generator
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
  --gpml-manifest-sha256 "$GPML_MANIFEST_SHA256" \
  --gpml-dir-name "$GPML_DIR_NAME" \
  --generator-git-commit "$GENERATOR_COMMIT"

# Optional: remove the clone to keep output directory smaller
rm -rf "$BRIDGEDB_DIR"

echo ""
echo "Done."
echo "  TSV : $OUT_TSV"
echo "  META: $OUT_META"
echo "  GPML manifest: $GPML_MANIFEST"