#!/bin/bash
# scripts/update_organisms.sh

# Directory containing GPML files (argument 1)
GPML_DIR=$1

if [ -z "$GPML_DIR" ]; then
    echo "Usage: $0 <gpml_directory>"
    echo "Example: $0 output/biocyc_pathways_20260201..."
    exit 1
fi

# Download latest organisms.tsv from BridgeDB
echo "Downloading latest organisms.tsv from BridgeDB..."
if command -v curl &> /dev/null; then
    curl -L -o organisms_bridgedb.tsv https://github.com/bridgedb/datasources/raw/main/organisms.tsv
elif command -v wget &> /dev/null; then
else
    echo "Error: curl or wget not found. Cannot download organisms.tsv."
    exit 1
fi

# Run the python script
echo "Updating organisms.tsv with local species from GPML files..."
# This script merges existing (BridgeDB) with new (local GPMLs)
python scripts/utils/generate_organisms_tsv.py "$GPML_DIR" \
    --existing organisms_bridgedb.tsv \
    --output organisms.tsv

# Clean up downloaded file
rm organisms_bridgedb.tsv

#Clean up mapping file
rm "$MAPPING_FILE"

echo "Done. Created organisms.tsv with merged data."