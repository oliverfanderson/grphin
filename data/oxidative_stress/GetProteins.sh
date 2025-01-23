#!/bin/bash

# USAGE: In the oxidative_stress folder run `bash GetProteins.sh`.

# List of TXIDs
TXIDS=("txid6239" "txid7227" "txid7955" "txid224308" "txid559292")

# Import directory
IMPORT_DIR="$HOME/neo4j/import"

# Function to download CSVs for each TXID
download_csvs() {
    local txid=$1
    local downloads_dir="./$txid"

    echo "Getting stress proteinss for $txid..."
    if [ ! -f "GetProteins.cypher" ]; then
        echo "Error: 'GetProteins.cypher' file not found."
        exit 1
    fi
    cat GetProteins.cypher | docker exec --interactive proteinweaver cypher-shell -u neo4j

    echo "Copying downloaded CSVs for $txid from $IMPORT_DIR to $downloads_dir..."
    mkdir -p "$downloads_dir"
    cp "$IMPORT_DIR"/"$txid"-stress-proteins.csv "$downloads_dir" 2>/dev/null

    if [ $? -eq 0 ]; then
        echo "CSVs for $txid successfully copied to $downloads_dir."
    else
        echo "Error: Failed to copy CSVs for $txid from $IMPORT_DIR."
        exit 1
    fi
}

# Loop through each TXID
for txid in "${TXIDS[@]}"; do
    echo "Processing TXID: $txid"
    download_csvs "$txid"
    merge_csv "$txid"
    echo "Completed processing for TXID: $txid"
done