#!/bin/bash

# Function to display usage
usage() {
    echo "Usage: $0 [-a | -s txid1,txid2,... | -n 1000]"
    echo "  -a       Run GRPhIN in graphlets_only mode for all species."
    echo "  -s       Run GRPhIN in graphlets_only mode for specific species using comma-separated TXIDs (e.g. 7227,7955)"
    echo "  -n       The number of randomized networks to run GRPhIN on."
    exit 1
}

# Parse command-line arguments
ALL_TXIDS=true
SPECIFIC_TXIDS=""
NUM_NETWORKS=50

while getopts ":as:n:" opt; do
    case $opt in
        a)
            ALL_TXIDS=true
            ;;
        s)
            ALL_TXIDS=false
            SPECIFIC_TXIDS=$OPTARG
            ;;
        n)
            NUM_NETWORKS=$OPTARG
            ;;
        \?)
            usage
            ;;
    esac
done

# Function to run GRPhIN by TXID
grphin_by_txid() {
    local txid=$1
    local i=$2
    local ppi="data/oxidative_stress/$txid/randomized_networks/stress_ppi$i.csv"
    local reg="data/oxidative_stress/$txid/randomized_networks/stress_reg$i.csv"
    local out="data/oxidative_stress/$txid/randomized_networks"
    echo "Running GRPhIN on TXID: $txid, Network: $i..."
    python3 grphin.py -u "$ppi" -d "$reg" -o "$out" -g True
    echo "TXID: $txid, Network: $i finished."
}

# Loop over randomized networks
for ((i=0; i<=(NUM_NETWORKS-1); i++)); do
    echo "Running GRPhIN on randomized network $i..."
    if [ "$ALL_TXIDS" = true ]; then
        for txid in 7227 224308 7955 6239 559292; do
            grphin_by_txid "txid$txid" "$i"
        done
    else
        IFS=',' read -ra TXID_ARRAY <<< "$SPECIFIC_TXIDS"
        for txid in "${TXID_ARRAY[@]}"; do
            grphin_by_txid "txid$txid" "$i"
        done
    fi
done