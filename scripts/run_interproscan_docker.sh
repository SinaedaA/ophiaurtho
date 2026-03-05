#!/bin/bash
# Wrapper script to run InterProScan via Docker
# This is called from Snakemake rules instead of directly calling InterProScan

set -e

# Parse arguments
INPUT_FILE=""
OUTPUT_FILE=""
APPLICATIONS=""
CPU=4

while [[ $# -gt 0 ]]; do
    case $1 in
        -i)
            INPUT_FILE="$2"
            shift 2
            ;;
        -o)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -appl)
            APPLICATIONS="$2"
            shift 2
            ;;
        --cpu)
            CPU="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

if [ -z "$INPUT_FILE" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Error: Input and output files are required"
    echo "Usage: $0 -i input.fasta -o output_dir -appl Pfam --cpu 4"
    exit 1
fi

# Get absolute paths
INPUT_ABS=$(readlink -f "$INPUT_FILE")
OUTPUT_DIR=$(dirname "$(readlink -f "$OUTPUT_FILE")")
OUTPUT_BASE=$(basename "$OUTPUT_FILE")

echo "Running InterProScan via Docker..."
echo "Input: $INPUT_ABS"
echo "Output: $OUTPUT_DIR/$OUTPUT_BASE"

# Run InterProScan container
# Mount input/output directories and database volume
docker run --rm \
    -v "$INPUT_ABS:/input/$(basename $INPUT_ABS)" \
    -v "$OUTPUT_DIR:/output" \
    -v ophiaurtho_interproscan_db:/opt/interproscan/data \
    interpro/interproscan:5.75-106.0 \
    /opt/interproscan/interproscan.sh \
    -i "/input/$(basename $INPUT_ABS)" \
    -d "/output" \
    -b "$OUTPUT_BASE" \
    -appl "$APPLICATIONS" \
    --cpu "$CPU" \
    -goterms \
    -iprlookup \
    -pa

echo "InterProScan complete!"
