#!/bin/bash
# download_sra.sh - Download FASTQ files from SRA accessions
# Usage: ./download_sra.sh samples.txt [output_dir]

set -euo pipefail

SAMPLES_FILE="${1:-../config/samples.txt}"
OUT_DIR="${2:-../data/fastq}"

if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file not found: $SAMPLES_FILE"
    echo "Create a file with one SRA accession per line (e.g., SRR1234567)"
    exit 1
fi

mkdir -p "$OUT_DIR"

echo "=== Downloading FASTQ from SRA ==="
echo "Samples file: $SAMPLES_FILE"
echo "Output directory: $OUT_DIR"

while IFS= read -r accession; do
    # Skip empty lines and comments
    [[ -z "$accession" || "$accession" =~ ^# ]] && continue
    
    echo "Downloading $accession..."
    fasterq-dump "$accession" \
        --outdir "$OUT_DIR" \
        --split-files \
        --threads 4 \
        --progress
    
    # Compress
    pigz "$OUT_DIR/${accession}_1.fastq" "$OUT_DIR/${accession}_2.fastq" 2>/dev/null || true
    
done < "$SAMPLES_FILE"

echo ""
echo "=== Download complete ==="
echo "FASTQ files in: $OUT_DIR"
echo ""
echo "To use with pipeline, create samples.tsv with:"
echo "sample    fastq_r1    fastq_r2"
echo "e.g.:"
echo "SRR1234567    data/fastq/SRR1234567_1.fastq.gz    data/fastq/SRR1234567_2.fastq.gz"