#!/bin/bash
# prepare_reference.sh - Download and prepare Aedes aegypti reference genome
# Usage: ./prepare_reference.sh [output_dir]

set -euo pipefail

OUT_DIR="${1:-../ref}"
mkdir -p "$OUT_DIR"

echo "=== Aedes aegypti Reference Genome Preparation ==="
echo "Output directory: $OUT_DIR"

# Reference: AaegL5.0 (latest high-quality assembly)
# NCBI Assembly: GCF_002204515.2
# VectorBase: AaegL5.0

REF_URL="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/204/515/GCF_002204515.2_AaegL5.0/GCF_002204515.2_AaegL5.0_genomic.fna.gz"
REF_NAME="AaegL5.0.fasta"

echo "Downloading reference genome..."
if [ ! -f "$OUT_DIR/$REF_NAME" ]; then
    wget -q --show-progress -O "$OUT_DIR/$REF_NAME.gz" "$REF_URL"
    gunzip -c "$OUT_DIR/$REF_NAME.gz" > "$OUT_DIR/$REF_NAME"
    rm "$OUT_DIR/$REF_NAME.gz"
else
    echo "Reference already exists, skipping download"
fi

echo "Creating BWA index..."
if [ ! -f "$OUT_DIR/$REF_NAME.bwt" ]; then
    bwa index -p "$OUT_DIR/AaegL5.0" "$OUT_DIR/$REF_NAME"
else
    echo "BWA index already exists"
fi

echo "Creating samtools faidx..."
if [ ! -f "$OUT_DIR/$REF_NAME.fai" ]; then
    samtools faidx "$OUT_DIR/$REF_NAME"
else
    echo "faidx already exists"
fi

echo "Creating sequence dictionary (for GATK)..."
if [ ! -f "$OUT_DIR/AaegL5.0.dict" ]; then
    gatk --java-options "-Xmx4g" CreateSequenceDictionary -R "$OUT_DIR/$REF_NAME" -O "$OUT_DIR/AaegL5.0.dict"
else
    echo "Sequence dictionary already exists"
fi

echo ""
echo "=== Reference preparation complete ==="
echo "Reference FASTA: $OUT_DIR/$REF_NAME"
echo "BWA index prefix: $OUT_DIR/AaegL5.0"
echo ""
echo "Update config/config.yaml with:"
echo "  reference:"
echo "    fasta: \"ref/$REF_NAME\""
echo "    index_prefix: \"ref/AaegL5.0\""