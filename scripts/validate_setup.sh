#!/bin/bash
# validate_setup.sh - Check pipeline setup before running

set -euo pipefail

CONFIG_FILE="${1:-config/config.yaml}"

echo "=== Pipeline Setup Validation ==="
echo "Config file: $CONFIG_FILE"
echo ""

# Check config file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Config file not found: $CONFIG_FILE"
    exit 1
fi

# Load config (using yq if available, otherwise basic checks)
echo "Checking configuration..."

# Check reference
REF_FASTA=$(grep -A1 "reference:" "$CONFIG_FILE" | grep "fasta:" | head -1 | sed 's/.*: *//' | sed 's/"//g')
if [ -n "$REF_FASTA" ] && [ "$REF_FASTA" != "" ]; then
    if [ -f "$REF_FASTA" ]; then
        echo "  ✓ Reference FASTA exists: $REF_FASTA"
    else
        echo "  ✗ Reference FASTA not found: $REF_FASTA"
        echo "    Run: ./scripts/prepare_reference.sh"
    fi
fi

REF_INDEX=$(grep -A2 "reference:" "$CONFIG_FILE" | grep "index_prefix:" | head -1 | sed 's/.*: *//' | sed 's/"//g')
if [ -n "$REF_INDEX" ] && [ "$REF_INDEX" != "" ]; then
    if [ -f "${REF_INDEX}.bwt" ]; then
        echo "  ✓ BWA index exists: ${REF_INDEX}.bwt"
    else
        echo "  ✗ BWA index not found: ${REF_INDEX}.bwt"
        echo "    Run: ./scripts/prepare_reference.sh"
    fi
fi

# Check samples
SAMPLES_SHEET=$(grep -A1 "samples:" "$CONFIG_FILE" | grep "sheet:" | head -1 | sed 's/.*: *//' | sed 's/"//g')
FASTQ_DIR=$(grep -A2 "samples:" "$CONFIG_FILE" | grep "fastq_dir:" | head -1 | sed 's/.*: *//' | sed 's/"//g')

if [ -n "$SAMPLES_SHEET" ] && [ "$SAMPLES_SHEET" != "" ] && [ -f "$SAMPLES_SHEET" ]; then
    echo "  ✓ Sample sheet exists: $SAMPLES_SHEET"
    # Count samples
    N_SAMPLES=$(tail -n +2 "$SAMPLES_SHEET" | grep -v "^#" | wc -l)
    echo "    Samples in sheet: $N_SAMPLES"
    # Check FASTQ files exist
    while IFS=$'\t' read -r sample r1 r2; do
        [[ "$sample" =~ ^# ]] && continue
        [[ -z "$sample" ]] && continue
        if [ ! -f "$r1" ]; then
            echo "    ✗ Missing: $r1"
        fi
        if [ ! -f "$r2" ]; then
            echo "    ✗ Missing: $r2"
        fi
    done < "$SAMPLES_SHEET"
elif [ -n "$FASTQ_DIR" ] && [ "$FASTQ_DIR" != "" ] && [ -d "$FASTQ_DIR" ]; then
    echo "  ✓ FASTQ directory exists: $FASTQ_DIR"
    N_R1=$(find "$FASTQ_DIR" -name "*_R1_*.fastq.gz" | wc -l)
    echo "    Paired FASTQ pairs found: $N_R1"
else
    echo "  ✗ No valid sample source configured"
    echo "    Set samples.sheet or samples.fastq_dir in config"
fi

# Check conda
echo ""
echo "Checking conda environment..."
if command -v mamba &> /dev/null; then
    echo "  ✓ mamba found"
elif command -v conda &> /dev/null; then
    echo "  ✓ conda found"
else
    echo "  ✗ Neither mamba nor conda found"
    echo "    Install miniforge/mambaforge first"
fi

# Check snakemake
if command -v snakemake &> /dev/null; then
    echo "  ✓ snakemake found: $(snakemake --version)"
else
    echo "  ✗ snakemake not in PATH (will use --use-conda)"
fi

# Check tools in conda env
echo ""
echo "Checking conda environment 'aedes-variant-calling'..."
if mamba env list | grep -q "aedes-variant-calling"; then
    echo "  ✓ Environment exists"
    # Quick tool checks
    mamba run -n aedes-variant-calling fastp --version 2>&1 | head -1 | sed 's/^/    /'
    mamba run -n aedes-variant-calling bwa 2>&1 | head -1 | sed 's/^/    /'
    mamba run -n aedes-variant-calling samtools --version 2>&1 | head -1 | sed 's/^/    /'
    mamba run -n aedes-variant-calling gatk --version 2>&1 | head -1 | sed 's/^/    /'
    mamba run -n aedes-variant-calling bcftools --version 2>&1 | head -1 | sed 's/^/    /'
    mamba run -n aedes-variant-calling multiqc --version 2>&1 | head -1 | sed 's/^/    /'
else
    echo "  ✗ Environment not found"
    echo "    Run: mamba env create -f config/environment.yaml"
fi

echo ""
echo "=== Validation Complete ==="
echo ""
echo "Next steps:"
echo "  1. Fix any ✗ items above"
echo "  2. Run dry-run: snakemake -n --use-conda"
echo "  3. Run pipeline: snakemake --cores 16 --use-conda"