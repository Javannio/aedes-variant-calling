# Aedes aegypti Dengue Variant Calling Pipeline

A portable, reproducible Snakemake workflow for variant calling in *Aedes aegypti* mosquitoes from Illumina whole-genome sequencing data. Adapted from the [tb-genomics-pipeline](https://github.com/Javannio/tb-genomics-pipeline) for mosquito genomics and dengue vector surveillance.

## Features

- **Portable**: Configuration-driven via `config.yaml`, works on any system with conda
- **Scalable**: Snakemake handles parallelization, resumption, and cluster execution
- **Reproducible**: Conda environments pin exact tool versions
- **Flexible**: Accepts FASTQ files directly or downloads from SRA
- **Standard QC**: fastp trimming, BWA-MEM alignment, samtools markdup, MultiQC
- **Best-practice variant calling**: GATK4 HaplotypeCaller (GVCF) → GenotypeGVCFs → hard filtering
- **Comprehensive reporting**: Per-sample and cohort-level QC, variant statistics, MultiQC dashboard

## Pipeline Overview

```
FASTQ (input)
    │
    ▼
fastp ──► Trimmed FASTQ + QC reports
    │
    ▼
BWA-MEM ──► SAM → Sort → Fixmate → Sort → MarkDuplicates ──► Markdup BAM
    │                                                    │
    │                                                    ▼
    │                                    samtools flagstat + coverage
    │                                                    │
    ▼                                                    ▼
GATK BaseRecalibrator (optional) ──► Recalibrated BAM
    │
    ▼
GATK HaplotypeCaller (GVCF mode) ──► per-sample GVCF
    │
    ▼
GATK GenotypeGVCFs ──► Cohort VCF
    │
    ▼
GATK VariantFiltration ──► Filtered Cohort VCF
    │
    ▼
bcftools split ──► Per-sample filtered VCFs
    │
    ▼
bcftools stats + MultiQC ──► Summary reports
```

## Quick Start

### 1. Clone and Setup

```bash
git clone <this-repo> aedes-variant-calling
cd aedes-variant-calling
```

### 2. Install Conda Environment

```bash
# Create environment (first time only)
mamba env create -f config/environment.yaml
mamba activate aedes-variant-calling
```

### 3. Prepare Reference Genome

```bash
# Download and index Aedes aegypti reference (AaegL5.0)
./scripts/prepare_reference.sh ref/
```

### 4. Configure Pipeline

Copy and edit the configuration:

```bash
cp config/config.yaml config/config.local.yaml
# Edit config.local.yaml with your paths
```

Key settings to modify:
- `reference.fasta`: Path to reference FASTA
- `reference.index_prefix`: BWA index prefix
- `samples.sheet`: Path to sample sheet (or use `fastq_dir` for auto-discovery)
- `output.base_dir`: Where results go

### 5. Prepare Samples

**Option A: Sample Sheet (recommended)**

Create `config/samples.tsv`:
```tsv
sample    fastq_r1    fastq_r2
AeAeg_001    data/fastq/AeAeg_001_R1.fastq.gz    data/fastq/AeAeg_001_R2.fastq.gz
AeAeg_002    data/fastq/AeAeg_002_R1.fastq.gz    data/fastq/AeAeg_002_R2.fastq.gz
```

**Option B: Auto-discovery from FASTQ directory**

Place paired FASTQ files in a directory (naming: `*_R1_*.fastq.gz`, `*_R2_*.fastq.gz`) and set `fastq_dir` in config.

**Option C: Download from SRA**

```bash
# Create samples.txt with SRA accessions (one per line)
echo -e "SRR1234567\nSRR1234568" > config/samples.txt
./scripts/download_sra.sh config/samples.txt data/fastq
```

### 6. Run Pipeline

```bash
# Dry run (check DAG)
snakemake -n --use-conda

# Run with 16 cores
snakemake --cores 16 --use-conda

# Run on cluster (SLURM example)
snakemake --cluster "sbatch --cpus-per-task={threads} --mem={resources.mem_mb}M" \
    --jobs 100 --use-conda
```

### 7. View Results

```bash
# MultiQC dashboard
open results/qc/multiqc_report.html

# Per-sample variants
ls results/variants/*.filtered.vcf.gz

# Summary tables
cat results/variants/variant_summary.tsv
cat results/qc/alignment_summary.tsv
```

## Configuration

All parameters are in `config/config.yaml`:

| Section | Description |
|---------|-------------|
| `reference` | Reference genome FASTA, index prefix, known variants VCF |
| `samples` | Sample sheet path or FASTQ directory for auto-discovery |
| `output` | Output directory structure |
| `tools` | Tool-specific parameters (threads, memory, algorithm settings) |
| `execution` | Conda usage, dry-run, intermediate file retention |
| `qc` | Quality thresholds for reporting |
| `reporting` | MultiQC, summary tables, VCF stats |

## Requirements

- **Conda/Mamba** (for environment management)
- **Snakemake ≥ 8.0**
- **Reference genome**: Aedes aegypti AaegL5.0 (~1.3 Gb)
- **Disk space**: ~50-100 GB per sample (intermediate files)
- **RAM**: ≥ 16 GB recommended (32 GB for joint genotyping)

## Tool Versions (pinned in environment.yaml)

| Tool | Version | Purpose |
|------|---------|---------|
| fastp | 0.23.4 | Adapter trimming, QC |
| BWA | 0.7.17 | Read alignment |
| samtools | 1.19 | BAM processing, QC |
| GATK | 4.5.0.0 | Variant calling (HaplotypeCaller, GenotypeGVCFs, VariantFiltration) |
| bcftools | 1.19 | VCF manipulation, stats |
| MultiQC | 1.18 | Report aggregation |

## Output Structure

```
results/
├── aligned/
│   ├── *.markdup.bam          # Marked-duplicate BAMs
│   ├── *.recal.bam            # BQSR-recalibrated BAMs (if known variants)
│   └── *.bai                  # BAM indexes
├── variants/
│   ├── *.g.vcf.gz             # Per-sample GVCFs
│   ├── cohort.raw.vcf.gz      # Joint-called cohort VCF
│   ├── cohort.filtered.vcf.gz # Hard-filtered cohort VCF
│   ├── *.filtered.vcf.gz      # Per-sample filtered VCFs
│   ├── *.variant_stats.txt    # bcftools stats per sample
│   └── variant_summary.tsv    # Cohort variant summary
├── qc/
│   ├── fastp/                 # fastp HTML/JSON reports
│   ├── flagstat/              # samtools flagstat
│   ├── coverage/              # samtools coverage
│   ├── multiqc_report.html    # Aggregated MultiQC report
│   └── alignment_summary.tsv  # Alignment metrics table
├── trimmed/                   # Trimmed FASTQ (temp, cleaned by default)
└── logs/                      # Per-rule logs
```

## Customization

### Adding Known Variants for BQSR

If you have a known variants VCF for *Aedes aegypti* (e.g., from VectorBase or a previous study):

```yaml
reference:
  known_variants: "ref/known_variants.vcf.gz"
```

### Adjusting Variant Filters

Modify hard filter thresholds in `config.yaml`:

```yaml
tools:
  bcftools:
    filter_expr: "QUAL<30 || DP<10 || MQ<40 || FS>60.0 || SOR>3.0 || QD<2.0 || ReadPosRankSum<-8.0"
```

### Resource Scaling

For large cohorts, adjust threads/memory per rule in the Snakefile or via command line:

```bash
snakemake --cores 32 --resources mem_mb=64000 --use-conda
```

## Portability Notes

- **No hardcoded paths**: All paths configurable via `config.yaml`
- **Conda environments**: Pinned versions ensure reproducibility across systems
- **Relative paths**: Works from any working directory
- **Cluster agnostic**: Use `--cluster` for SLURM, LSF, PBS, etc.
- **Cloud ready**: Can run on AWS Batch, Google Life Sciences, etc. via Snakemake executors

## Troubleshooting

**Out of memory during GenotypeGVCFs**:
- Increase `java_mem` in config (e.g., `"32g"`)
- Use `--resources mem_mb=64000` with cluster execution

**BWA index not found**:
- Run `./scripts/prepare_reference.sh` first
- Verify `reference.index_prefix` matches index files

**No samples found**:
- Check `samples.sheet` path and format (TSV with headers: sample, fastq_r1, fastq_r2)
- Or set `samples.fastq_dir` for auto-discovery

**Variant calling slow**:
- Increase `threads` in config
- Use `--jobs` with cluster for parallel sample processing

## Citation

If you use this pipeline, please cite:
- Original TB pipeline: [tb-genomics-pipeline](https://github.com/Javannio/tb-genomics-pipeline)
- GATK: McKenna et al., Genome Research (2010)
- BWA: Li & Durbin, Bioinformatics (2009)
- fastp: Chen et al., Bioinformatics (2018)
- Snakemake: Köster & Rahmann, Bioinformatics (2012)

## License

MIT License - see LICENSE file.

## Related Work

This pipeline is a variation of the [tb-genomics-pipeline](https://github.com/Javannio/tb-genomics-pipeline) adapted for:
- *Aedes aegypti* reference genome (AaegL5.0)
- Standard GATK variant calling (instead of TB-profiler)
- Dengue vector surveillance applications
- Portable, configuration-driven design

For dengue virus phylogenetics, see the main [dengue-caribbean-phylo](https://github.com/Javannio/dengue-caribbean-phylo) pipeline.