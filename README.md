# Aedes aegypti Variant Calling Pipeline

[![CI](https://github.com/Javannio/aedes-variant-calling/actions/workflows/ci.yml/badge.svg)](https://github.com/Javannio/aedes-variant-calling/actions/workflows/ci.yml)

A portable, reproducible Snakemake workflow for germline variant calling in *Aedes aegypti* mosquitoes from Illumina whole-genome sequencing data. Adapted from the [tb-genomics-pipeline](https://github.com/Javannio/tb-genomics-pipeline) for mosquito genomics and dengue vector surveillance.

## Features

- **Portable**: Configuration-driven via `config.yaml` + optional local overrides, works on any system with conda
- **GPU-accelerated**: Optional [NVIDIA Parabricks](https://docs.nvidia.com/clara/parabricks) DeepVariant calling with GLnexus joint genotyping; GATK4 CPU fallback included
- **Reproducible**: Conda environments pin exact tool versions
- **Flexible**: Accepts FASTQ files directly or downloads from SRA
- **Standard QC**: fastp trimming (both modes), read-group-tagged alignment, duplicate marking, MultiQC
- **Best-practice hard filtering**: SNP/indel split with type-specific thresholds; per-genotype DP/GQ masking; optional repeat-exclusion BED
- **Comprehensive reporting**: Per-sample and cohort-level QC, variant statistics, MultiQC dashboard

## Pipeline Overview

```
FASTQ (input)
    │
    ▼
trimming ──► trimmed FASTQ + QC          (fastp in both modes — fq2bam does not trim)
    │
    ▼
alignment ──► Markdup BAM                (BWA-MEM + samtools markdup, or pbrun fq2bam;
    │             │                       @RG read groups set in both modes)
    │             ▼
    │      samtools flagstat + coverage
    ▼
BQSR (optional) ──► Recalibrated BAM     (skipped unless known-sites VCF provided)
    │
    ▼
variant calling ──► per-sample GVCF
    │
    ├─ GPU mode:   pbrun deepvariant --gvcf
    │                 → GLnexus (DeepVariant config; true joint genotyping)
    │                 → cohort VCF
    │
    └─ CPU mode:   gatk HaplotypeCaller --emit-ref-confidence GVCF
                      (scattered: balanced whole-contig chunks, gathered back)
                      → gatk GenomicsDBImport (contig intervals merged into
                        ~N GenomicsDB arrays — slow on fragmented refs otherwise)
                      → gatk GenotypeGVCFs (gendb:// workspace)
    │
    ▼
bcftools hard filtering ──► Filtered cohort VCF
    │   - SNPs and indels filtered separately
    │   - low-DP/low-GQ genotypes masked to ./. (not site-removed)
    │   - optional: regions in exclude_bed (e.g. RepeatMasker) skipped
    ▼
bcftools split ──► Per-sample filtered VCFs
    │
    ▼
bcftools stats + MultiQC ──► Summary reports
```

### Why two genotyping paths?

Parabricks DeepVariant GVCFs do not contain the `<NON_REF>` allele that GATK's
`GenotypeGVCFs`/`GenomicsDBImport` require, so GPU mode joint-calls with
[GLnexus](https://github.com/dnanexus-rnd/GLnexus) (the standard DeepVariant
cohort caller) instead. GPU mode is therefore true joint genotyping, not
per-sample calling + merge — genotypes at sites private to another sample are
correctly called hom-ref rather than no-call.

## Quick Start

### 1. Clone and Setup

```bash
git clone https://github.com/Javannio/aedes-variant-calling.git
cd aedes-variant-calling
```

### 2. Install Conda Environment

```bash
# Create environment (first time only)
mamba env create -f config/environment.yaml
mamba activate aedes-variant-calling
```

For GPU mode you also need Docker Desktop (with WSL2 integration) and the
Parabricks image:

```bash
docker pull nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1
```

### 3. Prepare Reference Genome

```bash
# Download and index Aedes aegypti reference (AaegL5.0)
./scripts/prepare_reference.sh ref/
```

### 4. Configure Pipeline

Defaults live in `config/config.yaml`. Put machine-specific paths in
`config/config.local.yaml` (gitignored) — it is deep-merged over the defaults:

```bash
cp config/config.local.yaml.example config/config.local.yaml
# Edit config.local.yaml with your paths
```

Key settings to modify:
- `reference.fasta`: Path to reference FASTA
- `reference.index_prefix`: BWA index prefix
- `reference.exclude_bed`: Optional BED of regions to skip (see Repeat Masking)
- `samples.sheet`: Path to sample sheet (or use `fastq_dir` for auto-discovery)
- `tools.parabricks.enabled`: `true` for GPU mode

### 5. Prepare Samples

**Option A: Sample Sheet (recommended)**

Create `config/samples.tsv`:
```tsv
sample	fastq_r1	fastq_r2
AeAeg_001	data/fastq/AeAeg_001_R1.fastq.gz	data/fastq/AeAeg_001_R2.fastq.gz
AeAeg_002	data/fastq/AeAeg_002_R1.fastq.gz	data/fastq/AeAeg_002_R2.fastq.gz
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

# Run on SLURM (Snakemake 8+: plugin-based executors)
pip install snakemake-executor-plugin-slurm
snakemake --jobs 100 --executor slurm --use-conda
```

### 7. View Results

```bash
# MultiQC dashboard
xdg-open qc/multiqc_report.html

# Per-sample variants
ls variants/*.filtered.vcf.gz

# Summary tables
cat variants/variant_summary.tsv
cat qc/alignment_summary.tsv
```

## Configuration

All parameters are in `config/config.yaml` (overridable in `config.local.yaml`).
Snakemake's own `configfile:` directives handle the deep merge, and the merged
config is validated against `workflow/schemas/config.schema.yaml` — type or
structure mistakes fail at startup instead of mid-run:

| Section | Description |
|---------|-------------|
| `reference` | Reference genome FASTA, index prefix, known variants VCF, exclude BED |
| `samples` | Sample sheet path or FASTQ directory for auto-discovery |
| `output` | Output directory structure |
| `tools` | Tool-specific parameters (threads, filters, Parabricks settings) |
| `execution` | Conda usage, dry-run, intermediate file retention |
| `qc` | Quality thresholds for reporting |

CPU-mode performance knobs (both under `tools.gatk`): `hc_scatter` controls
how many balanced whole-contig HaplotypeCaller chunks run in parallel (1 =
disabled), and `genomicsdb_num_partitions` controls how many GenomicsDB
arrays the contig intervals are merged into (0 = disabled; important for the
~1,900-contig AaegL5 reference).

## Variant Filtering

Hard filtering runs in `bcftools`, split by variant type, with **per-caller**
defaults (selected automatically from `tools.bcftools.deepvariant` vs
`tools.bcftools.gatk`):

| Type | DeepVariant (GPU) | GATK (CPU) |
|------|-------------------|------------|
| SNPs | `QUAL<30` | `QUAL<30 \|\| MQ<40 \|\| QD<2.0 \|\| FS>60.0 \|\| SOR>3.0 \|\| ReadPosRankSum<-8.0` |
| Indels | `QUAL<30` | `QUAL<30 \|\| QD<2.0 \|\| FS>200.0 \|\| ReadPosRankSum<-20.0` |
| Genotypes | `FMT/DP<10 \|\| FMT/GQ<20 \|\| VAF<0.2` → masked to `./.` | `FMT/DP<10 \|\| FMT/GQ<20` → masked to `./.` |

GATK-only annotations (`QD`, `MQ`, `FS`, `SOR`, `ReadPosRankSum`) appear only
in the GATK set because GenotypeGVCFs guarantees them; GLnexus output lacks
them and missing tags make bcftools filter expressions silently pass — so the
two callers keep separate expressions (see comments in `config.yaml`).

Note: cohort-level `INFO/DP` sums depth across all samples and scales with
cohort size — depth is therefore enforced at the **genotype** level, so a
weak sample's call is masked without deleting a site that other samples
support. Multiallelic sites are excluded from the filtered output.

## Repeat Masking and Mappability

*AaegL5* is large (~1.3 Gb), highly repetitive, and heterozygous — repeats
and low-mappability regions generate clusters of false-positive calls even
at high mapping quality. Before interpreting results:

1. Download RepeatMasker annotations for *Ae. aegypti* (VectorBase/Ensembl Metazoa)
2. Build an exclusion BED of high-copy repeats (optionally merged with
   low-mappability regions from a mappability track)
3. Set `reference.exclude_bed` in config — the filter step will skip these
   regions (`bcftools view -T ^bed`)

Without this, expect inflated SNP counts in pericentromeric and transposable
element-rich regions.

## Requirements

- **Conda/Mamba** (for environment management)
- **Snakemake ≥ 8.0**
- **GPU mode**: Docker + NVIDIA Container Toolkit + Parabricks license/image
- **Reference genome**: *Aedes aegypti* AaegL5.0 (~1.3 Gb, BWA-indexed)
- **Disk space**: ~50-100 GB per sample (intermediate files)
- **RAM**: ≥ 16 GB (32 GB+ recommended for large cohorts / GenomicsDBImport)

## Tool Versions (pinned in environment.yaml)

| Tool | Version | Purpose |
|------|---------|---------|
| fastp | 0.23.4 | Adapter trimming, QC (both modes) |
| BWA | 0.7.17 | Read alignment (CPU mode) |
| samtools | 1.19 | BAM processing, QC |
| GATK | 4.5.0.0 | HaplotypeCaller / GenotypeGVCFs (CPU mode) |
| bcftools | 1.19 | Filtering, merging, stats |
| GLnexus | 1.4.1 | Joint genotyping of DeepVariant gVCFs (GPU mode) |
| MultiQC | 1.18 | Report aggregation |
| Parabricks | 4.3.1 | fq2bam / DeepVariant (GPU mode) |

## Output Structure

```
.
├── aligned/
│   ├── *.markdup.bam          # Marked-duplicate BAMs
│   ├── *.recal.bam            # BQSR-recalibrated BAMs (if known variants)
│   └── *.bai                  # BAM indexes
├── variants/
│   ├── *.g.vcf.gz             # Per-sample GVCFs (DeepVariant/HaplotypeCaller)
│   ├── cohort.raw.vcf.gz      # Merged cohort VCF
│   ├── cohort.filtered.vcf.gz # Hard-filtered cohort VCF
│   ├── *.filtered.vcf.gz      # Per-sample filtered VCFs
│   ├── *.variant_stats.txt    # bcftools stats per sample
│   └── variant_summary.tsv    # Cohort variant summary
├── qc/
│   ├── fastp/                 # fastp HTML/JSON reports (both modes)
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

Without this, BQSR is skipped entirely — the variant callers read the
markdup BAM directly (no copy is made).

### Adjusting Variant Filters

```yaml
tools:
  bcftools:
    snp_filter: "QUAL<30"
    indel_filter: "QUAL<30"
    genotype_filter: "FMT/DP<10 || FMT/GQ<20"
```

### Resource Scaling

```bash
snakemake --cores 32 --resources mem_mb=64000 --use-conda
```

## Portability Notes

- **No hardcoded paths**: All paths configurable via config files
- **Config layering**: `config.local.yaml` deep-merges over `config.yaml`
- **Conda environments**: Pinned versions ensure reproducibility across systems
- **Cluster agnostic**: Snakemake 8 executor plugins (`--executor slurm`, lsf, etc.)
- **Cloud ready**: Any Snakemake-supported executor plugin (AWS Batch, k8s, ...)

## Testing

CI runs the **full pipeline end-to-end** on every push against a tiny
paired-end fixture dataset (`test/`: two samples, 40 read pairs each, 300 bp
reference, one planted heterozygous SNP). The run covers trimming →
alignment → marking → HaplotypeCaller → GenomicsDBImport → GenotypeGVCFs →
filtering → stats → MultiQC, and fails if the planted SNP is not called.
To reproduce locally:

```bash
bash test/make_fixtures.sh
bwa index test/ref/tiny.fa && samtools faidx test/ref/tiny.fa
gatk CreateSequenceDictionary -R test/ref/tiny.fa
cp config/config.ci.yaml config/config.local.yaml
snakemake --cores 2
```

## Provenance: Changes from tb-genomics-pipeline

This workflow started as a copy of my [TB genomics pipeline](https://github.com/Javannio/tb-genomics-pipeline)
and was rebuilt for mosquito genomics:

| Area | TB pipeline | This pipeline | Why |
|------|-------------|---------------|-----|
| Reference | *M. tuberculosis* H37Rv (~4.4 Mb) | *Ae. aegypti* AaegL5.0 (~1.3 Gb) | Different organism; 300× larger, repetitive genome |
| Calling | TB-Profiler (lineage/AMR) | DeepVariant/HaplotypeCaller germline calling | Need population variant calls, not catalogued resistance mutations |
| Joint genotyping | N/A (isolate) | GLnexus (GPU) / GenomicsDBImport (CPU) | Cohort-level analysis of related mosquitoes |
| Filtering | Minimal (TB has low diversity) | SNP/indel split, genotype-level DP, repeat exclusion | Repetitive, heterozygous genome generates false positives |
| Execution | CPU only | Optional Parabricks GPU mode | 1.3 Gb genome makes CPU calling slow |
| Config | Mixed hardcoding | Deep-merged layered config | Portability across machines |

## Validation Results

A full 2-sample run on the real AaegL5.0 reference (GPU/Parabricks + GLnexus)
passes sanity checks: **Ti/Tv 1.76** and **het/hom 0.59** on the PASS set
(329,239 sites). See [`docs/validation.md`](docs/validation.md) for all
metrics and the commands to reproduce them; the earlier pre-GLnexus run is in
[`docs/run_results.md`](docs/run_results.md).

## Troubleshooting

**Out of memory during GenotypeGVCFs / GenomicsDBImport**:
- Increase `java_mem` in config (e.g., `"32g"`)
- Use `--resources mem_mb=64000` with cluster execution

**BWA index not found**:
- Run `./scripts/prepare_reference.sh` first
- Verify `reference.index_prefix` matches index files

**No samples found**:
- Check `samples.sheet` path and format (TSV with headers: sample, fastq_r1, fastq_r2)
- Or set `samples.fastq_dir` for auto-discovery

**GPU mode: docker/gpus errors**:
- Ensure Docker Desktop WSL2 integration is enabled and `docker run --gpus all` works
- `--gpus all` belongs on `docker run`, not on `pbrun` itself

**Inflated SNP counts**:
- Set `reference.exclude_bed` (see Repeat Masking)

## Citations

If you use this pipeline, please cite:
- Original TB pipeline: [tb-genomics-pipeline](https://github.com/Javannio/tb-genomics-pipeline)
- Snakemake: Köster & Rahmann, Bioinformatics (2012); Mölder et al., eLife (2021)
- GATK: McKenna et al., Genome Research (2010); Van der Auwera & O'Connor, O'Reilly (2020)
- DeepVariant: Poplin et al., Nature Biotechnology (2018)
- NVIDIA Parabricks: https://docs.nvidia.com/clara/parabricks
- BWA: Li & Durbin, Bioinformatics (2009)
- samtools: Li et al., Bioinformatics (2009); Danecek et al., GigaScience (2021)
- bcftools: Danecek et al., GigaScience (2021)
- fastp: Chen et al., Bioinformatics (2018)
- MultiQC: Ewels et al., Bioinformatics (2016)

## License

MIT License - see [LICENSE](LICENSE) file.

## Related Work

- [tb-genomics-pipeline](https://github.com/Javannio/tb-genomics-pipeline) — the TB pipeline this was adapted from
- [dengue-caribbean-phylo](https://github.com/Javannio/dengue-caribbean-phylo) — dengue *virus* phylogenetics (this pipeline calls *mosquito* variants)
