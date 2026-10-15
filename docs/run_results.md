# Run Results — 2026-10-14

First end-to-end run of the pipeline in GPU (Parabricks) mode.

> **Note:** genotyping in this run used the since-replaced
> `pbrun genotypegvcf` + `bcftools merge` approach. The pipeline now uses
> **GLnexus** (`glnexus_cli --config DeepVariant`) for joint genotyping; the
> numbers below will be refreshed after the next full run.

## Environment

- Hardware: NVIDIA RTX 3060 (12 GB VRAM), driver 616.92, CUDA 13.4; WSL2 Ubuntu, ~15 GB RAM allocated to WSL
- Docker Desktop with WSL2 integration; image `nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1`
- Reference: AaegL5.0 (`GCF_002204515.2_AaegL5.0_genomic.fna`)
- Samples: SRR14609304, SRR14609308 (pre-trimmed paired FASTQs)

## Pipeline mode

- Alignment: `pbrun fq2bam` (GPU)
- BQSR: stub (copy + index; no known-sites resource)
- Variant calling: `pbrun deepvariant --gvcf` (GPU)
- Genotyping: `pbrun genotypegvcf` per sample, then `bcftools merge`
  (GATK joint genotyping is incompatible with DeepVariant GVCFs — no `<NON_REF>` allele)
- Filtering, splitting, stats, MultiQC: CPU (GATK/bcftools/multiqc via conda)

## Runtime

| Step | Time |
|---|---|
| fq2bam (per sample) | ~13 min |
| DeepVariant (per sample) | ~9 min |
| GenotypeGVCF + merge | ~1 min |
| Full run (16 jobs, 2 samples) | ~1 hour wall clock |

## Results

### Alignment

| Sample | Mapping rate | Properly paired |
|---|---|---|
| SRR14609304 | 99.48% | 94.33% |
| SRR14609308 | 99.58% | 94.85% |

### Variants (filtered)

Counts are identical per sample because `bcftools view -s` subsets genotypes
but keeps all cohort sites (union callset of 150,433 records).

| Metric | Value |
|---|---|
| Records | 150,433 |
| SNPs | 142,945 |
| Indels | 7,496 |
| Multiallelic sites | 222 |

## Key outputs

- `variants/cohort.raw.vcf.gz` — merged cohort VCF
- `variants/cohort.filtered.vcf.gz` — filtered cohort VCF
- `variants/SRR*.filtered.vcf.gz` — per-sample VCFs
- `variants/variant_summary.tsv` — per-sample variant counts
- `qc/alignment_summary.tsv` — mapping/pairing rates
- `qc/multiqc_report.html` — MultiQC report

## Known limitations / follow-ups

- BQSR is a stub; provide a known-sites resource (e.g. dbSNP for Ae. aegypti) for real recalibration.
- Per-sample VCFs contain all cohort sites; add `--trim-alt-alleles -c1` to `split_vcf` for sample-specific callsets.
- `bcftools merge` (used for the numbers above) produced `./.` no-calls where hom-ref was correct; GLnexus replaces it.
- Conda channel priority is not strict (`conda config --set channel_priority strict` recommended).
