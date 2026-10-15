# Validation — 2-sample real-data run (2026-10-15)

## Run context

| Item | Value |
|------|-------|
| Mode | GPU — Parabricks 4.3.1 (`pbrun fq2bam` → `pbrun deepvariant --gvcf`) |
| Joint genotyping | GLnexus 1.4.1 (`DeepVariant` preset) |
| Reference | AaegL5.0, `GCF_002204515.2_AaegL5.0_genomic.fna` |
| Samples | SRR14609304, SRR14609308 (paired-end, ~90 M read pairs each) |
| Filtering | `bcftools` per-caller defaults (DeepVariant set): QUAL<30 per type, genotype `DP<10 \|\| GQ<20 \|\| VAF<0.2` masked to `./.` |
| Tooling | bcftools 1.19, snakemake 8.15 |

## Cohort callset

**Raw (GLnexus)** — `variants/cohort.raw.vcf.gz`:

| Metric | Value |
|--------|-------|
| Records | 739,006 |
| SNPs | 665,768 |
| Indels | 75,106 |
| Multiallelic sites | 6,397 |
| Ti/Tv (all sites) | 1.58 |

**After filtering** — `variants/cohort.filtered.vcf.gz`:

| FILTER | Records |
|--------|---------|
| PASS | 329,239 |
| SNP_HARDFAIL | 355,376 |
| INDEL_HARDFAIL | 47,973 |
| MONOALLELIC | 11 |
| **total** | **732,599** |

PASS set: 306,250 SNPs, 22,989 indels, 0 multiallelic sites.

## PASS-set quality metrics

| Metric | Value | Interpretation |
|--------|-------|----------------|
| Ti/Tv | **1.76** (195,165 Ti / 111,084 Tv) | Healthy for insect whole-genome data (typical range ~1.5–2.0); low Ti/Tv would indicate residual false positives |
| het/hom ratio | **0.59** (het 185,362, hom-alt 128,471) | Moderate heterozygosity, consistent with outbred wild-type mosquitoes |
| Called variants per sample | SRR14609304: 178,805 · SRR14609308: 174,795 | Sites with a non-reference genotype for that sample (of the 329,239 PASS sites) |

For reference, the unfiltered all-site Ti/Tv is 1.59 (including hard-filter
failures) and the unfiltered het ratio is 0.62 — hard filtering removes
low-quality calls but does not distort the site spectrum.

## MultiQC

`qc/multiqc_report.html` aggregates fastp, samtools flagstat/coverage and
bcftools stats per sample. In CI it is uploaded as the `multiqc-report`
artifact on every run.

## Reproducing the metrics

```bash
export PATH=$CONDA_PREFIX/bin:$PATH
cd <repo>

# Ti/Tv and counts on PASS sites only
bcftools view -f PASS variants/cohort.filtered.vcf.gz | bcftools stats - | grep -E '^(SN|TSTV)'

# het/hom ratio (PASS sites)
bcftools view -f PASS variants/cohort.filtered.vcf.gz \
  | bcftools query -f '[%GT\t]\n' | tr '\t' '\n' | grep -v '^$' | awk '
    $1=="0/1" || $1=="1/0" || $1=="0|1" || $1=="1|0" {het++}
    $1=="1/1" || $1=="1|1" {hom++}
    END {printf "het=%d hom=%d het_ratio=%.4f\n", het, hom, het/(het+hom)}'

# Per-sample called (non-ref) records — note: subset BEFORE filtering,
# bcftools view applies -i expressions before -s subsetting
for S in SRR14609304 SRR14609308; do
  printf '%s ' "$S"
  bcftools view -s "$S" -Ou variants/cohort.filtered.vcf.gz \
    | bcftools view -H -i 'GT="alt"' | wc -l
done
```

## Caveats

- Per-sample `*.filtered.vcf.gz` files keep the **full PASS site set** of the
  cohort; the genotype column belongs to that sample only. Sample-specific
  record counts therefore come from `GT="alt"` subsetting (above), not from
  `bcftools stats` on the split file, which yields identical numbers for both
  samples by construction.
- Genotypes masked by the genotype-level filter keep `./.` while the site
  remains (depth is enforced per genotype so one weak sample does not delete
  a site other samples support).
- The pre-GLnexus run from 2026-10-14 is documented separately in
  `docs/run_results.md`.
