# Validation — 2-sample real-data run (2026-10-15)

This document is an **end-to-end pipeline test on real public data**: raw
FASTQs → fastp → Parabricks DeepVariant → GLnexus → per-caller filtering →
stats/MultiQC, plus a targeted positive control at the kdr loci. The pooled,
low-coverage, exome-capture data was chosen because it is openly available and
stress-tests the pipeline; the variant counts and genotype-based metrics below
are **pipeline diagnostics, not population-genetic results**.

## Run context

| Item | Value |
|------|-------|
| Mode | GPU — Parabricks 4.3.1 (`pbrun fq2bam` → `pbrun deepvariant --gvcf`) |
| Joint genotyping | GLnexus 1.4.1 (`DeepVariant` preset) |
| Reference | AaegL5.0, `GCF_002204515.2_AaegL5.0_genomic.fna` (~1.28 Gb, 2,310 sequences) |
| Samples | SRR14609304 + SRR14609308 (see provenance below) |
| Filtering | `bcftools` per-caller defaults (DeepVariant set): QUAL<30 per type, genotype `DP<10 \|\| GQ<20 \|\| VAF<0.2` masked to `./.` |
| Tooling | snakemake 8.15, bcftools 1.19, GATK 4.5 (unused in this run) |

## Sample provenance

Both runs are from BioProject [PRJNA731165](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA731165)
(study SRP320751, "Aedes aegypti strain:Colinas"): mosquitoes from Tapachula,
Chiapas, Mexico (collected 2018-02) exposed to permethrin in bottle bioassays.
Libraries are **exome-capture** (Illumina HiSeq 2500, 151 bp paired-end,
"RANDOM PCR" selection).

| Accession | Sample | Phenotype | Read pairs | Raw bases |
|-----------|--------|-----------|-----------:|----------:|
| SRR14609304 | Colinas pool, replicate 2 | 25 mosquitoes **dead** after permethrin | 23,828,607 | 7.20 Gbp |
| SRR14609308 | Colinas pool, replicate 2 | 25 mosquitoes **knockdown-resistant (kdr)** | 19,368,667 | 5.85 Gbp |

The local FASTQs are the **full deposited runs — no downsampling** (read
counts in flagstat, 47.7 M / 38.9 M reads, match the SRA spot counts exactly).

## Coverage and depth

| Metric | SRR14609304 | SRR14609308 |
|--------|------------:|------------:|
| Nominal depth (raw bases / 1.28 Gb) | 5.6× | 4.6× |
| Mapped | 99.48% | 99.58% |
| Duplicates | **70.0%** (33.4 M / 47.7 M) | **69.5%** (27.0 M / 38.9 M) |
| Depth after duplicate removal (`samtools coverage`) | 1.56× | 1.30× |
| Genome breadth at ≥1× (non-dup) | 15.1% | 13.7% |

The high duplication is expected for exome-capture + PCR libraries. The 15%
breadth at ~1.5× — far below the ~80% you would see at uniform 1.5× — reflects
the capture design: reads concentrate on baited targets, most of the genome is
essentially uncovered. Duplicate-marked reads are excluded by both DeepVariant
and HaplotypeCaller, so ~1.5× is the effective depth the caller saw.

## Why ~175k called variants per sample is expected (not a bug)

Typical individual *Ae. aegypti* whole-genome data (20–30×) yields millions of
variants. Three features of this dataset push the count down hard:

1. **Effective ~1.5× coverage on ~15% of the genome.** At this depth DeepVariant
   emits few confident genotypes, and the callset additionally requires
   `GQ≥20`/`DP≥10` per genotype (sites at <10× are masked to `./.`).
2. **Pools of 25 mosquitoes, not individuals.** A variant heterozygous in a
   single mosquito appears in the pool at ~2% allele frequency; a diploid
   caller can only confidently call variants common across most pool members
   (i.e. standing/lineage variants). The `VAF≥0.2` genotype mask removes the
   low-fraction tail as well.
3. **Both samples are the same strain (Colinas)** from the same experiment, so
   they share most of their variation — which is why their PASS site counts
   are nearly identical (178,805 vs 174,795 called-ALT records).

Hard filtering additionally removed 355,376 `SNP_HARDFAIL` and 47,973
`INDEL_HARDFAIL` sites — typical at low depth where QUAL/QD annotations are
unstable. If you need cohort-scale variant counts, use **individual** (not
pooled) libraries at ≥15× coverage.

## Hardware and runtimes

Hardware: Windows 11 PC, WSL2 (Ubuntu), NVIDIA GeForce RTX 3060 (12 GB,
CUDA driver 13.4 → Parabricks container CUDA 12.3), Docker Desktop engine
29.8.2, ~15 GB RAM for WSL. Image `nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1`.

Wall-clock per job on 2026-10-15 (jobs ran sequentially, snakemake `--cores 4`):

| Step | SRR14609304 | SRR14609308 |
|------|------------:|------------:|
| fastp | 1.9 min | 1.6 min |
| pbrun fq2bam (align + markdup) | 33.0 min | 27.3 min |
| pbrun deepvariant | 9.8 min | 9.0 min |
| GLnexus joint genotyping (both samples) | 6.0 min | — |
| filtering, split, stats, MultiQC | < 1 min total | |

Total successful wall-clock ≈ **1 h 50 min** end-to-end for both samples
(one additional partial `fq2bam` attempt was interrupted by a client timeout
and resumed).

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
| Ti/Tv | **1.76** (195,165 Ti / 111,084 Tv) | Healthy for insect data (typical range ~1.5–2.0); low Ti/Tv would indicate residual false positives. Exome-capture coding targets sit near the upper end, off-target reads near the lower end |
| het/hom ratio | **0.59** (het 185,362, hom-alt 128,471) | Diagnostic only: these are pools of 25, so diploid `het`/`hom` labels are the caller's model output, not individual genotypes. A pool "het" call reflects an intermediate read-level allele fraction across pool members, not Mendelian heterozygosity |
| Called variants per sample | SRR14609304: 178,805 · SRR14609308: 174,795 | Sites with a non-reference genotype for that sample (of the 329,239 PASS sites) |

For reference, the unfiltered all-site Ti/Tv is 1.59 (including hard-filter
failures) and the unfiltered het ratio is 0.62 — hard filtering removes
low-quality calls but does not distort the site spectrum.

## kdr positive control — V410L, V1016I, F1534C

Target-site check on the sodium-channel gene *para* (LOC5567355,
`NC_035109.1:315,926,360-316,405,639`, minus strand). Codon coordinates were
derived from the MSc primer pairs in `kdr_primers.csv` (in-silico PCR against
the reference: 0 mismatches, exact product sizes) intersected with the RefSeq
isoform XM_021852340.1 CDS frame; translation was verified base-for-base
against the NCBI protein XP_021708032.1. Note the numbering offset: the
literature sites correspond to RefSeq X1 residues 408 and 1012 (different
isoform N-termini/indels).

Because these are pools of 25, results are reported as **read-level AD allele
frequencies**, not diploid genotypes:

| Site | Position (NC_035109.1) | REF>ALT | Codon (transcript) | AA change | Dead pool AD → AF | kdr pool AD → AF | FILTER |
|------|------------------------|---------|--------------------|-----------|-------------------|------------------|--------|
| V410L | 316080722 | C>A | GTA>TTA | Val→Leu | 85,13 → **13.3%** | 4,56 → **93.3%** | PASS |
| V1016I | 315983763 | C>T | GTA>ATA | Val→Ile | 45,9 → **16.7%** | 0,25 → **100%** | PASS |
| F1534C | **covered — allele absent** (not a capture gap) | no ALT (REF only) | all 11 Phe codons REF | Phe retained, no Cys call | REF blocks, 51–68× | REF blocks, 38–44× | n/a (no record) |

- Both resistance alleles are strongly enriched in the knockdown-resistant
  pool versus the dead pool — the biological positive control passes.
- Site depth (AD totals) is 98×/60× (V410) and 54×/25× (V1016); mean
  amplicon depths are 77×/43× (V410), 45×/18× (V1016) for the dead/kdr
  pools — far above the ~1.5× genome-wide average, i.e. inside the capture
  targets, so these are genuine calls, not coverage artifacts.
- **F1534C: covered, allele absent — not a capture gap.** The F1534C
  amplicon has mean depths of 64× (dead pool) and 40× (kdr pool), with every
  one of the 11 candidate phenylalanine codons at 51–68× (dead) / 38–44×
  (kdr) per base, and both gVCFs are pure reference blocks across the
  amplicon. The pipeline confidently calls reference there — the Cys allele
  is not present at detectable frequency in these two Colinas pools.

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

# kdr positive control: read-level AD allele frequencies at the two sites
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%FILTER[\t%AD]\n' \
  -r NC_035109.1:316080722,NC_035109.1:315983763 variants/cohort.filtered.vcf.gz
```

## Caveats

- These are **pools of 25 mosquitoes** run through a diploid germline caller:
  `FMT/DP`, `GQ` and `VAF` genotype-level filters, het/hom ratios and
  `GT` labels are model outputs, not population-genetic quantities. The
  interpretable per-sample quantity is the **read-level AD allele frequency**
  (as used for the kdr sites above). Ti/Tv and site-level QUAL/FILTER metrics
  are ploidy-independent and remain meaningful.
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
