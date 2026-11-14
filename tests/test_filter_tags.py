#!/usr/bin/env python3
"""Unit tests for the filter-tag preflight.

    python tests/test_filter_tags.py
"""

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

from check_filter_tags import check_expressions, extract_tags, missing_tags, parse_header

GATK_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=QD,Number=1,Type=Float,Description="RMS QUAL">
##INFO=<ID=MQ,Number=1,Type=Float,Description="RMS MQ">
##INFO=<ID=FS,Number=1,Type=Float,Description="Fisher strand">
##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand odds ratio">
##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description="Z score">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
"""

GLNEXUS_HEADER = """\
##fileformat=VCFv4.2
##INFO=<ID=DP,Number=1,Type=Integer,Description="Depth">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
"""

GATK_SNP = "QUAL<30 || MQ<40 || QD<2.0 || FS>60.0 || SOR>3.0 || ReadPosRankSum<-8.0"
GATK_INDEL = "QUAL<30 || QD<2.0 || FS>200.0 || ReadPosRankSum<-20.0"
GATK_GENO = "FMT/DP<10 || FMT/GQ<20"
DV_SNP = "QUAL<30"
DV_GENO = (
    "FMT/DP<10 || FMT/GQ<20 || "
    "(FMT/AD[0:1]/(FMT/AD[0:0]+FMT/AD[0:1]))<0.2"
)


def main_tests():
    # 1. tag extraction from expressions
    info, fmt = extract_tags(GATK_SNP)
    assert info == {"MQ", "QD", "FS", "SOR", "ReadPosRankSum"}, info
    assert fmt == set()
    info, fmt = extract_tags(GATK_GENO)
    assert info == set() and fmt == {"DP", "GQ"}, (info, fmt)
    info, fmt = extract_tags(DV_GENO)
    assert fmt == {"DP", "GQ", "AD"}, fmt
    info, fmt = extract_tags(DV_SNP)
    assert info == set() and fmt == set(), (info, fmt)  # QUAL is builtin
    print("ok: tag extraction (bare INFO, FMT/, QUAL builtin, AD brackets)")

    # 2. GATK filters against a GATK header: no problems
    problems = check_expressions(
        [("snp_filter", GATK_SNP), ("indel_filter", GATK_INDEL),
         ("genotype_filter", GATK_GENO)],
        GATK_HEADER,
    )
    assert problems == [], problems
    print("ok: GATK filter set passes a GATK header")

    # 3. GATK filters against GLnexus/DeepVariant header: must flag the
    #    GATK-only annotations (the silent-pass bug this check closes)
    problems = check_expressions(
        [("snp_filter", GATK_SNP), ("indel_filter", GATK_INDEL),
         ("genotype_filter", GATK_GENO)],
        GLNEXUS_HEADER,
    )
    assert len(problems) == 2  # snp + indel; genotype tags exist
    labels = {p[0] for p in problems}
    assert labels == {"snp_filter", "indel_filter"}
    miss_i = set(problems[0][2])
    assert {"QD", "FS", "SOR", "ReadPosRankSum"} <= miss_i, miss_i
    assert "MQ" in miss_i
    print("ok: GATK filters fail loudly on a DeepVariant/GLnexus header")

    # 4. DeepVariant filter set passes the GLnexus header
    problems = check_expressions(
        [("snp_filter", DV_SNP), ("indel_filter", DV_SNP),
         ("genotype_filter", DV_GENO)],
        GLNEXUS_HEADER,
    )
    assert problems == [], problems
    print("ok: DeepVariant filter set passes a DeepVariant header")

    # 5. missing tags: INFO and FORMAT reported separately
    problems = check_expressions([("x", "QD<2.0 && FMT/GQ<20")], GLNEXUS_HEADER)
    assert problems and problems[0][2] == ["QD"] and problems[0][3] == [], problems
    print("ok: INFO and FORMAT missing sets reported separately")

    # 6. end-to-end via main(): write VCFs, expect exit codes
    import gzip
    from check_filter_tags import main

    tmp = Path(tempfile.mkdtemp())
    gatk_vcf = tmp / "gatk.vcf.gz"
    with gzip.open(gatk_vcf, "wt") as fh:
        fh.write(GATK_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
    dv_vcf = tmp / "dv.vcf.gz"
    with gzip.open(dv_vcf, "wt") as fh:
        fh.write(GLNEXUS_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\n")
    rc = main(["--vcf", str(gatk_vcf), "--caller", "gatk",
               "--snp-filter", GATK_SNP, "--indel-filter", GATK_INDEL,
               "--geno-filter", GATK_GENO])
    assert rc == 0, "GATK filters on GATK header should pass"
    rc = main(["--vcf", str(dv_vcf), "--caller", "deepvariant",
               "--snp-filter", GATK_SNP, "--indel-filter", GATK_INDEL,
               "--geno-filter", GATK_GENO])
    assert rc == 1, "GATK filters on GLnexus header must fail"
    rc = main(["--vcf", str(dv_vcf), "--caller", "deepvariant",
               "--snp-filter", DV_SNP, "--indel-filter", DV_SNP,
               "--geno-filter", DV_GENO])
    assert rc == 0, "DeepVariant filters on GLnexus header should pass"
    print("ok: main() exit codes (0 pass / 1 missing tags / 0 caller-correct)")

    print("\nAll filter-tag preflight tests passed.")


if __name__ == "__main__":
    main_tests()
