#!/usr/bin/env python3
"""Fail loudly when a filter expression references tags the VCF lacks.

bcftools evaluates a missing INFO/FORMAT tag as false, so a GATK filter
expression applied to DeepVariant/GLnexus output silently does nothing.
This preflight parses the configured expressions, compares them against
the VCF header, and exits non-zero with an actionable message.
"""

import argparse
import re
import subprocess
import sys

# Bare tags before a comparison: QD<2.0, MQ<40, ReadPosRankSum<-8.0, QUAL<30
TAG_RE = re.compile(r"(?<![\w.])([A-Za-z_][A-Za-z0-9_]*)\s*(?:<=|>=|==|<|>)")
FMT_RE = re.compile(r"FMT/([A-Za-z0-9_]+)")
INFO_ID_RE = re.compile(r"##INFO=<ID=([^,>]+)")
FORMAT_ID_RE = re.compile(r"##FORMAT=<ID=([^,>]+)")
BUILTIN = {"QUAL"}  # always present; not a header tag


def extract_tags(expr):
    """(info_tags, format_tags) referenced by a bcftools filter expression."""
    fmt_tags = set(FMT_RE.findall(expr))
    bare = FMT_RE.sub(" ", expr)  # strip FMT/x before scanning bare tags
    info_tags = set()
    for m in TAG_RE.finditer(bare):
        tag = m.group(1)
        if tag not in BUILTIN:
            info_tags.add(tag)
    return info_tags, fmt_tags


def parse_header(header_text):
    """(info_ids, format_ids) declared in a VCF header."""
    return set(INFO_ID_RE.findall(header_text)), set(FORMAT_ID_RE.findall(header_text))


def missing_tags(expr, info_ids, format_ids):
    info_tags, fmt_tags = extract_tags(expr)
    return sorted(info_tags - info_ids), sorted(fmt_tags - format_ids)


def check_expressions(expressions, header_text):
    """[(label, expr, missing_info, missing_fmt), ...] for every expression."""
    info_ids, format_ids = parse_header(header_text)
    problems = []
    for label, expr in expressions:
        miss_i, miss_f = missing_tags(expr, info_ids, format_ids)
        if miss_i or miss_f:
            problems.append((label, expr, miss_i, miss_f))
    return problems


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--vcf", required=True)
    ap.add_argument("--caller", default="gatk", help="active caller (for the hint)")
    ap.add_argument("--snp-filter", required=True)
    ap.add_argument("--indel-filter", required=True)
    ap.add_argument("--geno-filter", required=True)
    a = ap.parse_args(argv)

    header = subprocess.run(
        ["bcftools", "view", "-h", a.vcf],
        capture_output=True, text=True, check=True,
    ).stdout

    problems = check_expressions(
        [
            ("snp_filter", a.snp_filter),
            ("indel_filter", a.indel_filter),
            ("genotype_filter", a.geno_filter),
        ],
        header,
    )
    if not problems:
        print(f"Filter tag preflight OK: {a.vcf}")
        return 0

    print(f"Filter tag preflight FAILED for {a.vcf}:")
    for label, expr, miss_i, miss_f in problems:
        missing = [f"INFO/{t}" for t in miss_i] + [f"FORMAT/{t}" for t in miss_f]
        print(f"  {label}: missing {', '.join(missing)}")
        print(f"    expression: {expr}")
    print(
        "  bcftools treats missing tags as false, so these filters would "
        "silently do nothing.\n"
        f"  Active caller: '{a.caller}'. Use the matching filter set "
        "(tools.bcftools.<caller>) or fix the expression — GATK-only "
        "annotations (QD/MQ/FS/SOR/ReadPosRankSum) are absent from "
        "DeepVariant/GLnexus output."
    )
    return 1


if __name__ == "__main__":
    sys.exit(main())
