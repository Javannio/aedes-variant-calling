#!/usr/bin/env python3
"""Per-sample QC gates — evaluated after alignment, before variant calling.

Reads a samtools flagstat and a samtools coverage report, compares them
against the thresholds in config `qc:`, writes a PASS/FAIL report TSV and
exits non-zero when any required gate fails (after writing the report).

Metrics and thresholds (fractions are 0-1, as in config):
  mapping_rate            = mapped / total            >= qc.min_mapping_rate
  duplicate_rate          = duplicates / total        <= qc.max_duplicate_rate
  mean_depth              = len-weighted mean depth   >= qc.min_coverage_depth
  genome_fraction_covered = covered bases / length    >= qc.min_genome_fraction_covered

With qc.enabled=false every gate is written as SKIP and the script exits 0
(report-only mode).
"""

import argparse
import json
import sys


def parse_flagstat(text):
    total = mapped = duplicates = None
    for line in text.splitlines():
        if "in total" in line and total is None:
            total = int(line.split()[0])
        elif "mapped (" in line and "primary" not in line and mapped is None:
            mapped = int(line.split()[0])
        elif "duplicates" in line and "primary" not in line and duplicates is None:
            duplicates = int(line.split()[0])
    if total is None or mapped is None or duplicates is None:
        raise ValueError(
            f"Could not parse flagstat (total={total}, mapped={mapped}, "
            f"duplicates={duplicates})"
        )
    return {
        "total_reads": total,
        "mapped_reads": mapped,
        "duplicate_reads": duplicates,
        "mapping_rate": mapped / total if total else 0.0,
        "duplicate_rate": duplicates / total if total else 0.0,
    }


def parse_coverage(text):
    total_len = covered = 0
    depth_weighted = 0.0
    for line in text.splitlines():
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        start, end = int(parts[1]), int(parts[2])
        length = max(end - start + 1, 1)
        covered += int(parts[4])
        depth_weighted += float(parts[6] if len(parts) > 6 else parts[5]) * length
        total_len += length
    if total_len == 0:
        raise ValueError("Coverage report has no contigs")
    return {
        "mean_depth": depth_weighted / total_len,
        "genome_fraction_covered": covered / total_len,
    }


def evaluate_gates(metrics, gates):
    """Return [(metric, value, op, threshold, passed), ...] in report order."""
    enabled = gates.get("enabled", True)
    checks = [
        ("mapping_rate", metrics["mapping_rate"], ">=", gates.get("min_mapping_rate", 0.80)),
        ("duplicate_rate", metrics["duplicate_rate"], "<=", gates.get("max_duplicate_rate", 0.50)),
        ("mean_depth", metrics["mean_depth"], ">=", gates.get("min_coverage_depth", 10)),
        (
            "genome_fraction_covered",
            metrics["genome_fraction_covered"],
            ">=",
            gates.get("min_genome_fraction_covered", 0.80),
        ),
    ]
    rows = []
    for metric, value, op, threshold in checks:
        if not enabled:
            status = "SKIP"
        else:
            passed = value >= threshold if op == ">=" else value <= threshold
            status = "PASS" if passed else "FAIL"
        rows.append((metric, value, op, threshold, status))
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--flagstat", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--gates", required=True, help="JSON of the config qc: block")
    ap.add_argument("--out", required=True, help="Report TSV to write")
    args = ap.parse_args(argv)

    gates = json.loads(args.gates)
    with open(args.flagstat) as fh:
        metrics = parse_flagstat(fh.read())
    with open(args.coverage) as fh:
        metrics.update(parse_coverage(fh.read()))

    rows = evaluate_gates(metrics, gates)

    header = "sample\tmetric\tvalue\top\tthreshold\tstatus\n"
    body = "".join(
        f"{args.sample}\t{m}\t{v:.4f}\t{op}\t{t}\t{s}\n" for m, v, op, t, s in rows
    )
    with open(args.out, "w") as fh:
        fh.write(header + body)

    report = f"QC gate report for '{args.sample}':\n" + header + body
    print(report, end="")

    failures = [(m, v, op, t) for m, v, op, t, s in rows if s == "FAIL"]
    if failures:
        # Snakemake deletes the report TSV of a failed job, so the full
        # report must go to the log (stderr) to remain inspectable.
        detail = "; ".join(f"{m}={v:.4f} {op} {t}" for m, v, op, t in failures)
        print(
            f"QC gate FAIL for sample '{args.sample}': {detail}. "
            f"Adjust the qc: thresholds in config if intentional. "
            f"Full report above (the report TSV is kept only when all gates pass).",
            file=sys.stderr,
        )
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
