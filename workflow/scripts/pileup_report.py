#!/usr/bin/env python3
"""Read-level allele counting for the pool_seq branch and the targeted-sites
report.

Subcommands
  counts   Genome-wide candidate allele counts for a pool BAM: streams
           `samtools mpileup` and emits positions with non-ref support
           >= --min-alt-reads and AF >= --min-alt-af (pools have no
           meaningful diploid GT, so counts are the product).
  targets  Targeted-sites report: for each configured target position,
           counts ref/alt/other bases from two mpilemp passes (all reads,
           and MAPQ >= --min-mapq-pass) so sites are always reported with
           and without ambiguous (MAPQ 0) reads — callers can emit
           reference-only there (F1534C lesson).

Both disable BAQ (-B), remove overlapping mate bases (-x) and filter bases
by quality (-Q), matching the kdr validation methodology.
"""

import argparse
import gzip
import json
import os
import subprocess
import sys
import tempfile

ACGT = ("A", "C", "G", "T")


def parse_pileup_bases(ref, bases):
    """Count bases from an mpileup base string.

    Handles ^X read-start anchors, $ ends, +/-N indels, . , ref matches,
    * # deletions, < > reference skips and N bases.
    """
    counts = {b: 0 for b in ACGT}
    counts.update({"ref": 0, "N": 0, "del": 0, "skip": 0})
    i, n = 0, len(bases)
    while i < n:
        c = bases[i]
        if c == "^":  # read start: next char is the MAPQ offset
            i += 2
            continue
        if c == "$":  # read end
            i += 1
            continue
        if c in "+-":  # indel: sign, length digits, then the inserted bases
            j = i + 1
            while j < n and bases[j].isdigit():
                j += 1
            if j == i + 1:  # malformed token: treat the sign as data
                i += 1
                continue
            length = int(bases[i + 1 : j])
            i = j + length
            continue
        if c in ".,=":
            counts["ref"] += 1
        elif c in "ACGT":
            counts[c] += 1
        elif c in "acgt":
            counts[c.upper()] += 1
        elif c in "Nn":
            counts["N"] += 1
        elif c in "*#":
            counts["del"] += 1
        elif c in "<>":
            counts["skip"] += 1
        i += 1
    return counts


def parse_pileup_line(line):
    """(chrom, pos, ref, depth, counts) from one mpileup line."""
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 5:
        raise ValueError(f"Malformed mpileup line: {line!r}")
    chrom, pos, ref, depth, bases = parts[0], int(parts[1]), parts[2].upper(), int(parts[3]), parts[4]
    return chrom, pos, ref, depth, parse_pileup_bases(ref, bases)


def emit_counts_row(chrom, pos, ref, depth, counts, min_alt_reads, min_alt_af):
    """Row tuple for a candidate polymorphic position, or None."""
    acgt_total = counts["ref"] + sum(counts[b] for b in ACGT)
    if acgt_total == 0:
        return None
    alt = {b: counts[b] for b in ACGT if b != ref}
    max_alt_base, max_alt_count = max(alt.items(), key=lambda kv: kv[1])
    if max_alt_count < min_alt_reads:
        return None
    af = max_alt_count / acgt_total
    if af < min_alt_af:
        return None
    return (
        chrom, pos, ref, acgt_total, counts["ref"],
        counts["A"], counts["C"], counts["G"], counts["T"],
        counts["N"], counts["del"],
        max_alt_base, max_alt_count, round(af, 6),
    )


def run_mpileup(args):
    """Run samtools mpileup; returns {(chrom, pos): (ref, depth, counts)}."""
    proc = subprocess.run(
        ["samtools", "mpileup", "-x", "-B"] + args,
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(
            f"samtools mpileup failed (exit {proc.returncode}): "
            f"{' '.join(proc.args)}\n{proc.stderr.strip()}"
        )
    rows = {}
    for line in proc.stdout.splitlines():
        if not line:
            continue
        chrom, pos, ref_base, depth, counts = parse_pileup_line(line)
        rows[(chrom, pos)] = (ref_base, depth, counts)
    return rows


def load_targets(path_or_json):
    if path_or_json.startswith("["):
        return json.loads(path_or_json)
    with open(path_or_json) as fh:
        return json.load(fh)


def assemble_target_rows(sample, targets, fai_lengths, q0, q1):
    """Build report rows for every target from both mpileup passes."""
    rows = []
    for t in targets:
        chrom, pos, ref, alt = t["chrom"], int(t["pos"]), t["ref"].upper(), t["alt"].upper()
        base = [sample, t["name"], chrom, str(pos), ref, alt]
        if chrom not in fai_lengths:
            rows.append(base + ["0"] * 8 + ["REF_MISSING"])
            continue
        if pos > fai_lengths[chrom]:
            rows.append(base + ["0"] * 8 + ["OUT_OF_RANGE"])
            continue
        r0, d0, c0 = q0.get((chrom, pos), (None, 0, None))
        r1, d1, c1 = q1.get((chrom, pos), (None, 0, None))
        if d0 == 0 or c0 is None:
            rows.append(base + ["0"] * 8 + ["NO_COVERAGE"])
            continue
        if r0 and r0 != ref:
            rows.append(base + ["0"] * 8 + ["REF_MISMATCH"])
            continue
        acgt_total = c0["ref"] + sum(c0[b] for b in ACGT)
        ref_count = c0["ref"]
        alt_count = c0[alt]
        other = acgt_total - ref_count - alt_count
        af_alt = alt_count / acgt_total if acgt_total else 0.0
        af_total = alt_count / d0 if d0 else 0.0
        callable_acgt = 0
        if c1 is not None:
            callable_acgt = c1["ref"] + sum(c1[b] for b in ACGT)
        mapq0_frac = 1.0 - (callable_acgt / acgt_total) if acgt_total else 0.0
        if callable_acgt == 0 and acgt_total > 0:
            status = "MAPQ0_ONLY"
        else:
            status = "OK"
        rows.append(
            base
            + [str(d0), str(acgt_total), str(ref_count), str(alt_count), str(other),
               f"{af_alt:.6f}", f"{af_total:.6f}", f"{mapq0_frac:.4f}", status]
        )
    return rows


TARGET_HEADER = (
    "sample\ttarget\tchrom\tpos\tref\talt\tdepth_total\tdepth_acgt\t"
    "ref_count\talt_count\tother_count\taf_alt\taf_alt_incl_ambiguous\t"
    "mapq0_fraction\tstatus\n"
)

COUNTS_HEADER = (
    "chrom\tpos\tref\tdepth_acgt\tref_count\tcount_A\tcount_C\tcount_G\tcount_T\t"
    "count_N\tcount_del\tmax_alt_base\tmax_alt_count\tmax_alt_af\n"
)


def cmd_counts(a):
    proc = subprocess.Popen(
        ["samtools", "mpileup", "-x", "-B", "-q", "0", "-Q", str(a.min_bq),
         "-d", str(a.max_depth), "-f", a.ref, a.bam],
        stdout=subprocess.PIPE, text=True,
    )
    opener = gzip.open if a.out.endswith(".gz") else open
    n_in = n_out = 0
    with opener(a.out, "wt") as fh:
        fh.write(COUNTS_HEADER)
        for line in proc.stdout:
            n_in += 1
            chrom, pos, ref, depth, counts = parse_pileup_line(line)
            row = emit_counts_row(chrom, pos, ref, depth, counts,
                                  a.min_alt_reads, a.min_alt_af)
            if row:
                fh.write("\t".join(str(x) for x in row) + "\n")
                n_out += 1
    proc.stdout.close()
    rc = proc.wait()
    if rc != 0:
        raise SystemExit(f"samtools mpileup failed with exit code {rc}")
    print(f"[{a.sample}] scanned {n_in} positions, emitted {n_out} candidates -> {a.out}")


def cmd_targets(a):
    targets = load_targets(a.targets)
    fai_lengths = {}
    with open(a.fai) as fh:
        for line in fh:
            parts = line.split("\t")
            fai_lengths[parts[0]] = int(parts[1])
    known = [t for t in targets if t["chrom"] in fai_lengths
             and int(t["pos"]) <= fai_lengths[t["chrom"]]]
    with tempfile.NamedTemporaryFile("w", suffix=".bed", delete=False) as bed:
        for t in known:
            pos = int(t["pos"])
            bed.write(f"{t['chrom']}\t{pos - 1}\t{pos}\n")
        bed_path = bed.name
    try:
        common = ["-Q", str(a.min_bq), "-d", str(a.max_depth),
                  "-f", a.ref, "-l", bed_path, a.bam]
        q0 = run_mpileup(common + ["-q", "0"]) if known else {}
        q1 = run_mpileup(common + ["-q", str(a.min_mapq_pass)]) if known else {}
    finally:
        os.unlink(bed_path)
    rows = assemble_target_rows(a.sample, targets, fai_lengths, q0, q1)
    with open(a.out, "w") as fh:
        fh.write(TARGET_HEADER)
        for row in rows:
            fh.write("\t".join(row) + "\n")
    print(f"[{a.sample}] {len(rows)} target row(s) -> {a.out}")
    for row in rows:
        print("  " + "\t".join(row))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("counts", help="genome-wide candidate allele counts")
    c.add_argument("--sample", required=True)
    c.add_argument("--bam", required=True)
    c.add_argument("--ref", required=True)
    c.add_argument("--out", required=True)
    c.add_argument("--max-depth", type=int, default=100000)
    c.add_argument("--min-bq", type=int, default=13)
    c.add_argument("--min-alt-reads", type=int, default=2)
    c.add_argument("--min-alt-af", type=float, default=0.01)
    c.set_defaults(func=cmd_counts)

    t = sub.add_parser("targets", help="targeted-sites report")
    t.add_argument("--sample", required=True)
    t.add_argument("--bam", required=True)
    t.add_argument("--ref", required=True)
    t.add_argument("--fai", required=True)
    t.add_argument("--targets", required=True, help="JSON list or path to JSON")
    t.add_argument("--out", required=True)
    t.add_argument("--max-depth", type=int, default=100000)
    t.add_argument("--min-bq", type=int, default=13)
    t.add_argument("--min-mapq-pass", type=int, default=1)
    t.set_defaults(func=cmd_targets)

    a = ap.parse_args(argv)
    a.func(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
