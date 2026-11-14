#!/usr/bin/env python3
"""Unit tests for the pileup allele-count / targeted-sites report.

    python tests/test_pileup_report.py
"""

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

import pileup_report as pr


def counts(**kw):
    c = {b: 0 for b in ("A", "C", "G", "T")}
    c.update({"ref": 0, "N": 0, "del": 0, "skip": 0})
    c.update(kw)
    return c


def main_tests():
    # 1. base string parsing: ref matches, anchors, indels, deletions, N
    c = pr.parse_pileup_bases("A", ".,.,..")
    assert c["ref"] == 6, c
    c = pr.parse_pileup_bases("A", "^]..$")
    assert c["ref"] == 2, c
    c = pr.parse_pileup_bases("A", ".+2AG..")
    assert c["ref"] == 3 and c["A"] == 0, c  # inserted bases consumed, not counted
    c = pr.parse_pileup_bases("C", ".,*#acgtnN<>")
    assert c["ref"] == 2 and c["del"] == 2 and c["A"] == 1 and c["C"] == 1
    assert c["G"] == 1 and c["T"] == 1 and c["N"] == 2 and c["skip"] == 2, c
    print("ok: mpileup base-string parsing (anchors, indels, dels, N)")

    # 2. candidate emission criteria
    chrom, pos, ref, depth = "chr1", 100, "A", 20
    strong = counts(ref=4, C=16)
    row = pr.emit_counts_row(chrom, pos, ref, depth, strong, 2, 0.01)
    assert row is not None and row[4] == 4 and row[11] == "C" and row[12] == 16
    assert abs(row[13] - 0.8) < 1e-9, row
    weak = counts(ref=18, C=1)
    assert pr.emit_counts_row(chrom, pos, ref, depth, weak, 2, 0.01) is None  # 1 alt read
    rare = counts(ref=998, C=2)
    assert pr.emit_counts_row(chrom, pos, ref, 1000, rare, 2, 0.05) is None  # AF 0.002
    assert pr.emit_counts_row(chrom, pos, ref, 0, counts(), 2, 0.01) is None  # no acgt
    print("ok: candidate emission (min reads / min AF / empty)")

    # 3. target rows: OK with mapq0 audit
    targets = [{"name": "kdr", "chrom": "chr1", "pos": 100, "ref": "A", "alt": "C"}]
    fai = {"chr1": 1000}
    q0 = {("chr1", 100): ("A", 12, counts(ref=4, C=6, N=2))}
    q1 = {("chr1", 100): ("A", 6, counts(ref=2, C=4))}
    rows = pr.assemble_target_rows("pool1", targets, fai, q0, q1)
    assert len(rows) == 1
    r = rows[0]
    assert len(r) == 15, r  # 6 id + 8 metrics + status
    assert r[6] == "12" and r[7] == "10" and r[8] == "4" and r[9] == "6" and r[10] == "0"
    assert abs(float(r[11]) - 0.6) < 1e-9      # af_alt = 6/10
    assert abs(float(r[12]) - 0.5) < 1e-9      # af incl ambiguous = 6/12
    assert abs(float(r[13]) - 0.4) < 1e-9      # mapq0 fraction = 1 - 6/10
    assert r[14] == "OK"
    print("ok: target row OK (af, ambiguous AF, mapq0 fraction)")

    # 4. MAPQ0_ONLY: all reads fail the MAPQ pass (F1534C situation)
    q1_empty = {}
    rows = pr.assemble_target_rows("pool1", targets, fai, q0, q1_empty)
    assert rows[0][14] == "MAPQ0_ONLY" and float(rows[0][13]) == 1.0
    print("ok: MAPQ0_ONLY when no read passes the MAPQ filter")

    # 5. NO_COVERAGE / REF_MISSING / OUT_OF_RANGE / REF_MISMATCH
    rows = pr.assemble_target_rows("s", targets, fai, {}, {})
    assert rows[0][14] == "NO_COVERAGE"
    rows = pr.assemble_target_rows("s", [{"name": "x", "chrom": "nope", "pos": 1,
                                          "ref": "A", "alt": "C"}], fai, q0, q1)
    assert rows[0][14] == "REF_MISSING"
    rows = pr.assemble_target_rows("s", [{"name": "x", "chrom": "chr1", "pos": 2000,
                                          "ref": "A", "alt": "C"}], fai, q0, q1)
    assert rows[0][14] == "OUT_OF_RANGE"
    q_bad = {("chr1", 100): ("G", 10, counts(ref=10))}
    rows = pr.assemble_target_rows("s", targets, fai, q_bad, q1)
    assert rows[0][14] == "REF_MISMATCH"
    print("ok: NO_COVERAGE / REF_MISSING / OUT_OF_RANGE / REF_MISMATCH")

    # 6. cmd_targets end-to-end with a stubbed samtools mpileup
    tmp = Path(tempfile.mkdtemp())
    fai_path = tmp / "ref.fai"
    fai_path.write_text("chr1\t1000\t0\t0\t0\n")
    targets_path = tmp / "targets.json"
    targets_path.write_text('[{"name": "planted", "chrom": "chr1", "pos": 151, "ref": "C", "alt": "G"}]')
    out = tmp / "report.tsv"
    calls = []

    def fake_run(args):
        calls.append(args)
        if "-q" in args and args[args.index("-q") + 1] == "0":
            return {("chr1", 151): ("C", 10, counts(ref=5, G=5))}
        return {("chr1", 151): ("C", 10, counts(ref=5, G=5))}

    orig = pr.run_mpileup
    pr.run_mpileup = fake_run
    try:
        ns = type("A", (), {})()
        ns.sample, ns.bam, ns.ref, ns.fai = "s1", "x.bam", "ref.fa", str(fai_path)
        ns.targets, ns.out = str(targets_path), str(out)
        ns.max_depth, ns.min_bq, ns.min_mapq_pass = 100000, 13, 1
        pr.cmd_targets(ns)
    finally:
        pr.run_mpileup = orig
    text = out.read_text()
    lines = text.strip().split("\n")
    assert lines[0] == pr.TARGET_HEADER.strip()
    assert len(lines) == 2 and "planted" in lines[1] and "OK" in lines[1]
    assert len(calls) == 2  # both MAPQ passes ran
    print("ok: cmd_targets writes the report (stubbed mpileup, 2 passes)")

    print("\nAll pileup report tests passed.")


if __name__ == "__main__":
    main_tests()
