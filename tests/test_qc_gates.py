#!/usr/bin/env python3
"""Unit tests for the QC gate evaluator (no pipeline run needed).

    python tests/test_qc_gates.py
"""

import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

from qc_gates import evaluate_gates, main, parse_coverage, parse_flagstat

FLAGSTAT_OK = """\
1000 + 0 in total (QC-passed reads + QC-failed reads)
1000 + 0 primary
0 + 0 secondary
0 + 0 supplementary
0 + 0 duplicates
950 + 0 primary duplicates
980 + 0 mapped (98.00% : N/A)
980 + 0 primary mapped (98.00% : N/A)
900 + 0 properly paired (90.00% : N/A)
"""

FLAGSTAT_LOW_MAP = FLAGSTAT_OK.replace("980 + 0 mapped (98.00% : N/A)", "500 + 0 mapped (50.00% : N/A)")

COVERAGE_OK = """\
#rname\tstartpos\tendpos\tnumreads\tcovbases\tcoverage\tmeandepth\tmeanbaseq\tmeanmapq
chr1\t1\t1000\t500\t950\t95.00\t20.00\t30.0\t60.0
chr2\t1\t1000\t400\t900\t90.00\t15.00\t30.0\t60.0
"""

GATES = {
    "enabled": True,
    "min_mapping_rate": 0.80,
    "min_coverage_depth": 10,
    "min_genome_fraction_covered": 0.80,
    "max_duplicate_rate": 0.50,
}


def status(rows, metric):
    return {m: s for m, _, _, _, s in rows}[metric]


def main_tests():
    # 1. flagstat parsing (modern format, duplicates before mapped lines)
    m = parse_flagstat(FLAGSTAT_OK)
    assert m["total_reads"] == 1000
    assert m["mapped_reads"] == 980
    assert m["duplicate_reads"] == 0
    assert abs(m["mapping_rate"] - 0.98) < 1e-9
    assert m["duplicate_rate"] == 0.0
    print("ok: flagstat parsing")

    # 2. coverage parsing: mean depth = len-weighted (20/15 -> 17.5), breadth 0.925
    c = parse_coverage(COVERAGE_OK)
    assert abs(c["mean_depth"] - 17.5) < 1e-9, c
    assert abs(c["genome_fraction_covered"] - 0.925) < 1e-9, c
    print("ok: coverage parsing")

    # 3. all gates pass
    rows = evaluate_gates({**m, **c}, GATES)
    assert all(s == "PASS" for *_, s in rows), rows
    print("ok: all gates pass on good metrics")

    # 4. low mapping rate fails only that gate
    low = dict(m, mapping_rate=0.50)
    rows = evaluate_gates({**low, **c}, GATES)
    assert status(rows, "mapping_rate") == "FAIL"
    assert status(rows, "mean_depth") == "PASS"
    print("ok: low mapping rate fails mapping gate only")

    # 5. high duplicate rate fails
    dup = dict(m, duplicate_rate=0.60)
    rows = evaluate_gates({**dup, **c}, GATES)
    assert status(rows, "duplicate_rate") == "FAIL"
    print("ok: high duplicate rate fails")

    # 6. shallow / sparse coverage fails
    shallow = dict(c, mean_depth=2.0, genome_fraction_covered=0.30)
    rows = evaluate_gates({**m, **shallow}, GATES)
    assert status(rows, "mean_depth") == "FAIL"
    assert status(rows, "genome_fraction_covered") == "FAIL"
    print("ok: shallow + sparse coverage fail")

    # 7. disabled gates never fail
    off = dict(GATES, enabled=False)
    rows = evaluate_gates({**low, **shallow}, off)
    assert all(s == "SKIP" for *_, s in rows), rows
    print("ok: disabled gates write SKIP, never FAIL")

    # 8. end-to-end via main(): pass writes report + exit 0, fail exits 1
    tmp = Path(tempfile.mkdtemp())
    fs, cov = tmp / "f.flagstat", tmp / "c.coverage"
    fs.write_text(FLAGSTAT_OK)
    cov.write_text(COVERAGE_OK)
    out = tmp / "gates.tsv"
    rc = main(["--sample", "s1", "--flagstat", str(fs), "--coverage", str(cov),
               "--gates", '{"min_mapping_rate": 0.80}', "--out", str(out)])
    assert rc == 0 and out.exists()
    assert "PASS" in out.read_text()
    fs.write_text(FLAGSTAT_LOW_MAP)
    out2 = tmp / "gates_fail.tsv"
    rc = main(["--sample", "s1", "--flagstat", str(fs), "--coverage", str(cov),
               "--gates", '{"min_mapping_rate": 0.80}', "--out", str(out2)])
    assert rc == 1
    text = out2.read_text()
    assert "FAIL" in text and "mapping_rate" in text
    print("ok: main() exits 0 on pass, 1 on fail, report always written")

    print("\nAll QC gate tests passed.")


if __name__ == "__main__":
    main_tests()
