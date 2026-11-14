#!/usr/bin/env python3
"""Schema and router tests for the typed sample sheet.

Runs without any sequencing data or pipeline execution:
    python tests/test_sample_sheet.py
"""

import os
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO / "workflow" / "scripts"))

from sample_sheet import load_samples, require_implemented_types, samples_by_type

SCHEMA = str(REPO / "workflow" / "schemas" / "samples.schema.yaml")
CI_SHEET = str(REPO / "config" / "ci_samples.tsv")

HEADER = "sample\ttype\tploidy\tpool_size\treference\tgroup\tfastq_r1\tfastq_r2\n"


def write_sheet(text):
    fd, path = tempfile.mkstemp(suffix=".tsv")
    with os.fdopen(fd, "w") as fh:
        fh.write(text)
    return path


def expect_fail(fn, needle, label):
    try:
        fn()
    except Exception as e:  # noqa: BLE001 - any loud failure is the contract
        msg = str(e)
        assert needle in msg, f"{label}: wanted {needle!r} in error, got:\n  {msg}"
        print(f"ok (fails loudly): {label}")
        return
    raise AssertionError(f"{label}: expected a failure, but it passed")


def main():
    # 1. The CI sheet loads and is typed
    samples = load_samples(CI_SHEET, "", SCHEMA)
    assert len(samples) == 2, "CI sheet should have 2 samples"
    assert all(s["type"] == "individual_wgs" for s in samples)
    assert all(s["ploidy"] == 2 for s in samples)
    print("ok: CI sheet loads as typed individual_wgs with ploidy default 2")

    # 2. Valid individual row (user's example shape)
    path = write_sheet(
        HEADER + "SRR8305757\tindividual_wgs\t2\t\tAaegL5\t\tfq/1_1.fq\tfq/1_2.fq\n"
    )
    s = load_samples(path, "", SCHEMA)[0]
    assert s["type"] == "individual_wgs" and s["ploidy"] == 2 and s["reference"] == "AaegL5"
    os.unlink(path)
    print("ok: valid individual_wgs row")

    # 3. Valid pool row (pool_size required, ploidy must stay empty)
    path = write_sheet(
        HEADER + "14566518\tpool_seq\t\t50\tAaegL5\tkdr\tfq/2_1.fq\tfq/2_2.fq\n"
    )
    s = load_samples(path, "", SCHEMA)[0]
    assert s["type"] == "pool_seq" and s["pool_size"] == 50 and s["ploidy"] is None
    os.unlink(path)
    print("ok: valid pool_seq row with pool_size=50")

    # 4. Minimal sheet: only required columns
    path = "sample\ttype\tfastq_r1\tfastq_r2\nmin\tindividual_wgs\tfq/a_1.fq\ta_2.fq\n"
    path = write_sheet(path)
    s = load_samples(path, "", SCHEMA)[0]
    assert s["ploidy"] == 2  # default applied
    os.unlink(path)
    print("ok: required columns only, ploidy defaults to 2")

    # 5. Missing type column
    path = write_sheet(
        "sample\tfastq_r1\tfastq_r2\nx\tfq/x_1.fq\tx_2.fq\n"
    )
    expect_fail(lambda: load_samples(path, "", SCHEMA), "type", "missing type column")
    os.unlink(path)

    # 6. Unknown type
    path = write_sheet(
        HEADER + "x\tpool\t\t25\t\t\tfq/x_1.fq\tx_2.fq\n"
    )
    expect_fail(lambda: load_samples(path, "", SCHEMA), "not one of", "unknown type enum")
    os.unlink(path)

    # 7. pool_seq without pool_size
    path = write_sheet(HEADER + "x\tpool_seq\t\t\t\t\tfq/x_1.fq\tx_2.fq\n")
    expect_fail(
        lambda: load_samples(path, "", SCHEMA), "pool_size", "pool_seq without pool_size"
    )
    os.unlink(path)

    # 8. pool_seq with ploidy set (pools have no diploid genotype)
    path = write_sheet(HEADER + "x\tpool_seq\t2\t25\t\t\tfq/x_1.fq\tx_2.fq\n")
    expect_fail(
        lambda: load_samples(path, "", SCHEMA), "is not one of", "pool_seq with ploidy"
    )
    os.unlink(path)

    # 9. individual_wgs with pool_size set
    path = write_sheet(HEADER + "x\tindividual_wgs\t2\t25\t\t\tfq/x_1.fq\tx_2.fq\n")
    expect_fail(
        lambda: load_samples(path, "", SCHEMA), "pool_seq", "individual with pool_size"
    )
    os.unlink(path)

    # 10. Duplicate sample names
    path = write_sheet(
        HEADER + "dup\tindividual_wgs\t2\t\t\t\tfq/a_1.fq\ta_2.fq\n"
        + "dup\tindividual_wgs\t2\t\t\t\tfq/b_1.fq\tb_2.fq\n"
    )
    expect_fail(lambda: load_samples(path, "", SCHEMA), "Duplicate", "duplicate sample names")
    os.unlink(path)

    # 11. Non-integer pool_size
    path = write_sheet(HEADER + "x\tpool_seq\t\tbig\t\t\tfq/x_1.fq\tx_2.fq\n")
    expect_fail(lambda: load_samples(path, "", SCHEMA), "integer", "non-integer pool_size")
    os.unlink(path)

    # 12. Router: unbuilt types fail loudly; implemented types pass
    mixed = load_samples(
        write_sheet(
            HEADER
            + "i1\tindividual_wgs\t2\t\t\t\tfq/i1_1.fq\ti1_2.fq\n"
            + "p1\tpool_seq\t\t25\t\tkdr\tfq/p1_1.fq\tp1_2.fq\n"
        ),
        "",
        SCHEMA,
    )
    by_type = samples_by_type(mixed)
    assert len(by_type["individual_wgs"]) == 1 and len(by_type["pool_seq"]) == 1
    try:
        require_implemented_types(mixed, {"individual_wgs"})
    except Exception as e:  # noqa: BLE001
        assert "pool_seq" in str(e) and "not implemented" in str(e), str(e)
        print("ok (fails loudly): unbuilt pool_seq type rejected by router")
    else:
        raise AssertionError("router accepted an unbuilt type")
    require_implemented_types(mixed[:1], {"individual_wgs"})  # individual passes
    print("ok: router accepts implemented types")

    # 13. fastq_dir discovery defaults to individual_wgs (with warning)
    tmp = tempfile.mkdtemp()
    for f in ("disc_R1_001.fastq.gz", "disc_R2_001.fastq.gz"):
        open(os.path.join(tmp, f), "wb").close()
    s = load_samples("", tmp, SCHEMA)[0]
    assert s["type"] == "individual_wgs" and s["ploidy"] == 2
    print("ok: fastq_dir discovery defaults type=individual_wgs")

    print("\nAll sample-sheet schema/router tests passed.")


if __name__ == "__main__":
    main()
