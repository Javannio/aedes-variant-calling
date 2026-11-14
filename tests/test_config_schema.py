#!/usr/bin/env python3
"""Config schema tests: shipped configs validate, drift fails loudly.

    python tests/test_config_schema.py
"""

import copy
import sys
from pathlib import Path

import yaml
from jsonschema import Draft7Validator, ValidationError

REPO = Path(__file__).resolve().parent.parent
SCHEMA_PATH = REPO / "workflow" / "schemas" / "config.schema.yaml"
BASE = REPO / "config" / "config.yaml"
CI = REPO / "config" / "config.ci.yaml"
LOCAL = REPO / "config" / "config.local.yaml"  # gitignored, may be absent


def deep_merge(base, over):
    """Snakemake-style recursive dict merge (over wins)."""
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load(path):
    with open(path) as fh:
        return yaml.safe_load(fh)


def expect_invalid(validator, config, needle, label):
    errors = list(validator.iter_errors(config))
    assert errors, f"{label}: expected a validation error, but config passed"
    msgs = "\n".join(str(e) for e in errors)
    assert needle in msgs, f"{label}: wanted {needle!r} in:\n  {msgs}"
    print(f"ok (fails loudly): {label}")


def main():
    schema = load(SCHEMA_PATH)
    validator = Draft7Validator(schema)

    # 1. Shipped base config validates
    base = load(BASE)
    validator.validate(base)
    print("ok: config/config.yaml validates against config.schema.yaml")

    # 2. CI overlay merged over base validates (snakemake merge order)
    ci = load(CI)
    merged_ci = deep_merge(base, ci)
    validator.validate(merged_ci)
    print("ok: config.ci.yaml merged over config.yaml validates")

    # 3. User local overlay merged in also validates (when present)
    if LOCAL.exists():
        merged_all = deep_merge(merged_ci, load(LOCAL))
        validator.validate(merged_all)
        print("ok: config.local.yaml merged in validates")
    else:
        print("skip: no config/config.local.yaml present")

    # 4. Blocks the workflow reads are actually present and typed
    qc = merged_ci.get("qc", {})
    for key in ("enabled", "min_mapping_rate", "min_coverage_depth",
                "min_genome_fraction_covered", "max_duplicate_rate"):
        assert key in qc, f"qc.{key} missing from base config"
    pool = merged_ci.get("pool_seq", {})
    for key in ("max_depth", "min_base_quality", "min_mapq_pass",
                "min_alt_reads", "min_alt_af"):
        assert key in pool, f"pool_seq.{key} missing from base config"
    filters = merged_ci["tools"]["bcftools"]
    for caller in ("gatk", "deepvariant"):
        for key in ("snp_filter", "indel_filter", "genotype_filter"):
            assert key in filters.get(caller, {}), \
                f"tools.bcftools.{caller}.{key} missing"
    # CI overlay replaces reporting.targets, so the kdr sites live in base
    targets = merged_ci["reporting"]["targets"]
    assert targets, "reporting.targets must not be empty"
    names = [t["name"] for t in targets]
    assert "planted_snp" in names, "CI target planted_snp missing"
    kdr = {t["name"]: t for t in base["reporting"]["targets"]}
    for site in ("V410L", "V1016I", "F1534C"):
        assert site in kdr, f"kdr target {site} missing from config.yaml"
        assert isinstance(kdr[site]["pos"], int) and kdr[site]["pos"] > 0
    print("ok: qc/pool_seq/filter-set/targets blocks present and typed")

    # 5. Out-of-range values fail (schema guards, not just documentation)
    bad = deep_merge(base, {"qc": {"min_mapping_rate": 1.5}})
    expect_invalid(validator, bad, "1.5", "qc.min_mapping_rate > 1")
    bad = deep_merge(base, {"pool_seq": {"min_alt_af": 2}})
    expect_invalid(validator, bad, "2", "pool_seq.min_alt_af > 1")
    bad = deep_merge(base, {"tools": {"threads": 0}})
    expect_invalid(validator, bad, "0", "tools.threads < 1")
    bad = deep_merge(base, {"reporting": {"targets": [{"name": "x", "chrom": "c",
                                                       "pos": 1, "ref": "A"}]}})
    expect_invalid(validator, bad, "alt", "target missing alt")
    bad = deep_merge(base, {"reporting": {"targets": [{"name": "x", "chrom": "c",
                                                       "pos": 0, "ref": "A",
                                                       "alt": "G"}]}})
    expect_invalid(validator, bad, "0", "target pos < 1")

    # 6. Unknown keys are allowed (documented: local configs carry extras)
    extra = deep_merge(base, {"my_local_tool": {"anything": True}})
    validator.validate(extra)
    print("ok: unknown keys allowed (local configs may carry extras)")

    # 7. Filter sets are plain bcftools expressions (non-empty strings)
    for caller, fs in filters.items():
        if caller in ("output_type",):
            continue
        if not isinstance(fs, dict) or "snp_filter" not in fs:
            continue
        for key in ("snp_filter", "indel_filter", "genotype_filter"):
            assert isinstance(fs[key], str) and fs[key].strip(), \
                f"tools.bcftools.{caller}.{key} must be a non-empty expression"
    print("ok: per-caller filter sets are non-empty strings")

    print("\nAll config-schema tests passed.")


if __name__ == "__main__":
    main()
