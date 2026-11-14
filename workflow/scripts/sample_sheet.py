"""Typed sample sheet loading and validation.

Every sample declares what it is (individual_wgs, pool_seq, mito, long_read).
The parsed sheet is validated against workflow/schemas/samples.schema.yaml
before any rule runs, and non-implemented types fail loudly at parse time.
"""

import os
from pathlib import Path

import pandas as pd
import yaml
from jsonschema import ValidationError, validate as jsonschema_validate
from snakemake.exceptions import WorkflowError

REQUIRED_COLUMNS = ["sample", "type", "fastq_r1", "fastq_r2"]
OPTIONAL_COLUMNS = ["ploidy", "pool_size", "reference", "group"]
SAMPLE_TYPES = ("individual_wgs", "pool_seq", "mito", "long_read")
DEFAULT_PLOIDY = {"individual_wgs": 2, "mito": 1, "long_read": 2}


def _clean(value):
    """Treat pandas NaN / empty / 'nan' cells as None."""
    if value is None:
        return None
    if isinstance(value, float) and pd.isna(value):
        return None
    text = str(value).strip()
    if not text or text.lower() in ("nan", "none", "-"):
        return None
    return text


def _as_int(value, column, sample):
    if value is None:
        return None
    try:
        return int(float(value))
    except (TypeError, ValueError):
        raise ValueError(
            f"Sample sheet: column '{column}' must be an integer "
            f"(got {value!r} for sample '{sample}')"
        )


def _row_to_sample(row, sheet_path):
    sample = _clean(row.get("sample"))
    if sample is None:
        raise ValueError(f"Sample sheet '{sheet_path}': a row is missing the 'sample' name")
    out = {
        "sample": sample,
        "type": _clean(row.get("type")),
        "fastq_r1": _clean(row.get("fastq_r1")),
        "fastq_r2": _clean(row.get("fastq_r2")),
        "ploidy": _as_int(_clean(row.get("ploidy")), "ploidy", sample),
        "pool_size": _as_int(_clean(row.get("pool_size")), "pool_size", sample),
        "reference": _clean(row.get("reference")),
        "group": _clean(row.get("group")),
    }
    # Keys used by the rest of the Snakefile
    out["r1"] = out["fastq_r1"]
    out["r2"] = out["fastq_r2"]
    return out


def _validate_samples(samples, schema_path, sheet_label):
    with open(schema_path) as fh:
        schema = yaml.safe_load(fh)
    try:
        jsonschema_validate(instance={"samples": samples}, schema=schema)
    except ValidationError as e:
        where = "/".join(str(p) for p in e.absolute_path) or "<sheet>"
        raise ValueError(
            f"Invalid sample sheet '{sheet_label}': {where}: {e.message}\n"
            f"  Required columns: {', '.join(REQUIRED_COLUMNS)}; "
            f"optional: {', '.join(OPTIONAL_COLUMNS)}. "
            f"Valid types: {', '.join(SAMPLE_TYPES)}."
        ) from e
    names = [s["sample"] for s in samples]
    dups = sorted({n for n in names if names.count(n) > 1})
    if dups:
        raise ValueError(
            f"Duplicate sample name(s) in '{sheet_label}': {', '.join(dups)}"
        )
    for s in samples:
        if s["type"] in DEFAULT_PLOIDY and s["ploidy"] is None:
            s["ploidy"] = DEFAULT_PLOIDY[s["type"]]


def load_samples(sheet_path, fastq_dir, schema_path):
    """Load and validate the sample sheet (or discover FASTQs in fastq_dir).

    Returns a list of dicts: sample, type, ploidy, pool_size, reference,
    group, fastq_r1/r1, fastq_r2/r2. Raises ValueError on any schema
    violation — nothing runs before the sheet is valid.
    """
    samples = []
    if sheet_path and os.path.exists(sheet_path):
        df = pd.read_csv(sheet_path, sep="\t")
        missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(
                f"Sample sheet '{sheet_path}' is missing required column(s): "
                f"{', '.join(missing)}. Required: {', '.join(REQUIRED_COLUMNS)}; "
                f"optional: {', '.join(OPTIONAL_COLUMNS)}"
            )
        if df.empty:
            raise ValueError(f"Sample sheet '{sheet_path}' has no sample rows")
        for _, row in df.iterrows():
            samples.append(_row_to_sample(row.to_dict(), sheet_path))
        _validate_samples(samples, schema_path, sheet_path)
    elif fastq_dir and os.path.exists(fastq_dir):
        fastq_files = sorted(Path(fastq_dir).glob("*_R1_*.fastq.gz"))
        for r1 in fastq_files:
            sample = r1.name.split("_R1_")[0]
            r2 = r1.parent / r1.name.replace("_R1_", "_R2_")
            if r2.exists():
                samples.append({
                    "sample": sample,
                    "type": "individual_wgs",
                    "fastq_r1": str(r1),
                    "fastq_r2": str(r2),
                    "ploidy": 2,
                    "pool_size": None,
                    "reference": None,
                    "group": None,
                    "r1": str(r1),
                    "r2": str(r2),
                })
            else:
                print(f"Warning: No matching R2 for {r1}")
        if not samples:
            raise ValueError(f"No samples found in fastq_dir '{fastq_dir}'")
        _validate_samples(samples, schema_path, fastq_dir)
        print(
            "Warning: fastq_dir auto-discovery used — defaulting type=individual_wgs "
            "for all samples (use a typed sample sheet to declare pool_seq/mito/long_read)"
        )
    else:
        raise ValueError("No samples found. Provide samples.tsv or fastq_dir in config.")
    return samples


def require_implemented_types(samples, implemented):
    """Fail loudly (parse time) for sample types whose branch is not built yet."""
    by_type = {}
    for s in samples:
        by_type.setdefault(s["type"], []).append(s["sample"])
    for stype in sorted(by_type):
        if stype not in implemented:
            names = ", ".join(by_type[stype])
            raise WorkflowError(
                f"Sample type '{stype}' is not implemented yet "
                f"(sample(s): {names}). Implemented types: "
                f"{', '.join(sorted(implemented))}. Planned branches: "
                f"pool_seq allele counting + targeted-sites report, mito haploid "
                f"calling + consensus FASTA, long_read alignment."
            )


def samples_by_type(samples):
    """Route: {type: [sample, ...]} for every type present in the sheet."""
    out = {t: [] for t in SAMPLE_TYPES}
    for s in samples:
        out[s["type"]].append(s)
    return out
