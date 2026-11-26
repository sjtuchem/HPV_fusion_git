#!/usr/bin/env python3
"""Validate the batch TSV and emit a normalized six-column jobs file."""

from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path
from typing import Dict, Iterable, List


REQUIRED = ("sample", "merged_fastq", "raw_r1_a", "raw_r2_a")
OPTIONAL = ("raw_r1_b", "raw_r2_b")
SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9._@+-]+$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-list", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    parser.add_argument("--no-check-files", action="store_true")
    return parser.parse_args()


def meaningful_lines(path: Path) -> Iterable[str]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        for raw in handle:
            if not raw.strip() or raw.lstrip().startswith("#"):
                continue
            yield raw


def normalize_optional(value: str) -> str:
    value = (value or "").strip()
    return "." if value in ("", ".", "NA", "na", "-") else value


def validate_path(value: str, label: str, sample: str, check_files: bool) -> str:
    if not value or value == ".":
        raise ValueError("sample %s: %s is required" % (sample, label))
    path = Path(value)
    if not path.is_absolute():
        raise ValueError("sample %s: %s must be an absolute path: %s" % (sample, label, value))
    if check_files and (not path.is_file() or path.stat().st_size == 0):
        raise FileNotFoundError("sample %s: missing or empty %s: %s" % (sample, label, value))
    return str(path)


def main() -> None:
    args = parse_args()
    source = Path(args.sample_list)
    if not source.is_file() or source.stat().st_size == 0:
        raise FileNotFoundError("sample.list is missing or empty: %s" % source)
    reader = csv.DictReader(meaningful_lines(source), delimiter="\t")
    if reader.fieldnames is None:
        raise ValueError("sample.list has no header")
    fields = [item.strip() for item in reader.fieldnames]
    reader.fieldnames = fields
    missing = [item for item in REQUIRED if item not in fields]
    unexpected = [item for item in fields if item not in REQUIRED + OPTIONAL]
    if missing or unexpected:
        raise ValueError(
            "sample.list columns must be %s (optional: %s); missing=%s unexpected=%s" %
            (", ".join(REQUIRED), ", ".join(OPTIONAL), missing, unexpected)
        )

    rows: List[Dict[str, str]] = []
    seen = set()
    for line_number, row in enumerate(reader, start=2):
        if None in row:
            raise ValueError("line %d: too many tab-separated columns" % line_number)
        sample = (row.get("sample") or "").strip()
        if not sample or sample in (".", "..") or not SAFE_SAMPLE.fullmatch(sample):
            raise ValueError(
                "line %d: unsafe sample name %r; allowed characters are letters, numbers, . _ @ + -" %
                (line_number, sample)
            )
        if sample in seen:
            raise ValueError("line %d: duplicate sample name: %s" % (line_number, sample))
        seen.add(sample)
        merged = validate_path(
            (row.get("merged_fastq") or "").strip(), "merged_fastq", sample,
            not args.no_check_files,
        )
        r1a = validate_path(
            (row.get("raw_r1_a") or "").strip(), "raw_r1_a", sample,
            not args.no_check_files,
        )
        r2a = validate_path(
            (row.get("raw_r2_a") or "").strip(), "raw_r2_a", sample,
            not args.no_check_files,
        )
        r1b = normalize_optional(row.get("raw_r1_b", "."))
        r2b = normalize_optional(row.get("raw_r2_b", "."))
        if (r1b == ".") != (r2b == "."):
            raise ValueError("line %d: raw_r1_b and raw_r2_b must both be paths or both be ." % line_number)
        if r1b != ".":
            r1b = validate_path(r1b, "raw_r1_b", sample, not args.no_check_files)
            r2b = validate_path(r2b, "raw_r2_b", sample, not args.no_check_files)
        paths = [merged, r1a, r2a] + ([] if r1b == "." else [r1b, r2b])
        if len(set(paths)) != len(paths):
            raise ValueError("line %d: a FASTQ path is reused within sample %s" % (line_number, sample))
        rows.append({
            "sample": sample, "merged_fastq": merged,
            "raw_r1_a": r1a, "raw_r2_a": r2a,
            "raw_r1_b": r1b, "raw_r2_b": r2b,
        })
    if not rows:
        raise ValueError("sample.list contains no active sample rows")

    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_name(destination.name + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=REQUIRED + OPTIONAL, delimiter="\t",
            lineterminator="\n",
        )
        for row in rows:
            writer.writerow(row)
    temporary.replace(destination)

    report = {
        "status": "PASS", "sample_count": len(rows),
        "samples": [row["sample"] for row in rows],
        "source": str(source.resolve()), "normalized_jobs": str(destination.resolve()),
        "file_existence_checked": not args.no_check_files,
    }
    report_path = Path(args.report)
    report_tmp = report_path.with_name(report_path.name + ".tmp")
    report_tmp.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report_tmp.replace(report_path)
    print("[OK] validated %d unique sample rows" % len(rows), file=sys.stderr)


if __name__ == "__main__":
    main()
