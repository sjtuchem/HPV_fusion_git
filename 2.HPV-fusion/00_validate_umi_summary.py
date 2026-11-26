#!/usr/bin/env python3
"""Fail fast when a targeted-library run has no usable UMI evidence."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--summary", required=True)
    parser.add_argument(
        "--required-tube", action="append", default=[], choices=("A", "B")
    )
    return parser.parse_args()


def integer_value(summary: Dict[str, object], key: str, errors: List[str]) -> int:
    value = summary.get(key, 0)
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        errors.append("%s is not numeric" % key)
        return 0
    integer = int(value)
    if integer < 0 or integer != value:
        errors.append("%s is not a non-negative integer" % key)
        return 0
    return integer


def main() -> None:
    args = parse_args()
    summary_path = Path(args.summary)
    try:
        with summary_path.open(encoding="utf-8") as handle:
            summary = json.load(handle)
    except (OSError, ValueError) as error:
        print("[ERROR] Cannot read UMI summary %s: %s" % (summary_path, error), file=sys.stderr)
        raise SystemExit(1)
    if not isinstance(summary, dict):
        print("[ERROR] UMI summary root must be a JSON object: %s" % summary_path, file=sys.stderr)
        raise SystemExit(1)

    errors: List[str] = []
    tube_report: Dict[str, Dict[str, object]] = {}
    required_tubes = list(dict.fromkeys(args.required_tube))
    for tube in required_tubes:
        raw_pairs = integer_value(summary, "raw_pairs_%s" % tube, errors)
        eligible = integer_value(summary, "raw_umi_eligible_%s" % tube, errors)
        if raw_pairs == 0:
            errors.append("raw tube %s contains zero read pairs" % tube)
        if eligible == 0:
            errors.append("raw tube %s contains zero UMI-eligible read pairs" % tube)
        tube_report[tube] = {
            "raw_pairs": raw_pairs,
            "umi_eligible_pairs": eligible,
            "umi_eligible_rate": eligible / float(raw_pairs) if raw_pairs else 0.0,
        }

    merged_reads = integer_value(summary, "merged_reads", errors)
    merged_eligible = integer_value(summary, "merged_umi_eligible_reads", errors)
    if merged_reads == 0:
        errors.append("merged FASTQ contains zero reads")
    if merged_eligible == 0:
        errors.append("merged FASTQ contains zero reads associated with an eligible UMI")

    report = {
        "status": "FAIL" if errors else "PASS",
        "sample": summary.get("sample", "."),
        "required_tubes": required_tubes,
        "tubes": tube_report,
        "merged_reads": merged_reads,
        "merged_umi_eligible_reads": merged_eligible,
        "merged_umi_eligible_rate": (
            merged_eligible / float(merged_reads) if merged_reads else 0.0
        ),
        "errors": errors,
    }
    print(json.dumps(report, ensure_ascii=False, sort_keys=True))
    if errors:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
