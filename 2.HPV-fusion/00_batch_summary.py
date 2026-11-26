#!/usr/bin/env python3
"""Summarize completion state and final counts for one batch invocation."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


FIELDS = (
    "sample", "status", "final_exact_events", "final_loci", "final_evidence_rows",
    "final_unique_sequences", "main_events", "candidate_events", "review_events",
    "runtime_seconds", "job_exitval", "final_directory",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", required=True)
    parser.add_argument("--result-dir", required=True)
    parser.add_argument("--joblog", default="")
    parser.add_argument("--output-tsv", required=True)
    parser.add_argument("--output-json", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with open(args.jobs, "r", encoding="utf-8", newline="") as handle:
        jobs = list(csv.DictReader(handle, fieldnames=(
            "sample", "merged_fastq", "raw_r1_a", "raw_r2_a", "raw_r1_b", "raw_r2_b"
        ), delimiter="\t"))
    exit_values = {}
    if args.joblog and Path(args.joblog).is_file():
        with open(args.joblog, "r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                exit_values[int(row["Seq"])] = int(row["Exitval"])
    rows = []
    for sequence_number, job in enumerate(jobs, start=1):
        sample = job["sample"]
        root = Path(args.result_dir) / sample
        final = root / "final"
        complete = root / "logs" / (sample + ".complete")
        summary_path = final / (sample + ".fusion_summary.final.json")
        validation_path = root / "logs" / "11_validate_results.log"
        status = "INCOMPLETE_OR_FAILED"
        summary = {}
        validation = {}
        if summary_path.is_file() and summary_path.stat().st_size:
            summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if validation_path.is_file() and validation_path.stat().st_size:
            try:
                validation = json.loads(validation_path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                validation = {}
        job_exitval = exit_values.get(sequence_number, "")
        if job_exitval != "" and job_exitval != 0:
            status = "FAILED_THIS_BATCH"
        elif complete.is_file() and validation.get("status") == "PASS":
            status = "COMPLETE_PASS"
        runtime = ""
        timing = root / "logs" / (sample + ".runtime_metrics.tsv")
        if timing.is_file():
            for line in timing.read_text(encoding="utf-8").splitlines():
                parts = line.split("\t")
                if len(parts) == 3 and parts[0] == "TOTAL":
                    runtime = parts[2]
        rows.append({
            "sample": sample, "status": status,
            "final_exact_events": summary.get("final_exact_events", ""),
            "final_loci": summary.get("final_loci", ""),
            "final_evidence_rows": summary.get("final_evidence_rows", ""),
            "final_unique_sequences": summary.get("final_unique_sequences", ""),
            "main_events": validation.get("main_events", ""),
            "candidate_events": validation.get("candidate_events", ""),
            "review_events": validation.get("review_events", ""),
            "runtime_seconds": runtime, "job_exitval": job_exitval,
            "final_directory": str(final),
        })
    destination = Path(args.output_tsv)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    counts = {status: sum(row["status"] == status for row in rows) for status in (
        "COMPLETE_PASS", "FAILED_THIS_BATCH", "INCOMPLETE_OR_FAILED"
    )}
    nonpass = len(rows) - counts["COMPLETE_PASS"]
    report = {"sample_count": len(rows), "status_counts": counts, "samples": rows}
    Path(args.output_json).write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print("[BATCH SUMMARY] samples=%d complete_pass=%d nonpass=%d" % (
        len(rows), counts["COMPLETE_PASS"], nonpass
    ))
    if nonpass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
