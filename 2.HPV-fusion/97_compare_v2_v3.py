#!/usr/bin/env python3
"""Fail-fast biological regression comparison between completed v2 and v3 results."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--v2-final-dir", required=True)
    parser.add_argument("--v3-final-dir", required=True)
    parser.add_argument("--output-json", default="")
    return parser.parse_args()


def read_tsv(path: Path) -> List[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "rt" if str(path).endswith(".gz") else "r"
    with opener(path, mode, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def indexed_events(root: Path, sample: str) -> Dict[str, dict]:
    path = root / (sample + ".fusion_events.final.tsv")
    rows = read_tsv(path)
    result = {row["event_id"]: row for row in rows}
    if len(result) != len(rows):
        raise AssertionError("Duplicate event ID in %s" % path)
    return result


def main() -> None:
    args = parse_args()
    v2_root = Path(args.v2_final_dir)
    v3_root = Path(args.v3_final_dir)
    v2 = indexed_events(v2_root, args.sample)
    v3 = indexed_events(v3_root, args.sample)
    missing = sorted(set(v2) - set(v3))
    added = sorted(set(v3) - set(v2))
    stable_fields = [
        "hpv_type", "hpv_ref", "hpv_ref_length", "hpv_breakpoint_1based",
        "hpv_strand", "hpv_breakend_side", "human_chr", "human_breakpoint_1based",
        "human_breakpoint_start0", "human_breakpoint_end0", "human_strand",
        "human_breakend_side", "query_orientation", "result_tier",
    ] + [
        "%s_%dbp" % (measure, threshold)
        for threshold in (30, 40, 50)
        for measure in (
            "raw_reads", "umi_eligible_reads", "unique_umi_families"
        )
    ]
    differences = []
    for event_id in sorted(set(v2) & set(v3)):
        for field in stable_fields:
            if v2[event_id].get(field) != v3[event_id].get(field):
                differences.append({
                    "event_id": event_id, "field": field,
                    "v2": v2[event_id].get(field), "v3": v3[event_id].get(field),
                })
    v2_sequences = read_tsv(v2_root / (args.sample + ".fusion_unique_sequences.final.tsv.gz"))
    v3_sequences = read_tsv(v3_root / (args.sample + ".fusion_unique_sequences.final.tsv.gz"))
    v2_evidence = read_tsv(v2_root / (args.sample + ".fusion_evidence.final.tsv.gz"))
    v3_evidence = read_tsv(v3_root / (args.sample + ".fusion_evidence.final.tsv.gz"))
    report = {
        "status": "PASS" if not missing and not added and not differences else "FAIL",
        "sample": args.sample,
        "v2_exact_events": len(v2), "v3_exact_events": len(v3),
        "missing_event_ids": missing, "added_event_ids": added,
        "stable_field_difference_count": len(differences),
        "first_stable_field_differences": differences[:20],
        "v2_unique_original_sequences": len(v2_sequences),
        "v3_unique_normalized_sequences": len(v3_sequences),
        "v2_event_sequence_evidence": len(v2_evidence),
        "v3_event_sequence_evidence": len(v3_evidence),
        "normalization_removed_sequence_rows": len(v2_sequences) - len(v3_sequences),
        "normalization_removed_evidence_rows": len(v2_evidence) - len(v3_evidence),
    }
    text = json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True)
    print(text)
    if args.output_json:
        output = Path(args.output_json)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(text + "\n", encoding="utf-8")
    if report["status"] != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
