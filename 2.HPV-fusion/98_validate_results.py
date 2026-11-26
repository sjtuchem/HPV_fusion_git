#!/usr/bin/env python3
"""Independent validator for normalized HPV-fusion final files."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Set, Tuple


csv.field_size_limit(1024 * 1024 * 1024)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--final-dir", required=True)
    parser.add_argument("--primary-threshold", type=int, default=30)
    parser.add_argument("--main-min-raw-reads", type=int, default=3)
    parser.add_argument("--main-min-umi-families", type=int, default=2)
    parser.add_argument("--candidate-min-raw-reads", type=int, default=2)
    parser.add_argument("--candidate-min-umi-families", type=int, default=1)
    parser.add_argument("--audit-sequence-metadata", default="")
    parser.add_argument("--old-sequence-csv", default="")
    return parser.parse_args()


def as_int(row: dict, key: str) -> int:
    return int(row.get(key, "0") or 0)


def read_tsv(path: Path) -> List[dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    mode = "rt" if str(path).endswith(".gz") else "r"
    with opener(path, mode, encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def expected_tier(row: dict, threshold: int, args: argparse.Namespace) -> str:
    raw = as_int(row, "raw_reads_%dbp" % threshold)
    umi = as_int(row, "unique_umi_families_%dbp" % threshold)
    if raw >= args.main_min_raw_reads and umi >= args.main_min_umi_families:
        return "MAIN"
    if raw >= args.candidate_min_raw_reads and umi >= args.candidate_min_umi_families:
        return "CANDIDATE"
    if umi == 0:
        return "REVIEW_UMI0"
    return "REVIEW_LOW_SUPPORT"


def expected_band(row: dict, primary: int) -> str:
    maximum = as_int(row, "max_supported_threshold")
    if maximum >= 50:
        return "50BP"
    if maximum >= 40:
        return "40BP_ONLY"
    return "%dBP_ONLY" % primary


def validate_old_sequences(metadata_path: str, old_csv: str) -> dict:
    with open(old_csv, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or "ori_seq" not in reader.fieldnames:
            raise ValueError("Old sequence CSV must contain ori_seq")
        old_sequences = {row["ori_seq"].strip().upper() for row in reader if row["ori_seq"].strip()}
    expected_ids = {"SEQ_" + hashlib.sha256(seq.encode("ascii")).hexdigest() for seq in old_sequences}
    found: Set[str] = set()
    with gzip.open(metadata_path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["sequence_id"] in expected_ids:
                found.add(row["sequence_id"])
    missing = expected_ids - found
    if missing:
        raise AssertionError("%d old exact sequences are missing from audit metadata" % len(missing))
    return {"old_distinct_sequences": len(old_sequences), "old_exact_sequences_recovered": len(found)}


def validate_partition(root: Path, prefix: str, kind: str, rows: Sequence[dict], id_key: str) -> None:
    expected = {
        "main": {row[id_key] for row in rows if row["result_tier"] == "MAIN"},
        "candidate": {row[id_key] for row in rows if row["result_tier"] == "CANDIDATE"},
        "review": {row[id_key] for row in rows if row["result_tier"].startswith("REVIEW_")},
    }
    observed = {
        label: {row[id_key] for row in read_tsv(root / (prefix + kind + "." + label + ".tsv"))}
        for label in ("main", "candidate", "review")
    }
    if observed != expected:
        raise AssertionError("%s partition files do not match final table" % kind)


def main() -> None:
    args = parse_args()
    root = Path(args.final_dir)
    prefix = args.sample + "."
    paths = {
        "events": root / (prefix + "fusion_events.final.tsv"),
        "loci": root / (prefix + "fusion_loci.final.tsv"),
        "evidence": root / (prefix + "fusion_evidence.final.tsv.gz"),
        "sequences": root / (prefix + "fusion_unique_sequences.final.tsv.gz"),
        "members": root / (prefix + "fusion_sequence_members.final.tsv.gz"),
        "exact_umi": root / (prefix + "fusion_umi_families.exact.tsv.gz"),
        "corrected_umi": root / (prefix + "fusion_umi_families.corrected.tsv.gz"),
        "summary": root / (prefix + "fusion_summary.final.json"),
    }
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("Missing or empty %s file: %s" % (label, path))

    events = read_tsv(paths["events"])
    loci = read_tsv(paths["loci"])
    evidence = read_tsv(paths["evidence"])
    sequences = read_tsv(paths["sequences"])
    members = read_tsv(paths["members"])
    exact_umi = read_tsv(paths["exact_umi"])
    corrected_umi = read_tsv(paths["corrected_umi"])
    event_by_id = {row["event_id"]: row for row in events}
    sequence_by_id = {row["sequence_id"]: row for row in sequences}
    locus_by_id = {row["locus_id"]: row for row in loci}
    if len(event_by_id) != len(events) or len(sequence_by_id) != len(sequences):
        raise AssertionError("Duplicate event or normalized sequence ID")
    if len(locus_by_id) != len(loci):
        raise AssertionError("Duplicate locus ID")

    member_keys: Set[Tuple[str, str]] = set()
    members_by_sequence: DefaultDict[str, List[dict]] = defaultdict(list)
    for row in members:
        key = (row["sequence_id"], row["original_sequence_id"])
        if key in member_keys:
            raise AssertionError("Duplicate normalized/original sequence membership: %s" % (key,))
        member_keys.add(key)
        if row["sequence_id"] not in sequence_by_id:
            raise AssertionError("Member references missing normalized sequence")
        members_by_sequence[row["sequence_id"]].append(row)

    for sequence_id, row in sequence_by_id.items():
        sequence = row["sequence"]
        expected_id = "SEQ_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()
        if expected_id != sequence_id:
            raise AssertionError("Normalized sequence hash mismatch: %s" % sequence_id)
        values = members_by_sequence[sequence_id]
        if as_int(row, "original_sequence_count") != len(values):
            raise AssertionError("Original-sequence member count mismatch: %s" % sequence_id)
        if as_int(row, "raw_read_count") != sum(as_int(item, "raw_read_count") for item in values):
            raise AssertionError("Normalized sequence raw-read sum mismatch: %s" % sequence_id)
        if as_int(row, "umi_eligible_read_count") != sum(
            as_int(item, "umi_eligible_read_count") for item in values
        ):
            raise AssertionError("Normalized sequence UMI-eligible sum mismatch: %s" % sequence_id)

    primary_flag = "pass_highconf_%dbp" % args.primary_threshold
    evidence_keys: Set[Tuple[str, str]] = set()
    evidence_by_event: DefaultDict[str, List[dict]] = defaultdict(list)
    events_by_sequence: DefaultDict[str, Set[str]] = defaultdict(set)
    for row in evidence:
        key = (row["event_id"], row["sequence_id"])
        if key in evidence_keys:
            raise AssertionError("Duplicate event-normalized-sequence evidence: %s" % (key,))
        evidence_keys.add(key)
        if as_int(row, primary_flag) != 1:
            raise AssertionError("Non-primary-passing evidence in final table")
        if row["event_id"] not in event_by_id or row["sequence_id"] not in sequence_by_id:
            raise AssertionError("Evidence foreign key is missing")
        evidence_by_event[row["event_id"]].append(row)
        events_by_sequence[row["sequence_id"]].add(row["event_id"])

    for sequence_id, row in sequence_by_id.items():
        listed = {item for item in row["event_ids"].split(";") if item}
        if listed != events_by_sequence[sequence_id] or as_int(row, "event_count") != len(listed):
            raise AssertionError("Sequence-to-event mapping mismatch: %s" % sequence_id)
        if as_int(row, "multi_event_sequence") != int(len(listed) > 1):
            raise AssertionError("multi_event_sequence flag mismatch: %s" % sequence_id)

    exact_keys: Set[Tuple[str, str, str]] = set()
    exact_by_event: DefaultDict[str, Set[Tuple[str, str]]] = defaultdict(set)
    for row in exact_umi:
        key = (row["event_id"], row["tube"], row["umi"])
        if key in exact_keys:
            raise AssertionError("Duplicate exact UMI family: %s" % (key,))
        exact_keys.add(key)
        exact_by_event[row["event_id"]].add((row["tube"], row["umi"]))
    corrected_keys: Set[Tuple[str, str, str]] = set()
    corrected_by_event: DefaultDict[str, Set[Tuple[str, str]]] = defaultdict(set)
    for row in corrected_umi:
        key = (row["event_id"], row["tube"], row["corrected_umi"])
        if key in corrected_keys:
            raise AssertionError("Duplicate corrected UMI family: %s" % (key,))
        corrected_keys.add(key)
        corrected_by_event[row["event_id"]].add((row["tube"], row["corrected_umi"]))
    if not set(corrected_by_event).issubset(set(event_by_id)):
        raise AssertionError("Corrected UMI references missing event")

    for event_id, event in event_by_id.items():
        if as_int(event, "human_breakpoint_start0") != as_int(event, "human_breakpoint_1based") - 1:
            raise AssertionError("Human start0 coordinate mismatch: %s" % event_id)
        if as_int(event, "human_breakpoint_end0") != as_int(event, "human_breakpoint_1based"):
            raise AssertionError("Human end0 coordinate mismatch: %s" % event_id)
        if not 1 <= as_int(event, "hpv_breakpoint_1based") <= as_int(event, "hpv_ref_length"):
            raise AssertionError("HPV coordinate outside reference: %s" % event_id)
        rows = evidence_by_event[event_id]
        expected_raw = sum(as_int(sequence_by_id[row["sequence_id"]], "raw_read_count") for row in rows)
        expected_eligible = sum(
            as_int(sequence_by_id[row["sequence_id"]], "umi_eligible_read_count") for row in rows
        )
        if as_int(event, "unique_sequences_%dbp" % args.primary_threshold) != len(rows):
            raise AssertionError("Event normalized-sequence count mismatch: %s" % event_id)
        if as_int(event, "raw_reads_%dbp" % args.primary_threshold) != expected_raw:
            raise AssertionError("Event raw-read count mismatch: %s" % event_id)
        if as_int(event, "umi_eligible_reads_%dbp" % args.primary_threshold) != expected_eligible:
            raise AssertionError("Event UMI-eligible count mismatch: %s" % event_id)
        if as_int(event, "unique_umi_families_%dbp" % args.primary_threshold) != len(exact_by_event[event_id]):
            raise AssertionError("Event exact UMI count mismatch: %s" % event_id)
        if as_int(event, "corrected_umi_families_%dbp" % args.primary_threshold) != len(
            corrected_by_event[event_id]
        ):
            raise AssertionError("Event corrected UMI count mismatch: %s" % event_id)
        if len(corrected_by_event[event_id]) > len(exact_by_event[event_id]):
            raise AssertionError("Corrected UMI count exceeds exact UMI count")
        if event["result_tier"] != expected_tier(event, args.primary_threshold, args):
            raise AssertionError("Event tier mismatch: %s" % event_id)
        if event["support_band"] != expected_band(event, args.primary_threshold):
            raise AssertionError("Event support-band mismatch: %s" % event_id)
        if event["report_class"] != event["result_tier"] + "_" + event["support_band"]:
            raise AssertionError("Event report-class mismatch: %s" % event_id)

    validate_partition(root, prefix, "fusion_events", events, "event_id")
    validate_partition(root, prefix, "fusion_loci", loci, "locus_id")

    event_locus_count: DefaultDict[str, int] = defaultdict(int)
    for locus in loci:
        event_ids = {item for item in locus["exact_event_ids"].split(";") if item}
        if as_int(locus, "exact_event_count") != len(event_ids):
            raise AssertionError("Locus exact-event count mismatch")
        if not event_ids.issubset(set(event_by_id)):
            raise AssertionError("Locus references missing exact event")
        if locus["representative_event_id"] not in event_ids:
            raise AssertionError("Locus representative is not a member event")
        for event_id in event_ids:
            event_locus_count[event_id] += 1
        locus_sequence_ids = {
            row["sequence_id"] for event_id in event_ids for row in evidence_by_event[event_id]
        }
        expected_locus_raw = sum(
            as_int(sequence_by_id[sequence_id], "raw_read_count")
            for sequence_id in locus_sequence_ids
        )
        expected_locus_eligible = sum(
            as_int(sequence_by_id[sequence_id], "umi_eligible_read_count")
            for sequence_id in locus_sequence_ids
        )
        expected_locus_families = set(
            family for event_id in event_ids for family in exact_by_event[event_id]
        )
        if as_int(locus, "unique_sequences_%dbp" % args.primary_threshold) != len(locus_sequence_ids):
            raise AssertionError("Locus normalized-sequence count mismatch")
        if as_int(locus, "raw_reads_%dbp" % args.primary_threshold) != expected_locus_raw:
            raise AssertionError("Locus raw-read count mismatch")
        if as_int(locus, "umi_eligible_reads_%dbp" % args.primary_threshold) != expected_locus_eligible:
            raise AssertionError("Locus UMI-eligible count mismatch")
        if as_int(locus, "unique_umi_families_%dbp" % args.primary_threshold) != len(
            expected_locus_families
        ):
            raise AssertionError("Locus exact UMI-family count mismatch")
        if as_int(locus, "corrected_umi_families_%dbp" % args.primary_threshold) > len(
            expected_locus_families
        ):
            raise AssertionError("Locus corrected UMI count exceeds exact count")
        if as_int(locus, "human_breakpoint_start0") != as_int(locus, "human_breakpoint_1based") - 1:
            raise AssertionError("Locus human start0 mismatch")
        if as_int(locus, "human_breakpoint_end0") != as_int(locus, "human_breakpoint_1based"):
            raise AssertionError("Locus human end0 mismatch")
        if locus["result_tier"] != expected_tier(locus, args.primary_threshold, args):
            raise AssertionError("Locus tier mismatch")
        if locus["support_band"] != expected_band(locus, args.primary_threshold):
            raise AssertionError("Locus support-band mismatch")
    if set(event_locus_count) != set(event_by_id) or any(value != 1 for value in event_locus_count.values()):
        raise AssertionError("Every exact event must belong to exactly one locus")

    with open(paths["summary"], "r", encoding="utf-8") as handle:
        summary = json.load(handle)
    expected_summary = {
        "final_exact_events": len(events), "final_loci": len(loci),
        "final_evidence_rows": len(evidence), "final_unique_sequences": len(sequences),
        "exact_umi_families": len(exact_umi), "corrected_umi_families": len(corrected_umi),
    }
    for key, value in expected_summary.items():
        if summary.get(key) != value:
            raise AssertionError("Summary mismatch for %s" % key)

    report = {
        "status": "PASS", "sample": args.sample,
        "final_exact_events": len(events), "final_loci": len(loci),
        "main_events": sum(row["result_tier"] == "MAIN" for row in events),
        "candidate_events": sum(row["result_tier"] == "CANDIDATE" for row in events),
        "review_events": sum(row["result_tier"].startswith("REVIEW_") for row in events),
        "final_evidence_rows": len(evidence), "final_unique_sequences": len(sequences),
        "original_sequences_represented": len({row["original_sequence_id"] for row in members}),
        "exact_umi_families": len(exact_umi), "corrected_umi_families": len(corrected_umi),
        "checks": [
            "normalized sequence hashes and one-to-many members",
            "event-sequence foreign keys and primary threshold",
            "human/HPV coordinate conventions", "raw-read and exact/corrected UMI aggregation",
            "event and locus tiers/partitions", "exact-event to locus membership", "summary counts",
        ],
    }
    if args.audit_sequence_metadata or args.old_sequence_csv:
        if not args.audit_sequence_metadata or not args.old_sequence_csv:
            raise ValueError("Both old-sequence regression inputs are required together")
        report["old_sequence_regression"] = validate_old_sequences(
            args.audit_sequence_metadata, args.old_sequence_csv
        )
    print(json.dumps(report, ensure_ascii=False, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
