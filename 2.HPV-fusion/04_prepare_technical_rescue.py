#!/usr/bin/env python3
"""Prepare terminally trimmed rescue sequences and collapse identical trimmed bases."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import Counter, defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, List, Sequence, Set, Tuple

from fusion_common import wrap_fasta


csv.field_size_limit(1024 * 1024 * 1024)


METADATA_COLUMNS = [
    "sample", "sequence_id", "sequence_length", "raw_read_count",
    "umi_eligible_read_count", "unique_umi_count", "unique_umi_family_count",
    "umi_family_counts", "tube_A_read_count", "tube_B_read_count",
    "unknown_tube_read_count", "missing_raw_id_count", "ambiguous_read_id_count",
    "technical_intervals0", "technical_bases", "sequence",
]

MAP_COLUMNS = [
    "rescue_sequence_id", "normalized_sequence_id", "original_sequence_id",
    "trim_start0", "trim_end0", "original_length", "trimmed_length",
    "removed_left_bases", "removed_right_bases", "trim_reason", "source_event_ids",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict-evidence", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--out-fasta", required=True)
    parser.add_argument("--out-metadata", required=True)
    parser.add_argument("--out-map", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--threshold", type=int, default=30)
    parser.add_argument("--min-human-mapq", type=int, default=20)
    parser.add_argument("--min-hpv-mapq", type=int, default=0)
    parser.add_argument("--terminal-slack", type=int, default=10)
    parser.add_argument("--min-trimmed-length", type=int, default=60)
    return parser.parse_args()


def as_int(row: dict, key: str) -> int:
    return int(row.get(key, "0") or 0)


def parse_intervals(value: str) -> List[Tuple[int, int]]:
    if not value or value == ".":
        return []
    intervals: List[Tuple[int, int]] = []
    for item in value.split(";"):
        start0, end0 = item.split("-", 1)
        intervals.append((int(start0), int(end0)))
    return sorted(intervals)


def parse_family_counts(value: str) -> Dict[str, int]:
    if not value or value == ".":
        return {}
    counts: Dict[str, int] = {}
    for item in value.split(";"):
        family, count = item.rsplit(":", 1)
        counts[family] = counts.get(family, 0) + int(count)
    return counts


def trim_bounds(
    length: int,
    intervals: Sequence[Tuple[int, int]],
    terminal_slack: int,
) -> Tuple[int, int, str]:
    """Return one contiguous biological segment; never splice across an internal cassette."""
    if not intervals:
        return 0, length, "NO_TECHNICAL_INTERVAL"
    trim_start0 = 0
    trim_end0 = length
    left_terminal = intervals[0][0] <= terminal_slack
    right_terminal = length - intervals[-1][1] <= terminal_slack
    if left_terminal:
        trim_start0 = intervals[0][1]
    if right_terminal:
        trim_end0 = intervals[-1][0]
    if not left_terminal and not right_terminal:
        return 0, length, "INTERNAL_TECHNICAL_INTERVAL"
    for start0, end0 in intervals:
        if max(trim_start0, start0) < min(trim_end0, end0):
            return 0, length, "INTERNAL_TECHNICAL_INTERVAL_REMAINS"
    label = (
        "TRIM_BOTH" if left_terminal and right_terminal
        else "TRIM_LEFT" if left_terminal else "TRIM_RIGHT"
    )
    return trim_start0, trim_end0, label


def rescue_candidates(args: argparse.Namespace) -> Tuple[Set[str], DefaultDict[str, Set[str]], dict]:
    selected: Set[str] = set()
    source_events: DefaultDict[str, Set[str]] = defaultdict(set)
    stats = {
        "strict_evidence_rows": 0,
        "length_passing_but_not_strict_rows": 0,
        "technical_only_rescue_rows": 0,
    }
    length_column = "pass_length_%dbp" % args.threshold
    highconf_column = "pass_highconf_%dbp" % args.threshold
    with gzip.open(args.strict_evidence, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sequence_id", "event_id", length_column, highconf_column,
            "full_human_like", "human_mapq", "human_nh", "hpv_mapq",
            "hpv_technical_overlap_bases", "human_technical_overlap_bases",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError("Strict evidence is missing required rescue columns")
        for row in reader:
            stats["strict_evidence_rows"] += 1
            if not as_int(row, length_column) or as_int(row, highconf_column):
                continue
            stats["length_passing_but_not_strict_rows"] += 1
            technical_overlap = (
                as_int(row, "hpv_technical_overlap_bases")
                + as_int(row, "human_technical_overlap_bases")
            )
            core_mapping_ok = (
                not as_int(row, "full_human_like")
                and as_int(row, "human_mapq") >= args.min_human_mapq
                and as_int(row, "human_nh") == 1
                and as_int(row, "hpv_mapq") >= args.min_hpv_mapq
            )
            if technical_overlap > 0 and core_mapping_ok:
                selected.add(row["sequence_id"])
                source_events[row["sequence_id"]].add(row["event_id"])
                stats["technical_only_rescue_rows"] += 1
    stats["technical_only_rescue_sequences"] = len(selected)
    return selected, source_events, stats


def aggregate_group_metadata(members: Sequence[dict], rescue_id: str, sequence: str) -> dict:
    family_counts: Counter[str] = Counter()
    for member in members:
        family_counts.update(parse_family_counts(member["metadata"]["umi_family_counts"]))
    unique_umis = set(family.split("|", 1)[-1] for family in family_counts)
    first = members[0]["metadata"]
    return {
        "sample": first["sample"],
        "sequence_id": rescue_id,
        "sequence_length": len(sequence),
        "raw_read_count": sum(as_int(member["metadata"], "raw_read_count") for member in members),
        "umi_eligible_read_count": sum(
            as_int(member["metadata"], "umi_eligible_read_count") for member in members
        ),
        "unique_umi_count": len(unique_umis),
        "unique_umi_family_count": len(family_counts),
        "umi_family_counts": ";".join(
            "%s:%d" % (family, count) for family, count in sorted(family_counts.items())
        ) or ".",
        "tube_A_read_count": sum(as_int(member["metadata"], "tube_A_read_count") for member in members),
        "tube_B_read_count": sum(as_int(member["metadata"], "tube_B_read_count") for member in members),
        "unknown_tube_read_count": sum(
            as_int(member["metadata"], "unknown_tube_read_count") for member in members
        ),
        "missing_raw_id_count": sum(
            as_int(member["metadata"], "missing_raw_id_count") for member in members
        ),
        "ambiguous_read_id_count": sum(
            as_int(member["metadata"], "ambiguous_read_id_count") for member in members
        ),
        "technical_intervals0": ".",
        "technical_bases": 0,
        "sequence": sequence,
    }


def main() -> None:
    args = parse_args()
    if args.threshold <= 0 or args.terminal_slack < 0 or args.min_trimmed_length <= 0:
        raise ValueError("Thresholds and lengths must be valid positive values")
    for output in (args.out_fasta, args.out_metadata, args.out_map, args.summary):
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    selected, source_events, stats = rescue_candidates(args)
    groups: DefaultDict[str, List[dict]] = defaultdict(list)
    found: Set[str] = set()
    trim_reasons: Counter[str] = Counter()

    with gzip.open(args.sequence_metadata, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not set(METADATA_COLUMNS).issubset(set(reader.fieldnames)):
            raise ValueError("Sequence metadata is missing required columns")
        for row in reader:
            original_id = row["sequence_id"]
            if original_id not in selected:
                continue
            found.add(original_id)
            sequence = row["sequence"].upper()
            intervals = parse_intervals(row["technical_intervals0"])
            trim_start0, trim_end0, reason = trim_bounds(
                len(sequence), intervals, args.terminal_slack
            )
            if reason.startswith("INTERNAL") or reason == "NO_TECHNICAL_INTERVAL":
                trim_reasons[reason] += 1
                continue
            trimmed = sequence[trim_start0:trim_end0]
            if len(trimmed) < args.min_trimmed_length:
                trim_reasons["TOO_SHORT_AFTER_TRIM"] += 1
                continue
            groups[trimmed].append({
                "metadata": row,
                "trim_start0": trim_start0,
                "trim_end0": trim_end0,
                "trim_reason": reason,
                "source_event_ids": ";".join(sorted(source_events[original_id])),
            })
            trim_reasons[reason] += 1

    missing = selected - found
    if missing:
        raise ValueError("%d rescue candidates are missing from metadata; first=%s" % (
            len(missing), sorted(missing)[0]
        ))

    fasta_tmp = args.out_fasta + ".tmp"
    metadata_tmp = args.out_metadata + ".tmp"
    map_tmp = args.out_map + ".tmp"
    summary_tmp = args.summary + ".tmp"
    original_written = 0
    with open(fasta_tmp, "w", encoding="ascii") as fasta_out, gzip.open(
        metadata_tmp, "wt", encoding="utf-8", newline=""
    ) as metadata_out, open(map_tmp, "w", encoding="utf-8", newline="") as map_out:
        metadata_writer = csv.DictWriter(
            metadata_out, fieldnames=METADATA_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        map_writer = csv.DictWriter(
            map_out, fieldnames=MAP_COLUMNS, delimiter="\t", lineterminator="\n"
        )
        metadata_writer.writeheader()
        map_writer.writeheader()
        for trimmed, members in sorted(
            groups.items(), key=lambda item: hashlib.sha256(item[0].encode("ascii")).hexdigest()
        ):
            digest = hashlib.sha256(trimmed.encode("ascii")).hexdigest()
            rescue_id = "RSEQ_" + digest
            normalized_id = "SEQ_" + digest
            metadata_writer.writerow(aggregate_group_metadata(members, rescue_id, trimmed))
            fasta_out.write(">%s normalized=%s members=%d\n%s\n" % (
                rescue_id, normalized_id, len(members), wrap_fasta(trimmed)
            ))
            for member in sorted(members, key=lambda item: item["metadata"]["sequence_id"]):
                row = member["metadata"]
                original_length = len(row["sequence"])
                map_writer.writerow({
                    "rescue_sequence_id": rescue_id,
                    "normalized_sequence_id": normalized_id,
                    "original_sequence_id": row["sequence_id"],
                    "trim_start0": member["trim_start0"],
                    "trim_end0": member["trim_end0"],
                    "original_length": original_length,
                    "trimmed_length": len(trimmed),
                    "removed_left_bases": member["trim_start0"],
                    "removed_right_bases": original_length - member["trim_end0"],
                    "trim_reason": member["trim_reason"],
                    "source_event_ids": member["source_event_ids"],
                })
                original_written += 1

    unique_trimmed = len(groups)
    stats.update({
        # Backward-compatible key: the number of sequences actually sent to BWA/STAR.
        "rescue_sequences_written": unique_trimmed,
        "rescue_original_sequences_written": original_written,
        "rescue_unique_trimmed_sequences_written": unique_trimmed,
        "collapsed_original_sequence_excess": original_written - unique_trimmed,
        "threshold": args.threshold,
        "terminal_slack": args.terminal_slack,
        "min_trimmed_length": args.min_trimmed_length,
        "trim_reason_counts": dict(sorted(trim_reasons.items())),
        "normalized_sequence_definition": "SEQ_ + sha256(terminally_trimmed_sequence)",
    })
    with open(summary_tmp, "w", encoding="utf-8") as handle:
        json.dump(stats, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(fasta_tmp, args.out_fasta)
    os.replace(metadata_tmp, args.out_metadata)
    os.replace(map_tmp, args.out_map)
    os.replace(summary_tmp, args.summary)
    print("[OK] technical rescue originals=%d unique_trimmed=%d" % (
        original_written, unique_trimmed
    ), flush=True)
    print("[OUT] %s" % args.out_fasta, flush=True)


if __name__ == "__main__":
    main()
