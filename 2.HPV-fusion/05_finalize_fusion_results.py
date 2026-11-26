#!/usr/bin/env python3
"""Build normalized exact-event and nearby-locus HPV-fusion result tables."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import Counter, defaultdict, deque
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Set, Tuple


csv.field_size_limit(1024 * 1024 * 1024)


EVENT_BASE_COLUMNS = [
    "sample", "event_id", "hpv_type", "hpv_ref", "hpv_ref_length",
    "hpv_breakpoint_1based", "hpv_strand", "hpv_breakend_side", "human_chr",
    "human_breakpoint_1based", "human_breakpoint_start0", "human_breakpoint_end0",
    "human_strand", "human_breakend_side", "genes", "query_orientation",
    "max_junction_gap_bases",
]

FINAL_EVIDENCE_BASE_COLUMNS = [
    "sample", "event_id", "sequence_id", "evidence_sources", "normalization_types",
    "aligned_sequence_ids", "original_sequence_count", "strict_original_sequence_count",
    "rescued_original_sequence_count", "multi_event_sequence", "hpv_type", "hpv_ref",
    "hpv_ref_length", "hpv_breakpoint_1based", "hpv_strand", "hpv_breakend_side",
    "human_chr", "human_breakpoint_1based", "human_breakpoint_start0",
    "human_breakpoint_end0", "human_strand", "human_breakend_side", "genes",
    "query_orientation", "junction_gap_bases", "hpv_match_bases", "human_match_bases",
    "hpv_effective_match_bases", "human_effective_match_bases", "hpv_mapq",
    "human_mapq", "human_nh", "hpv_cigar", "human_cigar",
    "full_human_aligned_fraction", "raw_read_count", "umi_eligible_read_count",
    "exact_umi_family_count",
]

FINAL_SEQUENCE_COLUMNS = [
    "sample", "sequence_id", "event_ids", "event_count", "multi_event_sequence",
    "evidence_sources", "normalization_types", "original_sequence_count",
    "sequence_length", "raw_read_count", "umi_eligible_read_count", "unique_umi_count",
    "exact_umi_family_count", "tube_A_read_count", "tube_B_read_count",
    "unknown_tube_read_count", "missing_raw_id_count", "ambiguous_read_id_count",
    "technical_bases_removed_sum", "sequence",
]

MEMBER_COLUMNS = [
    "sample", "sequence_id", "original_sequence_id", "normalization_type",
    "evidence_sources", "trim_start0", "trim_end0", "original_length",
    "normalized_length", "removed_left_bases", "removed_right_bases", "raw_read_count",
    "umi_eligible_read_count", "exact_umi_family_count", "technical_intervals0",
    "technical_bases",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--strict-evidence", required=True)
    parser.add_argument("--rescue-evidence", default="")
    parser.add_argument("--rescue-map", default="")
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--thresholds", default="20,30,40,50")
    parser.add_argument("--primary-threshold", type=int, default=30)
    parser.add_argument("--min-human-mapq", type=int, default=20)
    parser.add_argument("--min-hpv-mapq", type=int, default=0)
    parser.add_argument("--main-min-raw-reads", type=int, default=3)
    parser.add_argument("--main-min-umi-families", type=int, default=2)
    parser.add_argument("--candidate-min-raw-reads", type=int, default=2)
    parser.add_argument("--candidate-min-umi-families", type=int, default=1)
    parser.add_argument("--locus-human-window", type=int, default=5)
    parser.add_argument("--locus-hpv-window", type=int, default=5)
    parser.add_argument("--umi-correction-distance", type=int, default=1)
    parser.add_argument("--pipeline-version", default="3.0")
    return parser.parse_args()


def as_int(row: dict, key: str) -> int:
    return int(row.get(key, "0") or 0)


def parse_family_counts(value: str) -> Dict[str, int]:
    if not value or value == ".":
        return {}
    counts: Counter[str] = Counter()
    for item in value.split(";"):
        family, count = item.rsplit(":", 1)
        counts[family] += int(count)
    return dict(counts)


def open_tsv_auto(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def read_evidence(path: str) -> List[dict]:
    with open_tsv_auto(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_rescue_map(path: str) -> Dict[str, List[dict]]:
    if not path:
        return {}
    grouped: DefaultDict[str, List[dict]] = defaultdict(list)
    seen: Set[Tuple[str, str]] = set()
    with open(path, "r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "rescue_sequence_id", "normalized_sequence_id", "original_sequence_id",
            "trim_start0", "trim_end0", "original_length", "trimmed_length",
            "removed_left_bases", "removed_right_bases", "trim_reason", "source_event_ids",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError("V3 rescue map is missing required columns")
        for row in reader:
            key = (row["rescue_sequence_id"], row["original_sequence_id"])
            if key in seen:
                raise ValueError("Duplicate rescue/original map row: %s" % (key,))
            seen.add(key)
            grouped[row["rescue_sequence_id"]].append(row)
    for rescue_id, rows in grouped.items():
        normalized = {row["normalized_sequence_id"] for row in rows}
        if len(normalized) != 1:
            raise ValueError("Rescue group maps to multiple normalized sequences: %s" % rescue_id)
    return dict(grouped)


def load_original_metadata(path: str, needed: Set[str]) -> Dict[str, dict]:
    result: Dict[str, dict] = {}
    if not needed:
        return result
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sample", "sequence_id", "sequence_length", "raw_read_count",
            "umi_eligible_read_count", "unique_umi_count", "unique_umi_family_count",
            "umi_family_counts", "tube_A_read_count", "tube_B_read_count",
            "unknown_tube_read_count", "missing_raw_id_count", "ambiguous_read_id_count",
            "technical_intervals0", "technical_bases", "sequence",
        }
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError("Sequence metadata is missing required columns")
        for row in reader:
            if row["sequence_id"] not in needed:
                continue
            row["family_counts"] = parse_family_counts(row["umi_family_counts"])
            result[row["sequence_id"]] = row
    missing = needed - set(result)
    if missing:
        raise ValueError("%d required original sequences are absent from metadata; first=%s" % (
            len(missing), sorted(missing)[0]
        ))
    return result


def rejection_reason(row: dict, threshold: int, args: argparse.Namespace) -> str:
    if as_int(row, "pass_highconf_%dbp" % threshold):
        return "STRICT_PASS"
    if not as_int(row, "pass_length_%dbp" % threshold):
        return "EFFECTIVE_ARM_BELOW_%dBP" % threshold
    if as_int(row, "full_human_like"):
        return "FULL_LENGTH_HUMAN_LIKE"
    if as_int(row, "human_nh") != 1:
        return "HUMAN_MULTIMAP_NH_NOT_1"
    if as_int(row, "human_mapq") < args.min_human_mapq:
        return "HUMAN_MAPQ_BELOW_%d" % args.min_human_mapq
    if as_int(row, "hpv_mapq") < args.min_hpv_mapq:
        return "HPV_MAPQ_BELOW_%d" % args.min_hpv_mapq
    if as_int(row, "hpv_technical_overlap_bases") or as_int(row, "human_technical_overlap_bases"):
        return "TERMINAL_TECHNICAL_OVERLAP_RESCUE_CANDIDATE"
    return "OTHER_MAPPING_FILTER"


def evidence_rank(row: dict, thresholds: Sequence[int], source: str) -> Tuple[int, int, int, int, int]:
    return (
        1 if source == "STRICT" else 0,
        sum(as_int(row, "pass_highconf_%dbp" % threshold) for threshold in thresholds),
        as_int(row, "human_mapq"),
        min(as_int(row, "hpv_effective_match_bases"), as_int(row, "human_effective_match_bases")),
        -as_int(row, "junction_gap_bases"),
    )


def normalize_evidence(
    strict_rows: Sequence[dict],
    rescue_rows: Sequence[dict],
    rescue_map: Dict[str, List[dict]],
    thresholds: Sequence[int],
    primary_threshold: int,
) -> Tuple[List[dict], Set[str]]:
    primary_column = "pass_highconf_%dbp" % primary_threshold
    grouped: Dict[Tuple[str, str], dict] = {}
    needed_originals: Set[str] = set()

    def consider(
        row: dict,
        source: str,
        normalized_id: str,
        original_ids: Set[str],
        normalization_type: str,
    ) -> None:
        if not as_int(row, primary_column):
            return
        needed_originals.update(original_ids)
        key = (row["event_id"], normalized_id)
        state = grouped.get(key)
        rank = evidence_rank(row, thresholds, source)
        if state is None:
            state = {
                "best": dict(row), "best_rank": rank, "sources": set(),
                "normalization_types": set(), "aligned_ids": set(),
                "original_ids": set(), "original_ids_by_source": defaultdict(set),
                "pass": dict((threshold, 0) for threshold in thresholds),
            }
            grouped[key] = state
        elif rank > state["best_rank"]:
            state["best"] = dict(row)
            state["best_rank"] = rank
        state["sources"].add(source)
        state["normalization_types"].add(normalization_type)
        state["aligned_ids"].add(row["sequence_id"])
        state["original_ids"].update(original_ids)
        state["original_ids_by_source"][source].update(original_ids)
        for threshold in thresholds:
            state["pass"][threshold] = max(
                state["pass"][threshold], as_int(row, "pass_highconf_%dbp" % threshold)
            )

    for row in strict_rows:
        consider(row, "STRICT", row["sequence_id"], {row["sequence_id"]}, "FULL_LENGTH")
    for row in rescue_rows:
        mappings = rescue_map.get(row["sequence_id"])
        if mappings is None:
            raise ValueError("Rescue evidence ID absent from V3 rescue map: %s" % row["sequence_id"])
        normalized_ids = {mapping["normalized_sequence_id"] for mapping in mappings}
        consider(
            row, "RESCUED_AFTER_TERMINAL_TRIM", next(iter(normalized_ids)),
            {mapping["original_sequence_id"] for mapping in mappings}, "TERMINAL_TRIMMED",
        )

    output: List[dict] = []
    for (event_id, normalized_id), state in sorted(grouped.items()):
        out = dict(state["best"])
        out["event_id"] = event_id
        out["sequence_id"] = normalized_id
        out["evidence_sources"] = ";".join(sorted(state["sources"]))
        out["normalization_types"] = ";".join(sorted(state["normalization_types"]))
        out["aligned_sequence_ids"] = ";".join(sorted(state["aligned_ids"]))
        out["original_sequence_count"] = len(state["original_ids"])
        out["strict_original_sequence_count"] = len(state["original_ids_by_source"]["STRICT"])
        out["rescued_original_sequence_count"] = len(
            state["original_ids_by_source"]["RESCUED_AFTER_TERMINAL_TRIM"]
        )
        for threshold in thresholds:
            out["pass_highconf_%dbp" % threshold] = state["pass"][threshold]
        out["_original_ids"] = state["original_ids"]
        out["_original_ids_by_source"] = state["original_ids_by_source"]
        output.append(out)
    return output, needed_originals


def normalized_sequence_bases(
    sequence_id: str,
    rows: Sequence[dict],
    original_metadata: Dict[str, dict],
    rescue_map: Dict[str, List[dict]],
) -> str:
    candidates: Set[str] = set()
    for row in rows:
        for original_id in row["_original_ids_by_source"]["STRICT"]:
            candidates.add(original_metadata[original_id]["sequence"].upper())
        for aligned_id in row["aligned_sequence_ids"].split(";"):
            for mapping in rescue_map.get(aligned_id, []):
                original = original_metadata[mapping["original_sequence_id"]]["sequence"].upper()
                candidates.add(original[as_int(mapping, "trim_start0"):as_int(mapping, "trim_end0")])
    if len(candidates) != 1:
        raise ValueError("Normalized sequence ID resolves to %d different base strings: %s" % (
            len(candidates), sequence_id
        ))
    sequence = next(iter(candidates))
    expected_id = "SEQ_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()
    if expected_id != sequence_id:
        raise ValueError("Normalized sequence hash mismatch: %s != %s" % (sequence_id, expected_id))
    return sequence


def aggregate_metadata(original_ids: Set[str], metadata: Dict[str, dict], sequence: str) -> dict:
    families: Counter[str] = Counter()
    for original_id in original_ids:
        families.update(metadata[original_id]["family_counts"])
    unique_umis = {family.split("|", 1)[-1] for family in families}
    first = metadata[next(iter(original_ids))]
    return {
        "sample": first["sample"],
        "sequence_length": len(sequence),
        "raw_read_count": sum(as_int(metadata[item], "raw_read_count") for item in original_ids),
        "umi_eligible_read_count": sum(
            as_int(metadata[item], "umi_eligible_read_count") for item in original_ids
        ),
        "unique_umi_count": len(unique_umis),
        "exact_umi_family_count": len(families),
        "family_counts": dict(families),
        "tube_A_read_count": sum(as_int(metadata[item], "tube_A_read_count") for item in original_ids),
        "tube_B_read_count": sum(as_int(metadata[item], "tube_B_read_count") for item in original_ids),
        "unknown_tube_read_count": sum(
            as_int(metadata[item], "unknown_tube_read_count") for item in original_ids
        ),
        "missing_raw_id_count": sum(
            as_int(metadata[item], "missing_raw_id_count") for item in original_ids
        ),
        "ambiguous_read_id_count": sum(
            as_int(metadata[item], "ambiguous_read_id_count") for item in original_ids
        ),
        "sequence": sequence,
    }


def build_sequence_metadata(
    evidence: Sequence[dict],
    original_metadata: Dict[str, dict],
    rescue_map: Dict[str, List[dict]],
) -> Tuple[Dict[str, dict], List[dict]]:
    rows_by_sequence: DefaultDict[str, List[dict]] = defaultdict(list)
    for row in evidence:
        rows_by_sequence[row["sequence_id"]].append(row)
    normalized: Dict[str, dict] = {}
    member_rows: List[dict] = []
    map_by_pair: Dict[Tuple[str, str], dict] = {}
    for mappings in rescue_map.values():
        for mapping in mappings:
            map_by_pair[(mapping["normalized_sequence_id"], mapping["original_sequence_id"])] = mapping

    for sequence_id, rows in sorted(rows_by_sequence.items()):
        original_ids = set(item for row in rows for item in row["_original_ids"])
        sequence = normalized_sequence_bases(sequence_id, rows, original_metadata, rescue_map)
        meta = aggregate_metadata(original_ids, original_metadata, sequence)
        meta["sequence_id"] = sequence_id
        normalized[sequence_id] = meta
        sources_by_original: DefaultDict[str, Set[str]] = defaultdict(set)
        for row in rows:
            for source, ids in row["_original_ids_by_source"].items():
                for original_id in ids:
                    sources_by_original[original_id].add(source)
        for original_id in sorted(original_ids):
            original = original_metadata[original_id]
            mapping = map_by_pair.get((sequence_id, original_id))
            if mapping is None:
                trim_start0, trim_end0 = 0, len(original["sequence"])
                normalization_type = "FULL_LENGTH"
                removed_left = removed_right = 0
            else:
                trim_start0 = as_int(mapping, "trim_start0")
                trim_end0 = as_int(mapping, "trim_end0")
                normalization_type = "TERMINAL_TRIMMED"
                removed_left = as_int(mapping, "removed_left_bases")
                removed_right = as_int(mapping, "removed_right_bases")
            member_rows.append({
                "sample": original["sample"],
                "sequence_id": sequence_id,
                "original_sequence_id": original_id,
                "normalization_type": normalization_type,
                "evidence_sources": ";".join(sorted(sources_by_original[original_id])),
                "trim_start0": trim_start0,
                "trim_end0": trim_end0,
                "original_length": len(original["sequence"]),
                "normalized_length": len(sequence),
                "removed_left_bases": removed_left,
                "removed_right_bases": removed_right,
                "raw_read_count": original["raw_read_count"],
                "umi_eligible_read_count": original["umi_eligible_read_count"],
                "exact_umi_family_count": len(original["family_counts"]),
                "technical_intervals0": original["technical_intervals0"],
                "technical_bases": original["technical_bases"],
            })
    return normalized, member_rows


def hamming_distance(left: str, right: str, max_distance: int) -> int:
    if len(left) != len(right):
        return max_distance + 1
    distance = 0
    for a, b in zip(left, right):
        if a != b:
            distance += 1
            if distance > max_distance:
                break
    return distance


def directional_umi_assignments(counts: Dict[str, int], max_distance: int) -> Dict[str, str]:
    """Count-aware directional clustering; exact UMI counts remain the primary result."""
    if max_distance <= 0:
        return dict((umi, umi) for umi in counts)
    ordered = sorted(counts, key=lambda umi: (-counts[umi], umi))
    assigned: Dict[str, str] = {}
    for root in ordered:
        if root in assigned:
            continue
        assigned[root] = root
        queue = deque([root])
        while queue:
            parent = queue.popleft()
            for child in ordered:
                if child in assigned:
                    continue
                if counts[parent] < 2 * counts[child] - 1:
                    continue
                if hamming_distance(parent, child, max_distance) <= max_distance:
                    assigned[child] = root
                    queue.append(child)
    return assigned


def corrected_family_count(family_counts: Dict[str, int], max_distance: int) -> int:
    by_tube: DefaultDict[str, Dict[str, int]] = defaultdict(dict)
    for family, count in family_counts.items():
        tube, umi = family.split("|", 1) if "|" in family else ("UNKNOWN", family)
        by_tube[tube][umi] = by_tube[tube].get(umi, 0) + count
    return sum(len(set(directional_umi_assignments(counts, max_distance).values())) for counts in by_tube.values())


def event_tier(raw_reads: int, exact_umi_families: int, args: argparse.Namespace) -> str:
    if raw_reads >= args.main_min_raw_reads and exact_umi_families >= args.main_min_umi_families:
        return "MAIN"
    if raw_reads >= args.candidate_min_raw_reads and exact_umi_families >= args.candidate_min_umi_families:
        return "CANDIDATE"
    if exact_umi_families == 0:
        return "REVIEW_UMI0"
    return "REVIEW_LOW_SUPPORT"


def support_band(max_threshold: int, primary: int) -> str:
    if max_threshold >= 50:
        return "50BP"
    if max_threshold >= 40:
        return "40BP_ONLY"
    return "%dBP_ONLY" % primary


def family_counts_for_rows(rows: Iterable[dict], sequence_metadata: Dict[str, dict]) -> Counter[str]:
    counts: Counter[str] = Counter()
    for row in rows:
        counts.update(sequence_metadata[row["sequence_id"]]["family_counts"])
    return counts


def aggregate_events(
    evidence_rows: Sequence[dict],
    sequence_metadata: Dict[str, dict],
    thresholds: Sequence[int],
    args: argparse.Namespace,
) -> Tuple[List[dict], List[dict], List[dict]]:
    grouped: DefaultDict[str, List[dict]] = defaultdict(list)
    events_by_sequence: DefaultDict[str, Set[str]] = defaultdict(set)
    for row in evidence_rows:
        grouped[row["event_id"]].append(row)
        events_by_sequence[row["sequence_id"]].add(row["event_id"])
    for row in evidence_rows:
        row["multi_event_sequence"] = int(len(events_by_sequence[row["sequence_id"]]) > 1)

    events: List[dict] = []
    exact_umi_rows: List[dict] = []
    corrected_umi_rows: List[dict] = []
    primary = args.primary_threshold
    for event_id, rows in grouped.items():
        first = rows[0]
        event = dict((key, first[key]) for key in EVENT_BASE_COLUMNS if key != "max_junction_gap_bases")
        genes = {gene for row in rows for gene in row["genes"].split("|") if gene}
        if len(genes) > 1:
            genes.discard("INTERGENIC")
        event["genes"] = "|".join(sorted(genes)) or "INTERGENIC"
        event["max_junction_gap_bases"] = max(as_int(row, "junction_gap_bases") for row in rows)
        sources = {source for row in rows for source in row["evidence_sources"].split(";") if source}
        event["evidence_sources"] = ";".join(sorted(sources))

        for threshold in thresholds:
            passing = [row for row in rows if as_int(row, "pass_highconf_%dbp" % threshold)]
            families = family_counts_for_rows(passing, sequence_metadata)
            event["unique_sequences_%dbp" % threshold] = len(passing)
            event["raw_reads_%dbp" % threshold] = sum(
                as_int(sequence_metadata[row["sequence_id"]], "raw_read_count") for row in passing
            )
            event["umi_eligible_reads_%dbp" % threshold] = sum(
                as_int(sequence_metadata[row["sequence_id"]], "umi_eligible_read_count")
                for row in passing
            )
            event["unique_umi_families_%dbp" % threshold] = len(families)
            event["corrected_umi_families_%dbp" % threshold] = corrected_family_count(
                families, args.umi_correction_distance
            )

        for source_name, source_value in (
            ("strict", "STRICT"), ("rescued", "RESCUED_AFTER_TERMINAL_TRIM")
        ):
            source_rows = [
                row for row in rows
                if row["_original_ids_by_source"][source_value]
                and as_int(row, "pass_highconf_%dbp" % primary)
            ]
            source_originals = set(
                original_id for row in source_rows
                for original_id in row["_original_ids_by_source"][source_value]
            )
            source_family_counts: Counter[str] = Counter()
            source_raw = 0
            for row in source_rows:
                for original_id in row["_original_ids_by_source"][source_value]:
                    # Original IDs belonging to the same normalized row are disjoint by construction.
                    source_raw += as_int(row["_original_metadata"][original_id], "raw_read_count")
                    source_family_counts.update(row["_original_metadata"][original_id]["family_counts"])
            event["%s_unique_sequences_%dbp" % (source_name, primary)] = len(source_rows)
            event["%s_original_sequences_%dbp" % (source_name, primary)] = len(source_originals)
            event["%s_raw_reads_%dbp" % (source_name, primary)] = source_raw
            event["%s_unique_umi_families_%dbp" % (source_name, primary)] = len(source_family_counts)

        max_supported = max(
            threshold for threshold in thresholds
            if as_int(event, "unique_sequences_%dbp" % threshold) > 0
        )
        event["max_supported_threshold"] = max_supported
        event["support_band"] = support_band(max_supported, primary)
        event["result_tier"] = event_tier(
            as_int(event, "raw_reads_%dbp" % primary),
            as_int(event, "unique_umi_families_%dbp" % primary), args,
        )
        event["report_class"] = "%s_%s" % (event["result_tier"], event["support_band"])
        event["ambiguous_sequence_count"] = sum(as_int(row, "multi_event_sequence") for row in rows)
        event["requires_manual_review"] = int(
            event["result_tier"].startswith("REVIEW_")
            or event["support_band"] == "%dBP_ONLY" % primary
            or event["ambiguous_sequence_count"] > 0
        )
        events.append(event)

        family_support: DefaultDict[str, dict] = defaultdict(
            lambda: {"read_count": 0, "sequence_ids": set(), "sources": set()}
        )
        for row in rows:
            meta = sequence_metadata[row["sequence_id"]]
            for family, count in meta["family_counts"].items():
                support = family_support[family]
                support["read_count"] += count
                support["sequence_ids"].add(row["sequence_id"])
                for source, original_ids in row["_original_ids_by_source"].items():
                    if any(family in row["_original_metadata"][item]["family_counts"] for item in original_ids):
                        support["sources"].add(source)
        by_tube: DefaultDict[str, Dict[str, int]] = defaultdict(dict)
        for family, support in family_support.items():
            tube, umi = family.split("|", 1) if "|" in family else ("UNKNOWN", family)
            by_tube[tube][umi] = support["read_count"]
        assignments = {
            (tube, umi): root
            for tube, counts in by_tube.items()
            for umi, root in directional_umi_assignments(counts, args.umi_correction_distance).items()
        }
        corrected_support: DefaultDict[Tuple[str, str], dict] = defaultdict(
            lambda: {"read_count": 0, "sequence_ids": set(), "sources": set(), "umis": set()}
        )
        for family, support in sorted(family_support.items()):
            tube, umi = family.split("|", 1) if "|" in family else ("UNKNOWN", family)
            corrected = assignments[(tube, umi)]
            exact_umi_rows.append({
                "sample": first["sample"], "event_id": event_id, "tube": tube, "umi": umi,
                "family_read_count": support["read_count"],
                "supporting_sequence_count": len(support["sequence_ids"]),
                "evidence_sources": ";".join(sorted(support["sources"])),
                "corrected_umi": corrected,
            })
            target = corrected_support[(tube, corrected)]
            target["read_count"] += support["read_count"]
            target["sequence_ids"].update(support["sequence_ids"])
            target["sources"].update(support["sources"])
            target["umis"].add(umi)
        for (tube, corrected), support in sorted(corrected_support.items()):
            corrected_umi_rows.append({
                "sample": first["sample"], "event_id": event_id, "tube": tube,
                "corrected_umi": corrected, "member_exact_umi_count": len(support["umis"]),
                "member_exact_umis": ";".join(sorted(support["umis"])),
                "family_read_count": support["read_count"],
                "supporting_sequence_count": len(support["sequence_ids"]),
                "evidence_sources": ";".join(sorted(support["sources"])),
            })

    events.sort(key=lambda row: (
        {"MAIN": 0, "CANDIDATE": 1, "REVIEW_LOW_SUPPORT": 2, "REVIEW_UMI0": 3}[row["result_tier"]],
        -as_int(row, "unique_umi_families_%dbp" % primary),
        -as_int(row, "raw_reads_%dbp" % primary), row["human_chr"],
        as_int(row, "human_breakpoint_1based"),
    ))
    return events, exact_umi_rows, corrected_umi_rows


def event_columns(thresholds: Sequence[int], primary: int) -> List[str]:
    source_counts = [
        "%s_%s_%dbp" % (source, measure, primary)
        for source in ("strict", "rescued")
        for measure in (
            "unique_sequences", "original_sequences", "raw_reads", "unique_umi_families"
        )
    ]
    counts = [
        name for threshold in thresholds
        for name in (
            "unique_sequences_%dbp" % threshold, "raw_reads_%dbp" % threshold,
            "umi_eligible_reads_%dbp" % threshold, "unique_umi_families_%dbp" % threshold,
            "corrected_umi_families_%dbp" % threshold,
        )
    ]
    return EVENT_BASE_COLUMNS + [
        "evidence_sources", "result_tier", "support_band", "report_class",
        "max_supported_threshold", "ambiguous_sequence_count", "requires_manual_review",
    ] + source_counts + counts


def circular_distance(left: int, right: int, length: int) -> int:
    direct = abs(left - right)
    return min(direct, length - direct) if length > 0 else direct


def cluster_events(events: Sequence[dict], args: argparse.Namespace) -> List[List[dict]]:
    grouped: DefaultDict[Tuple[str, ...], List[dict]] = defaultdict(list)
    for event in events:
        key = (
            event["hpv_type"], event["hpv_ref"], event["hpv_strand"], event["hpv_breakend_side"],
            event["human_chr"], event["human_strand"], event["human_breakend_side"],
            event["query_orientation"],
        )
        grouped[key].append(event)
    clusters: List[List[dict]] = []
    for values in grouped.values():
        local: List[List[dict]] = []
        for event in sorted(values, key=lambda row: (
            as_int(row, "human_breakpoint_1based"), as_int(row, "hpv_breakpoint_1based"), row["event_id"]
        )):
            eligible: List[Tuple[int, int]] = []
            for index, cluster in enumerate(local):
                if all(
                    abs(as_int(event, "human_breakpoint_1based") - as_int(other, "human_breakpoint_1based"))
                    <= args.locus_human_window
                    and circular_distance(
                        as_int(event, "hpv_breakpoint_1based"),
                        as_int(other, "hpv_breakpoint_1based"),
                        as_int(event, "hpv_ref_length"),
                    ) <= args.locus_hpv_window
                    for other in cluster
                ):
                    representative = cluster[0]
                    distance = (
                        abs(as_int(event, "human_breakpoint_1based") - as_int(representative, "human_breakpoint_1based"))
                        + circular_distance(
                            as_int(event, "hpv_breakpoint_1based"),
                            as_int(representative, "hpv_breakpoint_1based"),
                            as_int(event, "hpv_ref_length"),
                        )
                    )
                    eligible.append((distance, index))
            if eligible:
                local[min(eligible)[1]].append(event)
            else:
                local.append([event])
        clusters.extend(local)
    return clusters


def aggregate_loci(
    events: Sequence[dict], evidence: Sequence[dict], sequence_metadata: Dict[str, dict],
    thresholds: Sequence[int], args: argparse.Namespace,
) -> List[dict]:
    evidence_by_event: DefaultDict[str, List[dict]] = defaultdict(list)
    for row in evidence:
        evidence_by_event[row["event_id"]].append(row)
    loci: List[dict] = []
    primary = args.primary_threshold
    for cluster in cluster_events(events, args):
        event_ids = {event["event_id"] for event in cluster}
        representative = max(cluster, key=lambda row: (
            as_int(row, "unique_umi_families_%dbp" % primary),
            as_int(row, "raw_reads_%dbp" % primary), row["event_id"],
        ))
        digest = hashlib.sha256(";".join(sorted(event_ids)).encode("ascii")).hexdigest()[:20]
        rows = [row for event_id in event_ids for row in evidence_by_event[event_id]]
        genes = {gene for event in cluster for gene in event["genes"].split("|") if gene}
        if len(genes) > 1:
            genes.discard("INTERGENIC")
        locus = {
            "sample": representative["sample"], "locus_id": "LOC_" + digest,
            "exact_event_ids": ";".join(sorted(event_ids)), "exact_event_count": len(event_ids),
            "representative_event_id": representative["event_id"],
            "hpv_type": representative["hpv_type"], "hpv_ref": representative["hpv_ref"],
            "hpv_ref_length": representative["hpv_ref_length"],
            "hpv_breakpoint_1based": representative["hpv_breakpoint_1based"],
            "hpv_breakpoint_min1": min(as_int(item, "hpv_breakpoint_1based") for item in cluster),
            "hpv_breakpoint_max1": max(as_int(item, "hpv_breakpoint_1based") for item in cluster),
            "hpv_strand": representative["hpv_strand"],
            "hpv_breakend_side": representative["hpv_breakend_side"],
            "human_chr": representative["human_chr"],
            "human_breakpoint_1based": representative["human_breakpoint_1based"],
            "human_breakpoint_min1": min(as_int(item, "human_breakpoint_1based") for item in cluster),
            "human_breakpoint_max1": max(as_int(item, "human_breakpoint_1based") for item in cluster),
            "human_breakpoint_start0": as_int(representative, "human_breakpoint_1based") - 1,
            "human_breakpoint_end0": representative["human_breakpoint_1based"],
            "human_strand": representative["human_strand"],
            "human_breakend_side": representative["human_breakend_side"],
            "genes": "|".join(sorted(genes)) or "INTERGENIC",
            "query_orientation": representative["query_orientation"],
            "evidence_sources": ";".join(sorted({
                source for item in cluster for source in item["evidence_sources"].split(";") if source
            })),
        }
        for threshold in thresholds:
            by_sequence: Dict[str, dict] = {}
            for row in rows:
                if not as_int(row, "pass_highconf_%dbp" % threshold):
                    continue
                by_sequence.setdefault(row["sequence_id"], row)
            passing = list(by_sequence.values())
            families = family_counts_for_rows(passing, sequence_metadata)
            locus["unique_sequences_%dbp" % threshold] = len(passing)
            locus["raw_reads_%dbp" % threshold] = sum(
                as_int(sequence_metadata[row["sequence_id"]], "raw_read_count") for row in passing
            )
            locus["umi_eligible_reads_%dbp" % threshold] = sum(
                as_int(sequence_metadata[row["sequence_id"]], "umi_eligible_read_count")
                for row in passing
            )
            locus["unique_umi_families_%dbp" % threshold] = len(families)
            locus["corrected_umi_families_%dbp" % threshold] = corrected_family_count(
                families, args.umi_correction_distance
            )
        max_supported = max(
            threshold for threshold in thresholds
            if as_int(locus, "unique_sequences_%dbp" % threshold) > 0
        )
        locus["max_supported_threshold"] = max_supported
        locus["support_band"] = support_band(max_supported, primary)
        locus["result_tier"] = event_tier(
            as_int(locus, "raw_reads_%dbp" % primary),
            as_int(locus, "unique_umi_families_%dbp" % primary), args,
        )
        locus["report_class"] = "%s_%s" % (locus["result_tier"], locus["support_band"])
        locus["ambiguous_sequence_count"] = len({
            row["sequence_id"] for row in rows if as_int(row, "multi_event_sequence")
        })
        locus["requires_manual_review"] = int(
            locus["result_tier"].startswith("REVIEW_")
            or locus["support_band"] == "%dBP_ONLY" % primary
            or locus["ambiguous_sequence_count"] > 0
        )
        loci.append(locus)
    loci.sort(key=lambda row: (
        {"MAIN": 0, "CANDIDATE": 1, "REVIEW_LOW_SUPPORT": 2, "REVIEW_UMI0": 3}[row["result_tier"]],
        -as_int(row, "unique_umi_families_%dbp" % primary),
        -as_int(row, "raw_reads_%dbp" % primary), row["human_chr"],
        as_int(row, "human_breakpoint_1based"),
    ))
    return loci


def locus_columns(thresholds: Sequence[int]) -> List[str]:
    base = [
        "sample", "locus_id", "exact_event_ids", "exact_event_count", "representative_event_id",
        "hpv_type", "hpv_ref", "hpv_ref_length", "hpv_breakpoint_1based",
        "hpv_breakpoint_min1", "hpv_breakpoint_max1", "hpv_strand", "hpv_breakend_side",
        "human_chr", "human_breakpoint_1based", "human_breakpoint_min1",
        "human_breakpoint_max1", "human_breakpoint_start0", "human_breakpoint_end0",
        "human_strand", "human_breakend_side", "genes", "query_orientation",
        "evidence_sources", "result_tier", "support_band", "report_class",
        "max_supported_threshold", "ambiguous_sequence_count", "requires_manual_review",
    ]
    counts = [
        name for threshold in thresholds
        for name in (
            "unique_sequences_%dbp" % threshold, "raw_reads_%dbp" % threshold,
            "umi_eligible_reads_%dbp" % threshold, "unique_umi_families_%dbp" % threshold,
            "corrected_umi_families_%dbp" % threshold,
        )
    ]
    return base + counts


def write_tsv(path: Path, rows: Iterable[dict], columns: Sequence[str], gzip_output: bool = False) -> None:
    tmp = Path(str(path) + ".tmp")
    opener = gzip.open if gzip_output else open
    mode = "wt" if gzip_output else "w"
    with opener(tmp, mode, encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=columns, delimiter="\t", lineterminator="\n", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def main() -> None:
    args = parse_args()
    thresholds = sorted(set(int(item) for item in args.thresholds.split(",") if item.strip()))
    if args.primary_threshold not in thresholds:
        raise ValueError("Primary threshold must be included in thresholds")
    if args.locus_human_window < 0 or args.locus_hpv_window < 0:
        raise ValueError("Locus windows cannot be negative")
    report_thresholds = [value for value in thresholds if value >= args.primary_threshold]
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    strict_rows = read_evidence(args.strict_evidence)
    rescue_rows = read_evidence(args.rescue_evidence) if args.rescue_evidence else []
    rescue_map = load_rescue_map(args.rescue_map) if args.rescue_map else {}
    evidence, needed_originals = normalize_evidence(
        strict_rows, rescue_rows, rescue_map, thresholds, args.primary_threshold
    )
    original_metadata = load_original_metadata(args.sequence_metadata, needed_originals)
    for row in evidence:
        row["_original_metadata"] = original_metadata
    sequence_metadata, member_rows = build_sequence_metadata(evidence, original_metadata, rescue_map)
    events, exact_umi_rows, corrected_umi_rows = aggregate_events(
        evidence, sequence_metadata, report_thresholds, args
    )
    loci = aggregate_loci(events, evidence, sequence_metadata, report_thresholds, args)

    events_by_sequence: DefaultDict[str, Set[str]] = defaultdict(set)
    sources_by_sequence: DefaultDict[str, Set[str]] = defaultdict(set)
    normalization_by_sequence: DefaultDict[str, Set[str]] = defaultdict(set)
    original_ids_by_sequence: DefaultDict[str, Set[str]] = defaultdict(set)
    for row in evidence:
        events_by_sequence[row["sequence_id"]].add(row["event_id"])
        sources_by_sequence[row["sequence_id"]].update(row["evidence_sources"].split(";"))
        normalization_by_sequence[row["sequence_id"]].update(row["normalization_types"].split(";"))
        original_ids_by_sequence[row["sequence_id"]].update(row["_original_ids"])

    removed_by_sequence: Counter[str] = Counter()
    for member in member_rows:
        removed_by_sequence[member["sequence_id"]] += (
            as_int(member, "removed_left_bases") + as_int(member, "removed_right_bases")
        )
    sequence_rows: List[dict] = []
    for sequence_id, meta in sorted(sequence_metadata.items()):
        event_ids = events_by_sequence[sequence_id]
        sequence_rows.append({
            "sample": meta["sample"], "sequence_id": sequence_id,
            "event_ids": ";".join(sorted(event_ids)), "event_count": len(event_ids),
            "multi_event_sequence": int(len(event_ids) > 1),
            "evidence_sources": ";".join(sorted(sources_by_sequence[sequence_id])),
            "normalization_types": ";".join(sorted(normalization_by_sequence[sequence_id])),
            "original_sequence_count": len(original_ids_by_sequence[sequence_id]),
            "sequence_length": meta["sequence_length"], "raw_read_count": meta["raw_read_count"],
            "umi_eligible_read_count": meta["umi_eligible_read_count"],
            "unique_umi_count": meta["unique_umi_count"],
            "exact_umi_family_count": meta["exact_umi_family_count"],
            "tube_A_read_count": meta["tube_A_read_count"],
            "tube_B_read_count": meta["tube_B_read_count"],
            "unknown_tube_read_count": meta["unknown_tube_read_count"],
            "missing_raw_id_count": meta["missing_raw_id_count"],
            "ambiguous_read_id_count": meta["ambiguous_read_id_count"],
            "technical_bases_removed_sum": removed_by_sequence[sequence_id],
            "sequence": meta["sequence"],
        })

    for row in evidence:
        meta = sequence_metadata[row["sequence_id"]]
        row["raw_read_count"] = meta["raw_read_count"]
        row["umi_eligible_read_count"] = meta["umi_eligible_read_count"]
        row["exact_umi_family_count"] = meta["exact_umi_family_count"]

    evidence_columns = FINAL_EVIDENCE_BASE_COLUMNS + [
        "pass_highconf_%dbp" % threshold for threshold in thresholds
    ]
    write_tsv(
        out_dir / (args.sample + ".fusion_evidence.final.tsv.gz"),
        evidence, evidence_columns, gzip_output=True,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_unique_sequences.final.tsv.gz"),
        sequence_rows, FINAL_SEQUENCE_COLUMNS, gzip_output=True,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_sequence_members.final.tsv.gz"),
        member_rows, MEMBER_COLUMNS, gzip_output=True,
    )

    ecolumns = event_columns(report_thresholds, args.primary_threshold)
    write_tsv(out_dir / (args.sample + ".fusion_events.final.tsv"), events, ecolumns)
    write_tsv(
        out_dir / (args.sample + ".fusion_events.main.tsv"),
        (row for row in events if row["result_tier"] == "MAIN"), ecolumns,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_events.candidate.tsv"),
        (row for row in events if row["result_tier"] == "CANDIDATE"), ecolumns,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_events.review.tsv"),
        (row for row in events if row["result_tier"].startswith("REVIEW_")), ecolumns,
    )

    lcolumns = locus_columns(report_thresholds)
    write_tsv(out_dir / (args.sample + ".fusion_loci.final.tsv"), loci, lcolumns)
    write_tsv(
        out_dir / (args.sample + ".fusion_loci.main.tsv"),
        (row for row in loci if row["result_tier"] == "MAIN"), lcolumns,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_loci.candidate.tsv"),
        (row for row in loci if row["result_tier"] == "CANDIDATE"), lcolumns,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_loci.review.tsv"),
        (row for row in loci if row["result_tier"].startswith("REVIEW_")), lcolumns,
    )

    write_tsv(
        out_dir / (args.sample + ".fusion_umi_families.exact.tsv.gz"), exact_umi_rows,
        ["sample", "event_id", "tube", "umi", "family_read_count",
         "supporting_sequence_count", "evidence_sources", "corrected_umi"], gzip_output=True,
    )
    # Compatibility alias: this remains the exact, uncorrected UMI table used for tiering.
    write_tsv(
        out_dir / (args.sample + ".fusion_umi_families.final.tsv.gz"), exact_umi_rows,
        ["sample", "event_id", "tube", "umi", "family_read_count",
         "supporting_sequence_count", "evidence_sources", "corrected_umi"], gzip_output=True,
    )
    write_tsv(
        out_dir / (args.sample + ".fusion_umi_families.corrected.tsv.gz"), corrected_umi_rows,
        ["sample", "event_id", "tube", "corrected_umi", "member_exact_umi_count",
         "member_exact_umis", "family_read_count", "supporting_sequence_count",
         "evidence_sources"], gzip_output=True,
    )

    reason_groups: DefaultDict[str, dict] = defaultdict(
        lambda: {"rows": 0, "raw_reads": 0, "sequences": set(), "events": set()}
    )
    for row in strict_rows:
        reason = rejection_reason(row, args.primary_threshold, args)
        group = reason_groups[reason]
        group["rows"] += 1
        group["raw_reads"] += as_int(row, "raw_read_count")
        group["sequences"].add(row["sequence_id"])
        group["events"].add(row["event_id"])
    filter_rows = [{
        "filter_class": reason, "evidence_rows": values["rows"],
        "unique_sequences": len(values["sequences"]), "events": len(values["events"]),
        "raw_read_sum": values["raw_reads"],
    } for reason, values in sorted(reason_groups.items())]
    write_tsv(
        out_dir / (args.sample + ".fusion_filter_summary.tsv"), filter_rows,
        ["filter_class", "evidence_rows", "unique_sequences", "events", "raw_read_sum"],
    )

    summary = {
        "sample": args.sample, "pipeline_version": args.pipeline_version,
        "primary_threshold": args.primary_threshold, "evaluated_thresholds": thresholds,
        "final_report_thresholds": report_thresholds,
        "strict_audit_evidence_rows": len(strict_rows),
        "rescue_audit_evidence_rows": len(rescue_rows),
        "final_evidence_rows": len(evidence), "final_unique_sequences": len(sequence_rows),
        "original_sequences_represented": len({row["original_sequence_id"] for row in member_rows}),
        "normalized_evidence_source_counts": dict(sorted(Counter(
            row["evidence_sources"] for row in evidence
        ).items())),
        "normalized_sequence_source_counts": dict(sorted(Counter(
            row["evidence_sources"] for row in sequence_rows
        ).items())),
        "normalization_type_counts": dict(sorted(Counter(
            row["normalization_type"] for row in member_rows
        ).items())),
        "final_exact_events": len(events), "final_loci": len(loci),
        "event_tier_counts": dict(sorted(Counter(row["result_tier"] for row in events).items())),
        "event_support_band_counts": dict(sorted(Counter(row["support_band"] for row in events).items())),
        "event_report_class_counts": dict(sorted(Counter(row["report_class"] for row in events).items())),
        "locus_tier_counts": dict(sorted(Counter(row["result_tier"] for row in loci).items())),
        "exact_umi_families": len(exact_umi_rows),
        "corrected_umi_families": len(corrected_umi_rows),
        "umi_family_definition": "sample+tube+exact_event+exact_umi",
        "umi_correction": {
            "primary_tiering_uses": "exact_umi",
            "method": "directional_count_aware_hamming",
            "max_hamming_distance": args.umi_correction_distance,
            "parent_rule": "parent_count >= 2*child_count - 1",
        },
        "support_thresholds": {
            "main_min_raw_reads": args.main_min_raw_reads,
            "main_min_umi_families": args.main_min_umi_families,
            "candidate_min_raw_reads": args.candidate_min_raw_reads,
            "candidate_min_umi_families": args.candidate_min_umi_families,
        },
        "locus_clustering": {
            "human_window_bases": args.locus_human_window,
            "hpv_window_bases": args.locus_hpv_window,
            "same_reference_strand_breakend_orientation_required": True,
            "exact_event_table_preserved": True,
        },
        "coordinate_conventions": {
            "breakpoint_1based": "1-based inclusive reference base adjacent to the junction",
            "human_breakpoint_start0_end0": "0-based half-open one-base BED interval",
        },
        "audit_sources": {
            "strict_evidence": str(Path(args.strict_evidence).resolve()),
            "rescue_evidence": str(Path(args.rescue_evidence).resolve()) if args.rescue_evidence else "",
        },
    }
    summary_path = out_dir / (args.sample + ".fusion_summary.final.json")
    tmp = Path(str(summary_path) + ".tmp")
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")
    os.replace(tmp, summary_path)
    print("[OK] exact_events=%d loci=%d evidence=%d normalized_sequences=%d" % (
        len(events), len(loci), len(evidence), len(sequence_rows)
    ), flush=True)
    print("[TIER] %s" % summary["event_tier_counts"], flush=True)
    print("[OUT] %s" % (out_dir / (args.sample + ".fusion_events.final.tsv")), flush=True)


if __name__ == "__main__":
    main()
