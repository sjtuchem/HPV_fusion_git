#!/usr/bin/env python3
"""Call corrected HPV-human breakpoints and aggregate reads, sequences and UMI families."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Sequence, Set, Tuple

from fusion_common import cigar_metrics, interval_overlap, parse_sam_tags, reference_boundary


# A highly recurrent sequence can carry many distinct UMI families.  Their compact
# semicolon-delimited representation may exceed Python's conservative 128 KiB CSV
# default even though the TSV record itself is valid.
csv.field_size_limit(1024 * 1024 * 1024)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--human-sam", default="-", help="STAR SAM or - for stdin")
    parser.add_argument("--hpv-alignments", required=True)
    parser.add_argument("--sequence-metadata", required=True)
    parser.add_argument("--hpv-type-map", required=True)
    parser.add_argument("--gene-bed", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--thresholds", default="20,30,40,50")
    parser.add_argument("--min-human-mapq", type=int, default=20)
    parser.add_argument("--min-hpv-mapq", type=int, default=0)
    parser.add_argument("--max-junction-gap", type=int, default=20)
    parser.add_argument("--full-human-fraction", type=float, default=0.98)
    parser.add_argument(
        "--audit-only",
        action="store_true",
        help="Write the complete evidence audit and summary only; final compact tables are made later.",
    )
    return parser.parse_args()


def read_hpv_alignments(path: str) -> Dict[str, dict]:
    rows: Dict[str, dict] = {}
    integer_fields = {
        "query_length", "hpv_ref_length", "flag", "mapq", "hpv_match_bases",
        "hpv_aligned_query_bases", "hpv_ref_start1", "hpv_ref_end1",
        "hpv_query_start0", "hpv_query_end0", "left_flank_bases", "right_flank_bases",
        "edit_distance", "alignment_score",
    }
    with open(path, "r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for key in integer_fields:
                row[key] = int(row[key])
            rows[row["sequence_id"]] = row
    return rows


def parse_umi_family_counts(value: str) -> Dict[str, int]:
    counts: Dict[str, int] = {}
    if not value or value == ".":
        return counts
    for item in value.split(";"):
        family, count = item.rsplit(":", 1)
        counts[family] = int(count)
    return counts


def read_sequence_metadata(path: str, needed: Set[str]) -> Dict[str, dict]:
    rows: Dict[str, dict] = {}
    if not needed:
        return rows
    with gzip.open(path, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            sequence_id = row["sequence_id"]
            if sequence_id not in needed:
                continue
            for key in (
                "sequence_length", "raw_read_count", "umi_eligible_read_count",
                "unique_umi_count", "unique_umi_family_count", "tube_A_read_count",
                "tube_B_read_count", "unknown_tube_read_count", "missing_raw_id_count",
                "ambiguous_read_id_count", "technical_bases",
            ):
                row[key] = int(row[key])
            row["family_counts"] = parse_umi_family_counts(row["umi_family_counts"])
            intervals = []
            if row["technical_intervals0"] and row["technical_intervals0"] != ".":
                for item in row["technical_intervals0"].split(";"):
                    start0, end0 = item.split("-", 1)
                    intervals.append((int(start0), int(end0)))
            row["technical_intervals"] = intervals
            rows[sequence_id] = row
    missing = sorted(needed - set(rows))
    if missing:
        raise ValueError("Missing candidate sequence metadata; first=%s" % missing[0])
    return rows


def load_hpv_types(path: str) -> Dict[str, str]:
    mapping: Dict[str, str] = {}
    with open(path, "r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"GenBank ID", "Virus name"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError("HPV type CSV must contain GenBank ID and Virus name")
        for row in reader:
            accession = row["GenBank ID"].strip()
            name = row["Virus name"].strip()
            if accession:
                mapping[accession] = name
                mapping[accession.split(".")[0]] = name
    return mapping


class GenePointIndex:
    """Small binned interval index for a four-column BED: chrom, start0, end0, gene."""

    def __init__(self, path: str, bin_size: int = 1000000):
        self.bin_size = bin_size
        self.bins: DefaultDict[Tuple[str, int], List[Tuple[int, int, str]]] = defaultdict(list)
        self.interval_count = 0
        with open(path, "r", encoding="utf-8") as handle:
            for line_number, raw_line in enumerate(handle, start=1):
                if not raw_line.strip() or raw_line.startswith("#"):
                    continue
                fields = raw_line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                try:
                    start0 = int(fields[1])
                    end0 = int(fields[2])
                except ValueError:
                    continue
                if start0 < 0 or end0 <= start0:
                    raise ValueError("Invalid gene BED interval at line %d" % line_number)
                chrom, gene = fields[0], fields[3]
                for bin_number in range(start0 // bin_size, (end0 - 1) // bin_size + 1):
                    self.bins[(chrom, bin_number)].append((start0, end0, gene))
                self.interval_count += 1

    def genes_at(self, chrom: str, position1: int) -> List[str]:
        position0 = position1 - 1
        chroms = [chrom]
        if chrom.startswith("chr"):
            chroms.append(chrom[3:])
        else:
            chroms.append("chr" + chrom)
        genes: Set[str] = set()
        for candidate_chrom in chroms:
            for start0, end0, gene in self.bins.get((candidate_chrom, position0 // self.bin_size), []):
                if start0 <= position0 < end0 and gene:
                    genes.add(gene)
        return sorted(genes)


def query_name_parts(qname: str) -> Tuple[str, str, int, int]:
    fields = qname.rsplit("|", 3)
    if len(fields) != 4 or fields[1] not in ("FULL", "L", "R"):
        raise ValueError("Unexpected STAR query name: %s" % qname)
    return fields[0], fields[1], int(fields[2]), int(fields[3])


def alignment_rank(row: dict) -> Tuple[int, int, int, int]:
    return (row["alignment_score"], row["mapq"], row["human_match_bases"], -row["edit_distance"])


def read_human_alignments(path: str, needed: Set[str]) -> Tuple[Dict[str, dict], Dict[str, List[dict]], dict]:
    full_best: Dict[str, dict] = {}
    flanks: DefaultDict[str, List[dict]] = defaultdict(list)
    stats = {"sam_records": 0, "primary_mapped": 0, "full_mapped": 0, "flank_mapped": 0}
    handle = sys.stdin if path == "-" else open(path, "r", encoding="utf-8")
    try:
        for raw_line in handle:
            if raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("Malformed human SAM record")
            stats["sam_records"] += 1
            flag = int(fields[1])
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            sequence_id, kind, fragment_start0, fragment_end0 = query_name_parts(fields[0])
            if sequence_id not in needed:
                continue
            fragment_length = fragment_end0 - fragment_start0
            if fragment_length <= 0:
                raise ValueError("Invalid STAR fragment coordinates for %s" % fields[0])
            reverse = bool(flag & 0x10)
            metrics = cigar_metrics(fields[5], fragment_length, reverse)
            ref_start1 = int(fields[3])
            ref_end1 = ref_start1 + metrics.ref_bases - 1
            tags = parse_sam_tags(fields)
            row = {
                "sequence_id": sequence_id,
                "kind": kind,
                "fragment_start0": fragment_start0,
                "fragment_end0": fragment_end0,
                "human_chr": fields[2],
                "flag": flag,
                "human_strand": "-" if reverse else "+",
                "reverse": reverse,
                "mapq": int(fields[4]),
                "cigar": fields[5],
                "human_match_bases": metrics.match_bases,
                "human_aligned_query_bases": metrics.aligned_query_bases,
                "human_ref_start1": ref_start1,
                "human_ref_end1": ref_end1,
                "human_query_start0": fragment_start0 + metrics.query_start0,
                "human_query_end0": fragment_start0 + metrics.query_end0,
                "nh": int(tags.get("NH", "1")),
                "edit_distance": int(tags.get("NM", "0")),
                "alignment_score": int(tags.get("AS", str(metrics.match_bases))),
            }
            stats["primary_mapped"] += 1
            if kind == "FULL":
                stats["full_mapped"] += 1
                previous = full_best.get(sequence_id)
                if previous is None or alignment_rank(row) > alignment_rank(previous):
                    full_best[sequence_id] = row
            else:
                stats["flank_mapped"] += 1
                flanks[sequence_id].append(row)
    finally:
        if handle is not sys.stdin:
            handle.close()
    return full_best, dict(flanks), stats


def normalized_hpv_position(position1: int, reference_length: int) -> int:
    if reference_length <= 0:
        return position1
    return ((position1 - 1) % reference_length) + 1


def event_identifier(key: Sequence[object]) -> str:
    encoded = "\t".join(str(value) for value in key).encode("utf-8")
    return "EVT_" + hashlib.sha256(encoded).hexdigest()[:20]


def bool_int(value: bool) -> int:
    return 1 if value else 0


def write_outputs(
    args: argparse.Namespace,
    thresholds: List[int],
    hpv_rows: Dict[str, dict],
    metadata: Dict[str, dict],
    type_map: Dict[str, str],
    gene_index: GenePointIndex,
    full_best: Dict[str, dict],
    flank_rows: Dict[str, List[dict]],
    stats: dict,
) -> None:
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    min_threshold = min(thresholds)
    evidence_by_key: Dict[Tuple[object, ...], dict] = {}

    for sequence_id, hpv in hpv_rows.items():
        full = full_best.get(sequence_id)
        full_human_fraction = 0.0
        if full is not None and hpv["query_length"] > 0:
            full_human_fraction = full["human_match_bases"] / float(hpv["query_length"])
        full_human_like = full_human_fraction >= args.full_human_fraction
        hpv_reverse = hpv["hpv_strand"] == "-"

        for human in flank_rows.get(sequence_id, []):
            side = human["kind"]
            if side == "L":
                junction_gap = hpv["hpv_query_start0"] - human["human_query_end0"]
                human_query_side = "end"
                hpv_query_side = "start"
                orientation = "human(L)->HPV"
            elif side == "R":
                junction_gap = human["human_query_start0"] - hpv["hpv_query_end0"]
                human_query_side = "start"
                hpv_query_side = "end"
                orientation = "HPV->human(R)"
            else:
                continue
            if junction_gap < 0 or junction_gap > args.max_junction_gap:
                continue
            if hpv["hpv_match_bases"] < min_threshold or human["human_match_bases"] < min_threshold:
                continue

            human_breakpoint1 = reference_boundary(
                human["human_ref_start1"], human["human_ref_end1"], human["reverse"], human_query_side
            )
            hpv_raw_breakpoint1 = reference_boundary(
                hpv["hpv_ref_start1"], hpv["hpv_ref_end1"], hpv_reverse, hpv_query_side
            )
            hpv_breakpoint1 = normalized_hpv_position(hpv_raw_breakpoint1, hpv["hpv_ref_length"])
            human_ref_side = (
                "ref_start" if human_breakpoint1 == human["human_ref_start1"] else "ref_end"
            )
            hpv_ref_side = "ref_start" if hpv_raw_breakpoint1 == hpv["hpv_ref_start1"] else "ref_end"
            hpv_accession = hpv["hpv_ref"]
            hpv_type = type_map.get(hpv_accession, type_map.get(hpv_accession.split(".")[0], hpv_accession))
            genes = gene_index.genes_at(human["human_chr"], human_breakpoint1)
            event_key = (
                hpv_type, hpv_accession, hpv_breakpoint1, hpv["hpv_strand"], hpv_ref_side,
                human["human_chr"], human_breakpoint1, human["human_strand"], human_ref_side,
                orientation,
            )
            event_id = event_identifier(event_key)
            meta = metadata[sequence_id]
            hpv_technical_overlap = interval_overlap(
                hpv["hpv_query_start0"], hpv["hpv_query_end0"], meta["technical_intervals"]
            )
            human_technical_overlap = interval_overlap(
                human["human_query_start0"], human["human_query_end0"], meta["technical_intervals"]
            )
            hpv_effective_match = max(0, hpv["hpv_match_bases"] - hpv_technical_overlap)
            human_effective_match = max(0, human["human_match_bases"] - human_technical_overlap)
            mapping_confident = (
                not full_human_like
                and human["mapq"] >= args.min_human_mapq
                and human["nh"] == 1
                and hpv["mapq"] >= args.min_hpv_mapq
                and hpv_technical_overlap == 0
                and human_technical_overlap == 0
            )
            row = {
                "sample": args.sample,
                "event_id": event_id,
                "sequence_id": sequence_id,
                "hpv_type": hpv_type,
                "hpv_ref": hpv_accession,
                "hpv_ref_length": hpv["hpv_ref_length"],
                "hpv_breakpoint_1based": hpv_breakpoint1,
                "hpv_strand": hpv["hpv_strand"],
                "hpv_breakend_side": hpv_ref_side,
                "human_chr": human["human_chr"],
                "human_breakpoint_1based": human_breakpoint1,
                "human_breakpoint_start0": human_breakpoint1 - 1,
                "human_breakpoint_end0": human_breakpoint1,
                "human_strand": human["human_strand"],
                "human_breakend_side": human_ref_side,
                "genes": "|".join(genes) if genes else "INTERGENIC",
                "query_orientation": orientation,
                "junction_gap_bases": junction_gap,
                "hpv_match_bases": hpv["hpv_match_bases"],
                "human_match_bases": human["human_match_bases"],
                "hpv_technical_overlap_bases": hpv_technical_overlap,
                "human_technical_overlap_bases": human_technical_overlap,
                "hpv_effective_match_bases": hpv_effective_match,
                "human_effective_match_bases": human_effective_match,
                "hpv_mapq": hpv["mapq"],
                "human_mapq": human["mapq"],
                "human_nh": human["nh"],
                "hpv_cigar": hpv["cigar"],
                "human_cigar": human["cigar"],
                "full_human_aligned_fraction": "%.6f" % full_human_fraction,
                "full_human_like": bool_int(full_human_like),
                "mapping_confident": bool_int(mapping_confident),
                "raw_read_count": meta["raw_read_count"],
                "umi_eligible_read_count": meta["umi_eligible_read_count"],
                "unique_umi_family_count": meta["unique_umi_family_count"],
            }
            for threshold in thresholds:
                length_pass = (
                    hpv_effective_match >= threshold and human_effective_match >= threshold
                )
                row["pass_length_%dbp" % threshold] = bool_int(length_pass)
                row["pass_highconf_%dbp" % threshold] = bool_int(length_pass and mapping_confident)

            evidence_key = (event_id, sequence_id)
            previous = evidence_by_key.get(evidence_key)
            row_rank = (row["mapping_confident"], human["alignment_score"], human["mapq"], human["human_match_bases"])
            if previous is None or row_rank > previous["_rank"]:
                row["_rank"] = row_rank
                evidence_by_key[evidence_key] = row

    evidence_rows = sorted(evidence_by_key.values(), key=lambda row: (row["event_id"], row["sequence_id"]))
    for row in evidence_rows:
        row.pop("_rank", None)

    base_evidence_columns = [
        "sample", "event_id", "sequence_id", "hpv_type", "hpv_ref", "hpv_ref_length",
        "hpv_breakpoint_1based", "hpv_strand", "hpv_breakend_side", "human_chr",
        "human_breakpoint_1based", "human_breakpoint_start0", "human_breakpoint_end0",
        "human_strand", "human_breakend_side", "genes", "query_orientation",
        "junction_gap_bases", "hpv_match_bases", "human_match_bases", "hpv_mapq",
        "hpv_technical_overlap_bases", "human_technical_overlap_bases",
        "hpv_effective_match_bases", "human_effective_match_bases",
        "human_mapq", "human_nh", "hpv_cigar", "human_cigar",
        "full_human_aligned_fraction", "full_human_like", "mapping_confident",
        "raw_read_count", "umi_eligible_read_count", "unique_umi_family_count",
    ]
    evidence_columns = base_evidence_columns + [
        name for threshold in thresholds
        for name in ("pass_length_%dbp" % threshold, "pass_highconf_%dbp" % threshold)
    ]
    evidence_path = out_dir / (args.sample + ".fusion_evidence.tsv.gz")
    with gzip.open(evidence_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=evidence_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(evidence_rows)

    if args.audit_only:
        event_ids = set(row["event_id"] for row in evidence_rows)
        sequence_ids = set(row["sequence_id"] for row in evidence_rows)
        stats.update({
            "candidate_sequences": len(hpv_rows),
            "fusion_evidence_rows": len(evidence_rows),
            "fusion_events_in_audit": len(event_ids),
            "fusion_unique_sequences_in_audit": len(sequence_ids),
            "gene_intervals_loaded": gene_index.interval_count,
            "thresholds": thresholds,
            "audit_only": True,
            "coordinate_conventions": {
                "breakpoint_1based": "1-based inclusive reference base adjacent to the junction",
                "human_breakpoint_start0_end0": "0-based half-open one-base BED interval",
                "query_intervals": "0-based half-open in original merged-read orientation",
            },
        })
        summary_path = out_dir / (args.sample + ".fusion_summary.json")
        with open(summary_path, "w", encoding="utf-8") as handle:
            json.dump(stats, handle, ensure_ascii=False, indent=2, sort_keys=True)
            handle.write("\n")
        print("[OK] audit evidence=%d events=%d unique_sequences=%d" % (
            len(evidence_rows), len(event_ids), len(sequence_ids)
        ), flush=True)
        print("[OUT] %s" % evidence_path, flush=True)
        print("[OUT] %s" % summary_path, flush=True)
        return

    grouped: DefaultDict[str, List[dict]] = defaultdict(list)
    for row in evidence_rows:
        grouped[row["event_id"]].append(row)
    event_rows: List[dict] = []
    for event_id, rows in grouped.items():
        first = rows[0]
        event = dict((key, first[key]) for key in (
            "sample", "event_id", "hpv_type", "hpv_ref", "hpv_ref_length",
            "hpv_breakpoint_1based", "hpv_strand", "hpv_breakend_side", "human_chr",
            "human_breakpoint_1based", "human_breakpoint_start0", "human_breakpoint_end0",
            "human_strand", "human_breakend_side", "query_orientation",
        ))
        all_genes: Set[str] = set()
        for row in rows:
            all_genes.update(row["genes"].split("|"))
        if len(all_genes) > 1:
            all_genes.discard("INTERGENIC")
        event["genes"] = "|".join(sorted(all_genes))
        event["max_junction_gap_bases"] = max(row["junction_gap_bases"] for row in rows)

        for threshold in thresholds:
            passing = [row for row in rows if row["pass_length_%dbp" % threshold]]
            highconf = [row for row in rows if row["pass_highconf_%dbp" % threshold]]
            families: Set[str] = set()
            highconf_families: Set[str] = set()
            for row in passing:
                families.update(metadata[row["sequence_id"]]["family_counts"])
            for row in highconf:
                highconf_families.update(metadata[row["sequence_id"]]["family_counts"])
            event["unique_sequences_%dbp" % threshold] = len(passing)
            event["raw_reads_%dbp" % threshold] = sum(row["raw_read_count"] for row in passing)
            event["umi_eligible_reads_%dbp" % threshold] = sum(
                row["umi_eligible_read_count"] for row in passing
            )
            event["unique_umi_families_%dbp" % threshold] = len(families)
            event["highconf_unique_sequences_%dbp" % threshold] = len(highconf)
            event["highconf_raw_reads_%dbp" % threshold] = sum(row["raw_read_count"] for row in highconf)
            event["highconf_unique_umi_families_%dbp" % threshold] = len(highconf_families)
        event["recommended_pass_30bp"] = bool_int(event.get("highconf_unique_sequences_30bp", 0) > 0)
        event_rows.append(event)

    event_rows.sort(key=lambda row: (
        -row.get("highconf_raw_reads_30bp", 0), row["human_chr"], row["human_breakpoint_1based"],
        row["hpv_type"], row["hpv_breakpoint_1based"],
    ))
    event_base_columns = [
        "sample", "event_id", "hpv_type", "hpv_ref", "hpv_ref_length",
        "hpv_breakpoint_1based", "hpv_strand", "hpv_breakend_side", "human_chr",
        "human_breakpoint_1based", "human_breakpoint_start0", "human_breakpoint_end0",
        "human_strand", "human_breakend_side", "genes", "query_orientation",
        "max_junction_gap_bases",
    ]
    event_count_columns = [
        name for threshold in thresholds
        for name in (
            "unique_sequences_%dbp" % threshold, "raw_reads_%dbp" % threshold,
            "umi_eligible_reads_%dbp" % threshold, "unique_umi_families_%dbp" % threshold,
            "highconf_unique_sequences_%dbp" % threshold, "highconf_raw_reads_%dbp" % threshold,
            "highconf_unique_umi_families_%dbp" % threshold,
        )
    ]
    event_columns = event_base_columns + event_count_columns + ["recommended_pass_30bp"]
    events_path = out_dir / (args.sample + ".fusion_events.tsv")
    with open(events_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=event_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(event_rows)

    sequence_to_events: DefaultDict[str, Set[str]] = defaultdict(set)
    for row in evidence_rows:
        sequence_to_events[row["sequence_id"]].add(row["event_id"])
    sequence_columns = [
        "sample", "sequence_id", "event_ids", "event_count", "sequence_length",
        "raw_read_count", "umi_eligible_read_count", "unique_umi_count",
        "unique_umi_family_count", "umi_family_counts", "tube_A_read_count",
        "tube_B_read_count", "unknown_tube_read_count", "missing_raw_id_count",
        "ambiguous_read_id_count", "sequence",
        "technical_intervals0", "technical_bases",
    ]
    sequences_path = out_dir / (args.sample + ".fusion_unique_sequences.tsv.gz")
    with gzip.open(sequences_path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=sequence_columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sequence_id in sorted(sequence_to_events):
            meta = metadata[sequence_id]
            out = dict((key, meta[key]) for key in sequence_columns if key not in ("event_ids", "event_count"))
            out["event_ids"] = ";".join(sorted(sequence_to_events[sequence_id]))
            out["event_count"] = len(sequence_to_events[sequence_id])
            writer.writerow(out)

    stats.update({
        "candidate_sequences": len(hpv_rows),
        "fusion_evidence_rows": len(evidence_rows),
        "fusion_events": len(event_rows),
        "fusion_unique_sequences": len(sequence_to_events),
        "gene_intervals_loaded": gene_index.interval_count,
        "thresholds": thresholds,
        "coordinate_conventions": {
            "breakpoint_1based": "1-based inclusive reference base adjacent to the junction",
            "human_breakpoint_start0_end0": "0-based half-open one-base BED interval",
            "query_intervals": "0-based half-open in original merged-read orientation",
        },
    })
    summary_path = out_dir / (args.sample + ".fusion_summary.json")
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(stats, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    print("[OK] fusion evidence=%d events=%d unique_sequences=%d" % (
        len(evidence_rows), len(event_rows), len(sequence_to_events)
    ), flush=True)
    print("[OUT] %s" % events_path, flush=True)
    print("[OUT] %s" % sequences_path, flush=True)


def main() -> None:
    args = parse_args()
    thresholds = sorted(set(int(value) for value in args.thresholds.split(",") if value.strip()))
    if not thresholds or any(value <= 0 for value in thresholds):
        raise ValueError("Thresholds must be positive integers")
    if 30 not in thresholds:
        raise ValueError("Threshold list must include 30 for recommended_pass_30bp")
    hpv_rows = read_hpv_alignments(args.hpv_alignments)
    metadata = read_sequence_metadata(args.sequence_metadata, set(hpv_rows))
    type_map = load_hpv_types(args.hpv_type_map)
    gene_index = GenePointIndex(args.gene_bed)
    full_best, flank_rows, stats = read_human_alignments(args.human_sam, set(hpv_rows))
    write_outputs(
        args, thresholds, hpv_rows, metadata, type_map, gene_index,
        full_best, flank_rows, stats,
    )


if __name__ == "__main__":
    main()
