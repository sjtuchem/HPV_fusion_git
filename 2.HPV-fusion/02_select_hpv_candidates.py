#!/usr/bin/env python3
"""Parse one BWA-to-HPV stream and create a single set of STAR queries."""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys
from pathlib import Path
from typing import Dict, Optional, Tuple

from fusion_common import cigar_metrics, iter_fasta, parse_sam_tags, wrap_fasta


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sam", default="-", help="BWA SAM, or - for stdin")
    parser.add_argument("--rotated-first", action="store_true",
                        help="Read midpoint-rotated SAM first, then original SAM after the stream delimiter")
    parser.add_argument("--unique-fasta", required=True)
    parser.add_argument("--out-alignments", required=True)
    parser.add_argument("--out-star-fasta", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--prefilter-min-hpv", type=int, default=15)
    parser.add_argument("--prefilter-min-flank", type=int, default=15)
    return parser.parse_args()


def candidate_rank(row: dict) -> Tuple[int, int, int, int]:
    return (row["alignment_score"], row["mapq"], row["hpv_match_bases"], -row["edit_distance"])


def main() -> None:
    args = parse_args()
    for output in (args.out_alignments, args.out_star_fasta, args.summary):
        Path(output).parent.mkdir(parents=True, exist_ok=True)

    sam_handle = sys.stdin if args.sam == "-" else open(args.sam, "r", encoding="utf-8")
    ref_lengths: Dict[str, int] = {}
    candidates: Dict[str, dict] = {}
    rotated_candidates: Dict[str, dict] = {}
    rotated_scores: Dict[str, int] = {}
    rotated = args.rotated_first
    counts = {
        "sam_records": 0,
        "primary_mapped_records": 0,
        "partial_prefilter_records": 0,
        "candidate_sequences": 0,
        "star_full_queries": 0,
        "star_left_flank_queries": 0,
        "star_right_flank_queries": 0,
    }

    try:
        for raw_line in sam_handle:
            if raw_line.rstrip("\r\n") == "@CO\tHPV_ORIGINAL_REFERENCE":
                rotated = False
                continue
            if raw_line.startswith("@SQ"):
                fields = raw_line.rstrip("\n").split("\t")
                values = dict(item.split(":", 1) for item in fields[1:] if ":" in item)
                if "SN" in values and "LN" in values:
                    ref_lengths[values["SN"]] = int(values["LN"])
                continue
            if raw_line.startswith("@"):
                continue
            fields = raw_line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError("Malformed SAM record with fewer than 11 fields")
            counts["sam_records"] += 1
            sequence_id = fields[0]
            flag = int(fields[1])
            if flag & 0x4 or flag & 0x100 or flag & 0x800:
                continue
            counts["primary_mapped_records"] += 1
            hpv_ref = fields[2]
            ref_start1 = int(fields[3])
            mapq = int(fields[4])
            cigar = fields[5]
            query_length = len(fields[9])
            reverse = bool(flag & 0x10)
            metrics = cigar_metrics(cigar, query_length, reverse)
            ref_length = ref_lengths.get(hpv_ref, 0)
            ref_end1 = ref_start1 + metrics.ref_bases - 1
            tags = parse_sam_tags(fields)
            score = int(tags.get("AS", str(metrics.match_bases)))
            if rotated:
                midpoint = ref_length // 2
                origin = ref_length - midpoint
                if not ref_start1 <= origin < ref_end1:
                    continue
                rotated_scores[sequence_id] = max(score, rotated_scores.get(sequence_id, score))
                # Keep an unwrapped interval so the caller can select the correct
                # end by strand, then fold that breakpoint into 1..ref_length.
                ref_start1 += midpoint
                ref_end1 += midpoint
            else:
                alternative_score = rotated_scores.get(sequence_id)
                if alternative_score is not None:
                    if score < alternative_score:
                        continue
                    rotated_candidates.pop(sequence_id, None)
            left_flank = metrics.query_start0
            right_flank = query_length - metrics.query_end0
            if metrics.match_bases < args.prefilter_min_hpv:
                continue
            if max(left_flank, right_flank) < args.prefilter_min_flank:
                continue
            counts["partial_prefilter_records"] += 1
            row = {
                "sequence_id": sequence_id,
                "query_length": query_length,
                "hpv_ref": hpv_ref,
                "hpv_ref_length": ref_length,
                "flag": flag,
                "hpv_strand": "-" if reverse else "+",
                "mapq": mapq,
                "cigar": cigar,
                "hpv_match_bases": metrics.match_bases,
                "hpv_aligned_query_bases": metrics.aligned_query_bases,
                "hpv_ref_start1": ref_start1,
                "hpv_ref_end1": ref_end1,
                "hpv_query_start0": metrics.query_start0,
                "hpv_query_end0": metrics.query_end0,
                "left_flank_bases": left_flank,
                "right_flank_bases": right_flank,
                "edit_distance": int(tags.get("NM", "0")),
                "alignment_score": score,
            }
            destination = rotated_candidates if rotated else candidates
            previous = destination.get(sequence_id)
            if previous is None or candidate_rank(row) > candidate_rank(previous):
                destination[sequence_id] = row
    finally:
        if sam_handle is not sys.stdin:
            sam_handle.close()

    if rotated:
        raise ValueError("Missing original-reference SAM stream")
    candidates.update(rotated_candidates)
    counts["circular_origin_candidates"] = len(rotated_candidates)
    counts["candidate_sequences"] = len(candidates)
    align_tmp = args.out_alignments + ".tmp"
    fasta_tmp = args.out_star_fasta + ".tmp"
    summary_tmp = args.summary + ".tmp"
    columns = [
        "sequence_id", "query_length", "hpv_ref", "hpv_ref_length", "flag", "hpv_strand",
        "mapq", "cigar", "hpv_match_bases", "hpv_aligned_query_bases",
        "hpv_ref_start1", "hpv_ref_end1", "hpv_query_start0", "hpv_query_end0",
        "left_flank_bases", "right_flank_bases", "edit_distance", "alignment_score",
    ]
    with open(align_tmp, "w", encoding="utf-8", newline="") as align_out:
        writer = csv.DictWriter(align_out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for sequence_id in sorted(candidates):
            writer.writerow(candidates[sequence_id])

    found = set()
    with open(fasta_tmp, "w", encoding="ascii") as fasta_out:
        for sequence_id, seq in iter_fasta(args.unique_fasta):
            row = candidates.get(sequence_id)
            if row is None:
                continue
            found.add(sequence_id)
            query_length = len(seq)
            if query_length != row["query_length"]:
                raise ValueError("Query length mismatch for %s" % sequence_id)

            fasta_out.write(">%s|FULL|0|%d\n%s\n" % (sequence_id, query_length, wrap_fasta(seq)))
            counts["star_full_queries"] += 1
            left_end = row["hpv_query_start0"]
            if left_end >= args.prefilter_min_flank:
                fasta_out.write(">%s|L|0|%d\n%s\n" % (
                    sequence_id, left_end, wrap_fasta(seq[:left_end])
                ))
                counts["star_left_flank_queries"] += 1
            right_start = row["hpv_query_end0"]
            if query_length - right_start >= args.prefilter_min_flank:
                fasta_out.write(">%s|R|%d|%d\n%s\n" % (
                    sequence_id, right_start, query_length, wrap_fasta(seq[right_start:])
                ))
                counts["star_right_flank_queries"] += 1

    missing = sorted(set(candidates) - found)
    if missing:
        raise ValueError("%d candidate IDs were absent from unique FASTA; first=%s" % (len(missing), missing[0]))

    with open(summary_tmp, "w", encoding="utf-8") as summary_out:
        json.dump(counts, summary_out, indent=2, sort_keys=True)
        summary_out.write("\n")
    os.replace(align_tmp, args.out_alignments)
    os.replace(fasta_tmp, args.out_star_fasta)
    os.replace(summary_tmp, args.summary)
    print("[OK] HPV candidate unique sequences: %d" % len(candidates), flush=True)


if __name__ == "__main__":
    main()
