#!/usr/bin/env python3
"""Shared FASTQ, UMI and SAM/CIGAR helpers for the HPV-fusion pipeline."""

from __future__ import annotations

import gzip
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, TextIO, Tuple


X7_FWD = "GTGCAAG"
X7_RC = "CTTGCAC"
UMI_LEN = 15
UMI_ALLOWED_BASES = frozenset("ACT")  # H15 in the original qBDA code.
UMI_ELIGIBLE_STATUSES = frozenset(("OK_R2_ONLY", "OK_R1_R2_CONCORDANT"))
COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")
CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def open_text_auto(path: str, mode: str = "rt"):
    if path == "-":
        return sys.stdin if "r" in mode else sys.stdout
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode, encoding="utf-8")


def normalize_read_id(header: str) -> str:
    """Normalize Illumina/PEAR FASTQ identifiers without losing internal punctuation."""
    read_id = header.strip().split()[0]
    if read_id.startswith("@"):
        read_id = read_id[1:]
    return re.sub(r"/[12]$", "", read_id)


def iter_fastq(path: str) -> Iterator[Tuple[str, str, str]]:
    with open_text_auto(path, "rt") as handle:
        record_number = 0
        while True:
            header = handle.readline()
            if not header:
                break
            record_number += 1
            seq = handle.readline().rstrip("\r\n")
            plus = handle.readline()
            qual = handle.readline().rstrip("\r\n")
            if not seq or not plus or not qual:
                raise ValueError("Incomplete FASTQ record %d in %s" % (record_number, path))
            if not header.startswith("@") or not plus.startswith("+"):
                raise ValueError("Malformed FASTQ record %d in %s" % (record_number, path))
            if len(seq) != len(qual):
                raise ValueError("Sequence/quality length mismatch in record %d of %s" % (record_number, path))
            yield normalize_read_id(header), seq.upper(), qual


def iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    name: Optional[str] = None
    chunks: List[str] = []
    with open_text_auto(path, "rt") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks).upper()
                name = line[1:].split()[0]
                chunks = []
            else:
                if name is None:
                    raise ValueError("FASTA sequence occurs before a header in %s" % path)
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks).upper()


def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1].upper()


def _find_all(seq: str, pattern: str) -> List[int]:
    positions: List[int] = []
    start = 0
    while True:
        pos = seq.find(pattern, start)
        if pos < 0:
            return positions
        positions.append(pos)
        start = pos + 1


def _is_legal_umi(umi: str) -> bool:
    return len(umi) == UMI_LEN and all(base in UMI_ALLOWED_BASES for base in umi.upper())


def extract_umi_candidates(seq: str, source_prefix: str) -> List[Tuple[str, str]]:
    """Implement the original H15-X7 extraction rule exactly."""
    seq = seq.upper()
    candidates: List[Tuple[str, str]] = []
    for pos in _find_all(seq, X7_FWD):
        if pos >= UMI_LEN:
            umi = seq[pos - UMI_LEN:pos]
            if _is_legal_umi(umi):
                candidates.append((umi, "%s:X7_FWD:before" % source_prefix))
    for pos in _find_all(seq, X7_RC):
        end = pos + len(X7_RC)
        if end + UMI_LEN <= len(seq):
            umi = revcomp(seq[end:end + UMI_LEN])
            if _is_legal_umi(umi):
                candidates.append((umi, "%s:X7_RC:after_rc" % source_prefix))
    return candidates


def detect_technical_intervals(seq: str) -> List[Tuple[int, int]]:
    """Locate complete H15-X7-N11 cassettes as 0-based half-open query intervals."""
    seq = seq.upper()
    intervals: List[Tuple[int, int]] = []
    for pos in _find_all(seq, X7_FWD):
        end = pos + len(X7_FWD)
        if pos >= UMI_LEN and end + 11 <= len(seq) and _is_legal_umi(seq[pos - UMI_LEN:pos]):
            intervals.append((pos - UMI_LEN, end + 11))
    for pos in _find_all(seq, X7_RC):
        end = pos + len(X7_RC)
        if pos >= 11 and end + UMI_LEN <= len(seq):
            if _is_legal_umi(revcomp(seq[end:end + UMI_LEN])):
                intervals.append((pos - 11, end + UMI_LEN))
    if not intervals:
        return []
    merged: List[Tuple[int, int]] = []
    for start0, end0 in sorted(intervals):
        if merged and start0 <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end0))
        else:
            merged.append((start0, end0))
    return merged


def interval_overlap(start0: int, end0: int, intervals: Sequence[Tuple[int, int]]) -> int:
    return sum(max(0, min(end0, right) - max(start0, left)) for left, right in intervals)


def classify_pair_umi(seq1: str, seq2: str) -> Tuple[str, str, str]:
    r1_candidates = extract_umi_candidates(seq1, "R1")
    r2_candidates = extract_umi_candidates(seq2, "R2")
    r1 = sorted(set(item[0] for item in r1_candidates))
    r2 = sorted(set(item[0] for item in r2_candidates))

    def sources(umi: str, candidates: Sequence[Tuple[str, str]]) -> str:
        return ";".join(sorted(set(src for value, src in candidates if value == umi)))

    if len(r2) == 1 and len(r1) == 0:
        return r2[0], "OK_R2_ONLY", sources(r2[0], r2_candidates)
    if len(r2) == 1 and len(r1) == 1:
        if r2[0] == r1[0]:
            return r2[0], "OK_R1_R2_CONCORDANT", ";".join(
                (sources(r2[0], r2_candidates), sources(r1[0], r1_candidates))
            )
        return ";".join(sorted(set(r1 + r2))), "DISCORDANT", ";".join(
            src for _, src in r2_candidates + r1_candidates
        )
    if len(r2) == 0 and len(r1) == 1:
        return r1[0], "R1_ONLY_RESCUE", sources(r1[0], r1_candidates)
    if len(r2) == 0 and len(r1) == 0:
        return "", "NOT_FOUND", ""
    if len(r2) > 1 and len(r1) <= 1:
        return ";".join(r2), "AMBIGUOUS_R2", ";".join(src for _, src in r2_candidates)
    if len(r1) > 1 and len(r2) <= 1:
        return ";".join(r1), "AMBIGUOUS_R1", ";".join(src for _, src in r1_candidates)
    return ";".join(sorted(set(r1 + r2))), "AMBIGUOUS_BOTH", ";".join(
        src for _, src in r2_candidates + r1_candidates
    )


def parse_cigar(cigar: str) -> List[Tuple[int, str]]:
    parts = [(int(length), op) for length, op in CIGAR_RE.findall(cigar)]
    if not parts or "".join("%d%s" % item for item in parts) != cigar:
        raise ValueError("Unsupported or malformed CIGAR: %s" % cigar)
    return parts


def _outer_clip(parts: Sequence[Tuple[int, str]], from_left: bool) -> int:
    iterable = parts if from_left else reversed(parts)
    total = 0
    for length, op in iterable:
        if op in ("S", "H"):
            total += length
        else:
            break
    return total


@dataclass(frozen=True)
class CigarMetrics:
    query_start0: int
    query_end0: int
    match_bases: int
    aligned_query_bases: int
    ref_bases: int
    left_clip: int
    right_clip: int


def cigar_metrics(cigar: str, query_length: int, reverse: bool) -> CigarMetrics:
    """Return half-open query coordinates in the original input-read orientation."""
    parts = parse_cigar(cigar)
    left_aligned_clip = _outer_clip(parts, True)
    right_aligned_clip = _outer_clip(parts, False)
    if left_aligned_clip + right_aligned_clip > query_length:
        raise ValueError("CIGAR clipping exceeds query length: %s, qlen=%d" % (cigar, query_length))

    if reverse:
        query_start0 = right_aligned_clip
        query_end0 = query_length - left_aligned_clip
    else:
        query_start0 = left_aligned_clip
        query_end0 = query_length - right_aligned_clip

    match_bases = sum(length for length, op in parts if op in ("M", "=", "X"))
    aligned_query_bases = sum(length for length, op in parts if op in ("M", "I", "=", "X"))
    ref_bases = sum(length for length, op in parts if op in ("M", "D", "N", "=", "X"))
    return CigarMetrics(
        query_start0=query_start0,
        query_end0=query_end0,
        match_bases=match_bases,
        aligned_query_bases=aligned_query_bases,
        ref_bases=ref_bases,
        left_clip=query_start0,
        right_clip=query_length - query_end0,
    )


def reference_boundary(ref_start1: int, ref_end1: int, reverse: bool, query_side: str) -> int:
    """Reference base at the start/end boundary of an alignment in original-query orientation."""
    if query_side not in ("start", "end"):
        raise ValueError("query_side must be start or end")
    if query_side == "start":
        return ref_end1 if reverse else ref_start1
    return ref_start1 if reverse else ref_end1


def parse_sam_tags(fields: Sequence[str]) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for field in fields[11:]:
        pieces = field.split(":", 2)
        if len(pieces) == 3:
            tags[pieces[0]] = pieces[2]
    return tags


def wrap_fasta(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def ensure_parent(path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
