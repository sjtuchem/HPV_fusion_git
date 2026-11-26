#!/usr/bin/env python3
"""Extract paired-read UMIs and collapse the merged FASTQ to unique sequences."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sqlite3
from collections import Counter
from itertools import zip_longest
from pathlib import Path
from typing import Dict, List, Tuple

from fusion_common import (
    UMI_ELIGIBLE_STATUSES,
    classify_pair_umi,
    detect_technical_intervals,
    iter_fastq,
    wrap_fasta,
)


STATUS_CODE = {
    "OK_R2_ONLY": "E2",
    "OK_R1_R2_CONCORDANT": "EC",
    "R1_ONLY_RESCUE": "R1",
    "DISCORDANT": "D",
    "NOT_FOUND": "N",
    "AMBIGUOUS_R2": "A2",
    "AMBIGUOUS_R1": "A1",
    "AMBIGUOUS_BOTH": "AB",
}
CODE_STATUS = dict((value, key) for key, value in STATUS_CODE.items())


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--merged-fastq", required=True)
    parser.add_argument("--r1-a", required=True)
    parser.add_argument("--r2-a", required=True)
    parser.add_argument("--r1-b", default="")
    parser.add_argument("--r2-b", default="")
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--keep-umi-index", action="store_true")
    parser.add_argument("--max-missing-raw-id-rate", type=float, default=0.001)
    return parser.parse_args()


def check_inputs(paths: List[str]) -> None:
    for value in paths:
        path = Path(value)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("Missing or empty input: %s" % path)


def configure_sqlite(connection: sqlite3.Connection) -> None:
    connection.execute("PRAGMA journal_mode=OFF")
    connection.execute("PRAGMA synchronous=OFF")
    connection.execute("PRAGMA temp_store=MEMORY")
    connection.execute("PRAGMA locking_mode=EXCLUSIVE")
    connection.execute("PRAGMA cache_size=-262144")  # Up to 256 MiB page cache.


def flush_index_batch(connection: sqlite3.Connection, batch: List[Tuple[str, str, str, str]]) -> None:
    connection.executemany(
        """
        INSERT INTO umi(read_id, tube, status_code, umi) VALUES (?, ?, ?, ?)
        ON CONFLICT(read_id) DO UPDATE SET tube='', status_code='C', umi=''
        """,
        batch,
    )
    connection.commit()
    batch.clear()


def build_umi_index(connection: sqlite3.Connection, pairs: List[Tuple[str, str, str]], stats: Counter) -> None:
    connection.execute(
        "CREATE TABLE umi(read_id TEXT PRIMARY KEY, tube TEXT NOT NULL, status_code TEXT NOT NULL, umi TEXT NOT NULL) WITHOUT ROWID"
    )
    batch: List[Tuple[str, str, str, str]] = []
    for tube, r1_path, r2_path in pairs:
        paired = zip_longest(iter_fastq(r1_path), iter_fastq(r2_path))
        for pair_number, pair in enumerate(paired, start=1):
            read1, read2 = pair
            if read1 is None or read2 is None:
                raise ValueError("R1/R2 counts differ for tube %s at pair %d" % (tube, pair_number))
            rid1, seq1, _ = read1
            rid2, seq2, _ = read2
            if rid1 != rid2:
                raise ValueError(
                    "R1/R2 ID mismatch for tube %s at pair %d: %s != %s" %
                    (tube, pair_number, rid1, rid2)
                )

            stats["raw_pairs_%s" % tube] += 1
            umi, status, _ = classify_pair_umi(seq1, seq2)
            stats["raw_status_%s_%s" % (tube, status)] += 1
            if status in UMI_ELIGIBLE_STATUSES:
                stats["raw_umi_eligible_%s" % tube] += 1
            else:
                umi = ""
            batch.append((rid1, tube, STATUS_CODE[status], umi))
            if len(batch) >= 100000:
                flush_index_batch(connection, batch)
            if pair_number % 1000000 == 0:
                print("[UMI] tube=%s raw_pairs=%d" % (tube, pair_number), flush=True)
    if batch:
        flush_index_batch(connection, batch)
    stats["raw_read_id_collisions"] = connection.execute(
        "SELECT count(*) FROM umi WHERE status_code='C'"
    ).fetchone()[0]


def lookup_umi_batch(connection: sqlite3.Connection, read_ids: List[str]) -> Dict[str, Tuple[str, str, str]]:
    lookup: Dict[str, Tuple[str, str, str]] = {}
    # 500 stays below SQLite's traditional 999-variable limit.
    for start in range(0, len(read_ids), 500):
        chunk = read_ids[start:start + 500]
        placeholders = ",".join("?" for _ in chunk)
        query = "SELECT read_id, tube, status_code, umi FROM umi WHERE read_id IN (%s)" % placeholders
        for read_id, tube, code, umi in connection.execute(query, chunk):
            lookup[read_id] = (tube, CODE_STATUS.get(code, "AMBIGUOUS_READ_ID"), umi)
    return lookup


def add_merged_batch(
    connection: sqlite3.Connection,
    batch: List[Tuple[str, str]],
    sequences: Dict[str, list],
    stats: Counter,
) -> None:
    lookup = lookup_umi_batch(connection, [read_id for read_id, _ in batch])
    for read_id, seq in batch:
        tube = "UNKNOWN"
        status = "MISSING_RAW_ID"
        umi = ""
        value = lookup.get(read_id)
        if value is None:
            stats["merged_missing_raw_id"] += 1
        else:
            tube, status, umi = value
            if status == "AMBIGUOUS_READ_ID":
                tube = "UNKNOWN"
                stats["merged_ambiguous_read_id"] += 1
            else:
                stats["merged_status_%s" % status] += 1

        entry = sequences.get(seq)
        if entry is None:
            entry = [0, 0, Counter(), Counter(), 0, 0]
            sequences[seq] = entry
        entry[0] += 1
        entry[3][tube] += 1
        if status in UMI_ELIGIBLE_STATUSES and umi:
            entry[1] += 1
            entry[2][(tube, umi)] += 1
            stats["merged_umi_eligible_reads"] += 1
        elif status == "MISSING_RAW_ID":
            entry[4] += 1
        elif status == "AMBIGUOUS_READ_ID":
            entry[5] += 1
    batch.clear()


def main() -> None:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if bool(args.r1_b) != bool(args.r2_b):
        raise ValueError("R1-B and R2-B must either both be supplied or both be omitted")
    inputs = [args.merged_fastq, args.r1_a, args.r2_a]
    if args.r1_b:
        inputs.extend((args.r1_b, args.r2_b))
    check_inputs(inputs)

    index_path = out_dir / (args.sample + ".umi_index.sqlite")
    fasta_path = out_dir / (args.sample + ".unique_sequences.fa")
    metadata_path = out_dir / (args.sample + ".unique_sequences.tsv.gz")
    summary_path = out_dir / (args.sample + ".prepare_summary.json")
    fasta_tmp = Path(str(fasta_path) + ".tmp")
    metadata_tmp = Path(str(metadata_path) + ".tmp")
    summary_tmp = Path(str(summary_path) + ".tmp")

    stats: Counter = Counter()
    pairs = [("A", args.r1_a, args.r2_a)]
    if args.r1_b:
        pairs.append(("B", args.r1_b, args.r2_b))
    print("[STEP] Building disk-backed UMI index from raw R1/R2", flush=True)
    if index_path.exists():
        index_path.unlink()
    connection = sqlite3.connect(str(index_path))
    configure_sqlite(connection)
    build_umi_index(connection, pairs, stats)

    # sequence -> [raw reads, eligible reads, family Counter, tube Counter, missing, collision]
    sequences: Dict[str, list] = {}
    print("[STEP] Streaming merged FASTQ and collapsing identical sequences", flush=True)
    merged_batch: List[Tuple[str, str]] = []
    for record_number, (read_id, seq, _) in enumerate(iter_fastq(args.merged_fastq), start=1):
        stats["merged_reads"] += 1
        merged_batch.append((read_id, seq))
        if len(merged_batch) >= 10000:
            add_merged_batch(connection, merged_batch, sequences, stats)
        if record_number % 1000000 == 0:
            print(
                "[COLLAPSE] merged_reads=%d unique_sequences=%d" % (record_number, len(sequences)),
                flush=True,
            )
    if merged_batch:
        add_merged_batch(connection, merged_batch, sequences, stats)
    connection.close()

    raw_pairs_total = sum(stats["raw_pairs_%s" % tube] for tube, _, _ in pairs)
    if stats["merged_reads"] > raw_pairs_total:
        raise ValueError(
            "Merged FASTQ has more records than total raw pairs: %d > %d" %
            (stats["merged_reads"], raw_pairs_total)
        )
    missing_rate = (
        stats["merged_missing_raw_id"] / float(stats["merged_reads"])
        if stats["merged_reads"] else 0.0
    )
    stats["merged_missing_raw_id_rate"] = missing_rate
    if missing_rate > args.max_missing_raw_id_rate:
        raise ValueError(
            "Merged-to-raw read ID missing rate %.6f exceeds allowed %.6f; UMI association is unsafe" %
            (missing_rate, args.max_missing_raw_id_rate)
        )

    stats["unique_sequences"] = len(sequences)
    header = [
        "sample", "sequence_id", "sequence_length", "raw_read_count",
        "umi_eligible_read_count", "unique_umi_count", "unique_umi_family_count",
        "umi_family_counts", "tube_A_read_count", "tube_B_read_count",
        "unknown_tube_read_count", "missing_raw_id_count", "ambiguous_read_id_count",
        "technical_intervals0", "technical_bases", "sequence",
    ]
    with open(fasta_tmp, "w", encoding="ascii") as fasta_out, gzip.open(
        metadata_tmp, "wt", encoding="utf-8", newline=""
    ) as metadata_out:
        writer = csv.writer(metadata_out, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        for seq, entry in sequences.items():
            raw_count, eligible_count, family_counts, tube_counts, missing_count, collision_count = entry
            digest = hashlib.sha256(seq.encode("ascii")).hexdigest()
            sequence_id = "SEQ_" + digest
            umi_counts = ";".join(
                "%s|%s:%d" % (tube, umi, count)
                for (tube, umi), count in sorted(family_counts.items())
            ) or "."
            unique_umis = len(set(umi for _, umi in family_counts))
            technical_intervals = detect_technical_intervals(seq)
            technical_intervals_text = ";".join(
                "%d-%d" % (start0, end0) for start0, end0 in technical_intervals
            ) or "."
            technical_bases = sum(end0 - start0 for start0, end0 in technical_intervals)
            writer.writerow([
                args.sample, sequence_id, len(seq), raw_count, eligible_count,
                unique_umis, len(family_counts), umi_counts,
                tube_counts.get("A", 0), tube_counts.get("B", 0), tube_counts.get("UNKNOWN", 0),
                missing_count, collision_count, technical_intervals_text, technical_bases, seq,
            ])
            fasta_out.write(">%s count=%d umi_families=%d\n%s\n" % (
                sequence_id, raw_count, len(family_counts), wrap_fasta(seq)
            ))

    summary = dict(sorted(stats.items()))
    summary.update({
        "sample": args.sample,
        "merged_fastq": str(Path(args.merged_fastq).resolve()),
        "umi_eligible_statuses": sorted(UMI_ELIGIBLE_STATUSES),
        "umi_family_definition": "sample+tube+event+exact_umi",
        "umi_error_correction": "none_exact_only",
        "max_missing_raw_id_rate": args.max_missing_raw_id_rate,
    })
    with open(summary_tmp, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2, sort_keys=True)
        handle.write("\n")

    os.replace(fasta_tmp, fasta_path)
    os.replace(metadata_tmp, metadata_path)
    os.replace(summary_tmp, summary_path)
    if not args.keep_umi_index and index_path.exists():
        index_path.unlink()

    print("[OK] unique_sequences=%d merged_reads=%d" % (len(sequences), stats["merged_reads"]), flush=True)
    print("[OUT] %s" % fasta_path, flush=True)
    print("[OUT] %s" % metadata_path, flush=True)


if __name__ == "__main__":
    main()
