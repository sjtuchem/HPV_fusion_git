#!/usr/bin/env python3
"""Dependency-free formal-pipeline checks for MD5, UMI, coordinates and collapse."""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path

from fusion_common import cigar_metrics, classify_pair_umi, reference_boundary, revcomp


HERE = Path(__file__).resolve().parent


def load_script_module(filename: str, module_name: str):
    spec = importlib.util.spec_from_file_location(module_name, HERE / filename)
    if spec is None or spec.loader is None:
        raise RuntimeError("Cannot load %s" % filename)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_umi() -> None:
    umi = "ACT" * 5
    empty = "G" * 40
    r2_forward = "TT" + umi + "GTGCAAG" + "N" * 11
    called, status, _ = classify_pair_umi(empty, r2_forward)
    assert (called, status) == (umi, "OK_R2_ONLY")

    r1_reverse = "TTCTTGCAC" + revcomp(umi) + "N" * 11
    called, status, _ = classify_pair_umi(r1_reverse, r2_forward)
    assert (called, status) == (umi, "OK_R1_R2_CONCORDANT")

    invalid = "G" * 15 + "GTGCAAG" + "N" * 11
    called, status, _ = classify_pair_umi(empty, invalid)
    assert called == "" and status == "NOT_FOUND"


def test_cigar_coordinates() -> None:
    plus = cigar_metrics("40S60M", 100, False)
    assert (plus.query_start0, plus.query_end0, plus.match_bases) == (40, 100, 60)
    reverse = cigar_metrics("40S60M", 100, True)
    assert (reverse.query_start0, reverse.query_end0, reverse.match_bases) == (0, 60, 60)
    assert reference_boundary(101, 160, False, "start") == 101
    assert reference_boundary(101, 160, False, "end") == 160
    assert reference_boundary(101, 160, True, "start") == 160
    assert reference_boundary(101, 160, True, "end") == 101


def test_terminal_rescue_bounds() -> None:
    rescue = load_script_module("04_prepare_technical_rescue.py", "technical_rescue_test_module")
    assert rescue.trim_bounds(150, [(117, 150)], 10) == (0, 117, "TRIM_RIGHT")
    assert rescue.trim_bounds(150, [(0, 33)], 10) == (33, 150, "TRIM_LEFT")
    assert rescue.trim_bounds(200, [(0, 33), (167, 200)], 10) == (33, 167, "TRIM_BOTH")
    assert rescue.trim_bounds(200, [(80, 113)], 10)[2] == "INTERNAL_TECHNICAL_INTERVAL"


def test_directional_umi_correction() -> None:
    finalize = load_script_module("05_finalize_fusion_results.py", "finalize_v3_test_module")
    counts = {"AAAAAAAAAAAAAAA": 10, "AAAAAAAAAAAAAAT": 1, "CCCCCCCCCCCCCCC": 1}
    assigned = finalize.directional_umi_assignments(counts, 1)
    assert assigned["AAAAAAAAAAAAAAA"] == "AAAAAAAAAAAAAAA"
    assert assigned["AAAAAAAAAAAAAAT"] == "AAAAAAAAAAAAAAA"
    assert assigned["CCCCCCCCCCCCCCC"] == "CCCCCCCCCCCCCCC"


def write_metadata(
    path: Path,
    sequence_id: str,
    sequence: str,
    technical_intervals0: str = ".",
    technical_bases: int = 0,
) -> None:
    columns = [
        "sample", "sequence_id", "sequence_length", "raw_read_count",
        "umi_eligible_read_count", "unique_umi_count", "unique_umi_family_count",
        "umi_family_counts", "tube_A_read_count", "tube_B_read_count",
        "unknown_tube_read_count", "missing_raw_id_count", "ambiguous_read_id_count",
        "technical_intervals0", "technical_bases", "sequence",
    ]
    row = {
        "sample": "Tissue1", "sequence_id": sequence_id, "sequence_length": len(sequence),
        "raw_read_count": 5, "umi_eligible_read_count": 4, "unique_umi_count": 2,
        "unique_umi_family_count": 2, "umi_family_counts": "A|AAAAAAAAAAAAAAA:3;B|CCCCCCCCCCCCCCC:1",
        "tube_A_read_count": 4, "tube_B_read_count": 1, "unknown_tube_read_count": 0,
        "missing_raw_id_count": 0, "ambiguous_read_id_count": 0,
        "technical_intervals0": technical_intervals0,
        "technical_bases": technical_bases,
        "sequence": sequence,
    }
    with gzip.open(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)


def write_fastq(path: Path, records) -> None:
    with gzip.open(path, "wt", encoding="ascii") as handle:
        for read_id, sequence in records:
            handle.write("@%s\n%s\n+\n%s\n" % (read_id, sequence, "I" * len(sequence)))


def test_md5_verifier() -> None:
    with tempfile.TemporaryDirectory(prefix="hpv_fusion_md5_selftest_") as tmp_value:
        tmp = Path(tmp_value)
        fastq = tmp / "tiny.fastq.gz"
        write_fastq(fastq, [("R1", "ACGT")])
        expected = hashlib.md5(fastq.read_bytes()).hexdigest()
        sidecar = tmp / "tiny.fastq.gz.md5"
        sidecar.write_text("MD5 (tiny.fastq.gz) = %s\n" % expected, encoding="ascii")
        result = subprocess.run([
            sys.executable, str(HERE / "00_verify_md5.py"),
            "--pair", str(fastq), str(sidecar),
        ], check=True, stdout=subprocess.PIPE, text=True)
        assert '"status": "PASS"' in result.stdout


def test_prepare_and_collapse() -> None:
    umi_a = "ACT" * 5
    umi_b = "CAT" * 5
    r1_sequence = "G" * 40
    merged_sequence = "AC" * 50
    with tempfile.TemporaryDirectory(prefix="hpv_fusion_prepare_selftest_") as tmp_value:
        tmp = Path(tmp_value)
        r1_a, r2_a = tmp / "A_R1.fastq.gz", tmp / "A_R2.fastq.gz"
        r1_b, r2_b = tmp / "B_R1.fastq.gz", tmp / "B_R2.fastq.gz"
        merged = tmp / "merged.fastq.gz"
        write_fastq(r1_a, [("A1/1", r1_sequence), ("A2/1", r1_sequence)])
        write_fastq(r2_a, [
            ("A1/2", umi_a + "GTGCAAG" + "N" * 11),
            ("A2/2", umi_b + "GTGCAAG" + "N" * 11),
        ])
        write_fastq(r1_b, [("B1/1", r1_sequence)])
        write_fastq(r2_b, [("B1/2", umi_a + "GTGCAAG" + "N" * 11)])
        write_fastq(merged, [("A1", merged_sequence), ("A2", merged_sequence), ("B1", merged_sequence)])
        out_dir = tmp / "work"
        subprocess.run([
            sys.executable, str(HERE / "01_prepare_unique_sequences.py"),
            "--sample", "Tissue1", "--merged-fastq", str(merged),
            "--r1-a", str(r1_a), "--r2-a", str(r2_a),
            "--r1-b", str(r1_b), "--r2-b", str(r2_b), "--out-dir", str(out_dir),
        ], check=True, stdout=subprocess.PIPE, text=True)
        with gzip.open(out_dir / "Tissue1.unique_sequences.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        row = rows[0]
        assert row["raw_read_count"] == "3"
        assert row["umi_eligible_read_count"] == "3"
        assert row["unique_umi_count"] == "2"
        assert row["unique_umi_family_count"] == "3"
        assert row["tube_A_read_count"] == "2" and row["tube_B_read_count"] == "1"


def test_prepare_rescue_io() -> None:
    umi = "ACT" * 5
    biological = "A" * 67
    sequence = biological + umi + "GTGCAAG" + "N" * 11
    sequence2 = biological + "C" * 33
    sequence_id = "SEQ_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()
    sequence_id2 = "SEQ_" + hashlib.sha256(sequence2.encode("ascii")).hexdigest()
    with tempfile.TemporaryDirectory(prefix="hpv_fusion_rescue_selftest_") as tmp_value:
        tmp = Path(tmp_value)
        metadata = tmp / "metadata.tsv.gz"
        write_metadata(metadata, sequence_id, sequence, "67-100", 33)
        with gzip.open(metadata, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = reader.fieldnames
            metadata_rows = list(reader)
        second = dict(metadata_rows[0])
        second.update({"sequence_id": sequence_id2, "sequence": sequence2})
        with gzip.open(metadata, "wt", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows([metadata_rows[0], second])
        evidence = tmp / "evidence.tsv.gz"
        evidence_columns = [
            "sequence_id", "event_id", "pass_length_30bp", "pass_highconf_30bp",
            "full_human_like", "human_mapq", "human_nh", "hpv_mapq",
            "hpv_technical_overlap_bases", "human_technical_overlap_bases",
        ]
        with gzip.open(evidence, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=evidence_columns, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerow({
                "sequence_id": sequence_id, "event_id": "EVT_TEST",
                "pass_length_30bp": 1, "pass_highconf_30bp": 0,
                "full_human_like": 0, "human_mapq": 60, "human_nh": 1,
                "hpv_mapq": 60, "hpv_technical_overlap_bases": 0,
                "human_technical_overlap_bases": 10,
            })
            second_evidence = {
                "sequence_id": sequence_id2, "event_id": "EVT_TEST",
                "pass_length_30bp": 1, "pass_highconf_30bp": 0,
                "full_human_like": 0, "human_mapq": 60, "human_nh": 1,
                "hpv_mapq": 60, "hpv_technical_overlap_bases": 0,
                "human_technical_overlap_bases": 10,
            }
            writer.writerow(second_evidence)
        subprocess.run([
            sys.executable, str(HERE / "04_prepare_technical_rescue.py"),
            "--strict-evidence", str(evidence), "--sequence-metadata", str(metadata),
            "--out-fasta", str(tmp / "rescue.fa"),
            "--out-metadata", str(tmp / "rescue.tsv.gz"),
            "--out-map", str(tmp / "rescue_map.tsv"),
            "--summary", str(tmp / "rescue.json"),
        ], check=True, stdout=subprocess.PIPE, text=True)
        with open(tmp / "rescue_map.tsv", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 2
        assert {row["original_sequence_id"] for row in rows} == {sequence_id, sequence_id2}
        assert len({row["rescue_sequence_id"] for row in rows}) == 1
        assert len({row["normalized_sequence_id"] for row in rows}) == 1
        assert all(row["trim_start0"] == "0" and row["trim_end0"] == "67" for row in rows)
        rescue_text = (tmp / "rescue.fa").read_text(encoding="ascii")
        assert rescue_text.count(">") == 1
        assert rescue_text.rstrip().endswith(biological)
        with gzip.open(tmp / "rescue.tsv.gz", "rt", encoding="utf-8", newline="") as handle:
            rescue_metadata = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rescue_metadata) == 1
        assert rescue_metadata[0]["raw_read_count"] == "10"
        with open(tmp / "rescue.json", encoding="utf-8") as handle:
            import json
            rescue_summary = json.load(handle)
        assert rescue_summary["rescue_original_sequences_written"] == 2
        assert rescue_summary["rescue_unique_trimmed_sequences_written"] == 1


def test_synthetic_breakpoint_call() -> None:
    sequence = "A" * 40 + "C" * 60
    sequence_id = "SEQ_" + hashlib.sha256(sequence.encode("ascii")).hexdigest()
    with tempfile.TemporaryDirectory(prefix="hpv_fusion_selftest_") as tmp_value:
        tmp = Path(tmp_value)
        unique_fasta = tmp / "unique.fa"
        unique_fasta.write_text(">%s\n%s\n" % (sequence_id, sequence), encoding="ascii")
        bwa_sam = tmp / "hpv.sam"
        bwa_sam.write_text(
            "@HD\tVN:1.6\n"
            "@SQ\tSN:K02718.1\tLN:7904\n"
            "%s\t0\tK02718.1\t101\t60\t40S60M\t*\t0\t0\t%s\t%s\tNM:i:0\tAS:i:60\n"
            % (sequence_id, sequence, "I" * len(sequence)),
            encoding="utf-8",
        )
        hpv_tsv = tmp / "hpv.tsv"
        star_fasta = tmp / "star.fa"
        subprocess.run([
            sys.executable, str(HERE / "02_select_hpv_candidates.py"),
            "--sam", str(bwa_sam), "--unique-fasta", str(unique_fasta),
            "--out-alignments", str(hpv_tsv), "--out-star-fasta", str(star_fasta),
            "--summary", str(tmp / "hpv.json"), "--prefilter-min-hpv", "15",
            "--prefilter-min-flank", "15",
        ], check=True, stdout=subprocess.PIPE, text=True)
        star_text = star_fasta.read_text(encoding="ascii")
        assert "|FULL|0|100" in star_text and "|L|0|40" in star_text

        metadata = tmp / "metadata.tsv.gz"
        write_metadata(metadata, sequence_id, sequence)
        human_sam = tmp / "human.sam"
        human_sam.write_text(
            "@HD\tVN:1.6\n"
            "@SQ\tSN:chr1\tLN:1000000\n"
            "%s|L|0|40\t0\tchr1\t1001\t60\t40M\t*\t0\t0\t%s\t%s\tNH:i:1\tNM:i:0\tAS:i:40\n"
            % (sequence_id, "A" * 40, "I" * 40),
            encoding="utf-8",
        )
        hpv_map = tmp / "hpv.csv"
        hpv_map.write_text("GenBank ID,Virus name\nK02718,HPV16\n", encoding="utf-8")
        gene_bed = tmp / "genes.bed"
        gene_bed.write_text("chr1\t900\t1100\tGENE1\n", encoding="utf-8")
        final = tmp / "final"
        subprocess.run([
            sys.executable, str(HERE / "03_call_and_aggregate_fusions.py"),
            "--sample", "Tissue1", "--human-sam", str(human_sam),
            "--hpv-alignments", str(hpv_tsv), "--sequence-metadata", str(metadata),
            "--hpv-type-map", str(hpv_map), "--gene-bed", str(gene_bed),
            "--out-dir", str(final), "--thresholds", "20,30,40,50",
        ], check=True, stdout=subprocess.PIPE, text=True)

        with open(final / "Tissue1.fusion_events.tsv", "r", encoding="utf-8", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        assert len(rows) == 1
        event = rows[0]
        assert event["hpv_type"] == "HPV16"
        assert event["hpv_breakpoint_1based"] == "101"
        assert event["human_breakpoint_1based"] == "1040"
        assert event["genes"] == "GENE1"
        assert event["raw_reads_30bp"] == "5"
        assert event["unique_umi_families_30bp"] == "2"
        assert event["unique_sequences_50bp"] == "0"
        assert event["recommended_pass_30bp"] == "1"

        audit = tmp / "audit_only"
        subprocess.run([
            sys.executable, str(HERE / "03_call_and_aggregate_fusions.py"),
            "--sample", "Tissue1", "--human-sam", str(human_sam),
            "--hpv-alignments", str(hpv_tsv), "--sequence-metadata", str(metadata),
            "--hpv-type-map", str(hpv_map), "--gene-bed", str(gene_bed),
            "--out-dir", str(audit), "--thresholds", "20,30,40,50", "--audit-only",
        ], check=True, stdout=subprocess.PIPE, text=True)
        assert (audit / "Tissue1.fusion_evidence.tsv.gz").is_file()
        assert (audit / "Tissue1.fusion_summary.json").is_file()
        assert not (audit / "Tissue1.fusion_events.tsv").exists()

        compact = tmp / "compact"
        subprocess.run([
            sys.executable, str(HERE / "05_finalize_fusion_results.py"),
            "--sample", "Tissue1",
            "--strict-evidence", str(final / "Tissue1.fusion_evidence.tsv.gz"),
            "--sequence-metadata", str(metadata),
            "--out-dir", str(compact), "--thresholds", "20,30,40,50",
            "--primary-threshold", "30",
        ], check=True, stdout=subprocess.PIPE, text=True)
        with open(compact / "Tissue1.fusion_events.main.tsv", newline="") as handle:
            compact_events = list(csv.DictReader(handle, delimiter="\t"))
        assert len(compact_events) == 1
        assert compact_events[0]["result_tier"] == "MAIN"
        assert compact_events[0]["raw_reads_30bp"] == "5"
        assert compact_events[0]["unique_umi_families_30bp"] == "2"
        assert compact_events[0]["corrected_umi_families_30bp"] == "2"
        assert compact_events[0]["support_band"] == "40BP_ONLY"
        assert compact_events[0]["report_class"] == "MAIN_40BP_ONLY"
        with gzip.open(
            compact / "Tissue1.fusion_unique_sequences.final.tsv.gz", "rt", newline=""
        ) as handle:
            compact_sequences = list(csv.DictReader(handle, delimiter="\t"))
        assert len(compact_sequences) == 1
        assert "umi_family_counts" not in compact_sequences[0]
        assert compact_sequences[0]["original_sequence_count"] == "1"
        with open(compact / "Tissue1.fusion_loci.final.tsv", newline="") as handle:
            loci = list(csv.DictReader(handle, delimiter="\t"))
        assert len(loci) == 1 and loci[0]["exact_event_count"] == "1"
        with gzip.open(
            compact / "Tissue1.fusion_sequence_members.final.tsv.gz", "rt", newline=""
        ) as handle:
            members = list(csv.DictReader(handle, delimiter="\t"))
        assert len(members) == 1 and members[0]["normalization_type"] == "FULL_LENGTH"
        validation = subprocess.run([
            sys.executable, str(HERE / "98_validate_results.py"),
            "--sample", "Tissue1", "--final-dir", str(compact),
            "--primary-threshold", "30",
        ], check=True, stdout=subprocess.PIPE, text=True)
        assert '"status": "PASS"' in validation.stdout


def main() -> None:
    test_md5_verifier()
    test_umi()
    test_cigar_coordinates()
    test_terminal_rescue_bounds()
    test_directional_umi_correction()
    test_prepare_and_collapse()
    test_prepare_rescue_io()
    test_synthetic_breakpoint_call()
    print("[OK] all synthetic tests passed")


if __name__ == "__main__":
    main()
