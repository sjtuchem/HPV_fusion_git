#!/usr/bin/env python3
"""Dependency-free regression tests for batch plumbing and zero-event output."""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
import subprocess
import sys
import tempfile
import os
from pathlib import Path
from typing import Optional


HERE = Path(__file__).resolve().parent


def run(*arguments: str, stdout_path: Optional[Path] = None) -> subprocess.CompletedProcess:
    result = subprocess.run(
        [sys.executable, *arguments], text=True, capture_output=True, check=False
    )
    if result.returncode:
        raise AssertionError(
            "command failed (%d): %s\nstdout:\n%s\nstderr:\n%s" %
            (result.returncode, " ".join(arguments), result.stdout, result.stderr)
        )
    if stdout_path is not None:
        stdout_path.write_text(result.stdout, encoding="utf-8")
    return result


def write_fastq(path: Path, read_id: str, sequence: str) -> None:
    with gzip.open(path, "wt", encoding="ascii") as handle:
        handle.write("@%s\n%s\n+\n%s\n" % (read_id, sequence, "I" * len(sequence)))


def test_single_pair_prepare(root: Path) -> tuple[Path, Path, Path]:
    r1 = root / "raw_R1.fastq.gz"
    r2 = root / "raw_R2.fastq.gz"
    merged = root / "merged.fastq.gz"
    write_fastq(r1, "read1", "A" * 100)
    write_fastq(r2, "read1", "C" * 100)
    write_fastq(merged, "read1", "A" * 80)
    output = root / "prepare"
    run(
        str(HERE / "01_prepare_unique_sequences.py"),
        "--sample", "S0", "--merged-fastq", str(merged),
        "--r1-a", str(r1), "--r2-a", str(r2), "--out-dir", str(output),
    )
    summary = json.loads((output / "S0.prepare_summary.json").read_text(encoding="utf-8"))
    assert summary["merged_reads"] == 1
    assert summary["raw_pairs_A"] == 1
    assert summary.get("raw_pairs_B", 0) == 0
    return merged, r1, r2


def test_shared_md5(root: Path, r1: Path, r2: Path) -> None:
    manifest = root / "MD5.txt"
    lines = []
    for path in (r1, r2):
        lines.append("%s  %s" % (hashlib.md5(path.read_bytes()).hexdigest(), path.name))
    manifest.write_text("\n".join(lines) + "\n", encoding="utf-8")
    report = run(
        str(HERE / "00_verify_md5.py"),
        "--pair", str(r1), str(manifest), "--pair", str(r2), str(manifest),
    )
    assert json.loads(report.stdout)["status"] == "PASS"


def test_targeted_umi_guard(root: Path) -> None:
    summary = root / "umi_summary.json"
    valid = {
        "sample": "TARGETED",
        "raw_pairs_A": 100,
        "raw_umi_eligible_A": 70,
        "raw_pairs_B": 120,
        "raw_umi_eligible_B": 80,
        "merged_reads": 200,
        "merged_umi_eligible_reads": 140,
    }
    summary.write_text(json.dumps(valid), encoding="utf-8")
    report = run(
        str(HERE / "00_validate_umi_summary.py"), "--summary", str(summary),
        "--required-tube", "A", "--required-tube", "B",
    )
    assert json.loads(report.stdout)["status"] == "PASS"

    invalid = dict(valid)
    invalid["raw_umi_eligible_B"] = 0
    summary.write_text(json.dumps(invalid), encoding="utf-8")
    rejected = subprocess.run(
        [sys.executable, str(HERE / "00_validate_umi_summary.py"),
         "--summary", str(summary), "--required-tube", "A", "--required-tube", "B"],
        text=True, capture_output=True, check=False,
    )
    assert rejected.returncode != 0
    assert "raw tube B contains zero UMI-eligible read pairs" in rejected.stdout


def test_sample_list_and_manifest(root: Path, merged: Path, r1: Path, r2: Path) -> Path:
    sample_list = root / "sample.list"
    sample_list.write_text(
        "sample\tmerged_fastq\traw_r1_a\traw_r2_a\traw_r1_b\traw_r2_b\n"
        "S0\t%s\t%s\t%s\t.\t.\n" % (merged, r1, r2), encoding="utf-8"
    )
    jobs = root / "jobs.tsv"
    run(
        str(HERE / "00_validate_sample_list.py"), "--sample-list", str(sample_list),
        "--output", str(jobs), "--report", str(root / "list_report.json"),
    )
    assert jobs.read_text(encoding="utf-8").count("\n") == 1
    duplicate = root / "duplicate.list"
    duplicate.write_text(
        "sample\tmerged_fastq\traw_r1_a\traw_r2_a\traw_r1_b\traw_r2_b\n"
        + ("S0\t%s\t%s\t%s\t.\t.\n" % (merged, r1, r2)) * 2,
        encoding="utf-8",
    )
    rejected = subprocess.run(
        [sys.executable, str(HERE / "00_validate_sample_list.py"),
         "--sample-list", str(duplicate), "--output", str(root / "bad_jobs.tsv"),
         "--report", str(root / "bad_report.json")],
        text=True, capture_output=True, check=False,
    )
    assert rejected.returncode != 0 and "duplicate sample name" in rejected.stderr
    half_pair = root / "half_pair.list"
    half_pair.write_text(
        "sample\tmerged_fastq\traw_r1_a\traw_r2_a\traw_r1_b\traw_r2_b\n"
        "S1\t%s\t%s\t%s\t%s\t.\n" % (merged, r1, r2, r1), encoding="utf-8",
    )
    rejected = subprocess.run(
        [sys.executable, str(HERE / "00_validate_sample_list.py"),
         "--sample-list", str(half_pair), "--output", str(root / "bad_jobs2.tsv"),
         "--report", str(root / "bad_report2.json")],
        text=True, capture_output=True, check=False,
    )
    assert rejected.returncode != 0 and "must both be paths or both be ." in rejected.stderr
    manifest = root / "run_manifest.json"
    first = run(
        str(HERE / "00_build_run_manifest.py"), "--sample", "S0",
        "--pipeline-version", "test", "--input", "merged=%s" % merged,
        "--reference", "raw=%s" % r1, "--script", "self=%s" % __file__,
        "--parameter", "threshold=30", "--output", str(manifest),
    ).stdout.strip()
    second = run(
        str(HERE / "00_build_run_manifest.py"), "--sample", "S0",
        "--pipeline-version", "test", "--input", "merged=%s" % merged,
        "--reference", "raw=%s" % r1, "--script", "self=%s" % __file__,
        "--parameter", "threshold=30", "--output", str(manifest),
    ).stdout.strip()
    assert first == second == json.loads(manifest.read_text(encoding="utf-8"))["signature"]
    return jobs


def test_zero_event_result(root: Path, jobs: Path) -> None:
    sample_root = root / "result" / "S0"
    final = sample_root / "final"
    logs = sample_root / "logs"
    final.mkdir(parents=True)
    logs.mkdir(parents=True)
    strict = root / "strict.tsv.gz"
    metadata = root / "metadata.tsv.gz"
    with gzip.open(strict, "wt", encoding="utf-8") as handle:
        handle.write("sample\n")
    with gzip.open(metadata, "wt", encoding="utf-8") as handle:
        handle.write("sample\n")
    run(
        str(HERE / "05_finalize_fusion_results.py"),
        "--sample", "S0", "--strict-evidence", str(strict),
        "--sequence-metadata", str(metadata), "--out-dir", str(final),
        "--pipeline-version", "test",
    )
    validation = logs / "11_validate_results.log"
    run(
        str(HERE / "98_validate_results.py"), "--sample", "S0",
        "--final-dir", str(final), stdout_path=validation,
    )
    assert json.loads(validation.read_text(encoding="utf-8"))["status"] == "PASS"
    (logs / "S0.complete").write_text("status=PASS\nsignature=test\n", encoding="utf-8")
    (logs / "S0.runtime_metrics.tsv").write_text(
        "stage\tstatus\tduration_seconds\nTOTAL\tRUN\t0\n", encoding="utf-8"
    )
    run(
        str(HERE / "00_batch_summary.py"), "--jobs", str(jobs),
        "--result-dir", str(root / "result"),
        "--output-tsv", str(root / "batch_summary.tsv"),
        "--output-json", str(root / "batch_summary.json"),
    )
    report = json.loads((root / "batch_summary.json").read_text(encoding="utf-8"))
    assert report["status_counts"]["COMPLETE_PASS"] == 1
    summary = json.loads((final / "S0.fusion_summary.final.json").read_text(encoding="utf-8"))
    assert summary["final_exact_events"] == 0 and summary["final_loci"] == 0


def test_zero_candidate_runner_and_skip(root: Path, merged: Path, r1: Path, r2: Path) -> None:
    fake_bin = root / "fake_bin"
    fake_bin.mkdir()
    bwa = fake_bin / "bwa"
    bwa.write_text("#!/usr/bin/env python3\nprint('@HD\\tVN:1.6')\n", encoding="utf-8")
    bwa.chmod(0o755)
    for name in ("STAR", "samtools"):
        path = fake_bin / name
        path.write_text("#!/usr/bin/env python3\nraise SystemExit(99)\n", encoding="utf-8")
        path.chmod(0o755)

    reference = root / "hpv.fa"
    reference.write_text(">HPV_TEST\n" + "A" * 100 + "\n", encoding="ascii")
    for suffix in ("amb", "ann", "bwt", "pac", "sa"):
        (root / ("hpv.fa." + suffix)).write_text("index\n", encoding="ascii")
    star_index = root / "star_index"
    star_index.mkdir()
    for name in ("Genome", "SA", "SAindex"):
        (star_index / name).write_text("index\n", encoding="ascii")
    gene_bed = root / "genes.bed"
    gene_bed.write_text("chr1\t0\t100\tTESTGENE\n", encoding="utf-8")
    run_root = root / "runner"
    config = root / "runner_config.sh"
    config.write_text(
        "source '%s'\n" % (HERE / "config.sh")
        + "CODE_DIR='%s'\n" % HERE
        + "RESULT_DIR='%s'\n" % (run_root / "result")
        + "HPV_REF='%s'\n" % reference
        + "STAR_INDEX='%s'\n" % star_index
        + "GENE_BED='%s'\n" % gene_bed
        + "HPV_TYPE_MAP='%s'\n" % (HERE / "hpv_GenBankID_fromIHPRC.csv")
        + "THREADS_PER_SAMPLE=1\nVERIFY_RAW_MD5=0\nREQUIRE_UMI=0\nFORCE=0\n",
        encoding="utf-8",
    )
    environment = os.environ.copy()
    environment["PATH"] = str(fake_bin) + os.pathsep + environment["PATH"]
    environment["CONFIG_FILE"] = str(config)
    arguments = [
        "bash", str(HERE / "01_run_one_sample.sh"), "NEG", str(merged),
        str(r1), str(r2), ".", ".",
    ]
    first = subprocess.run(arguments, text=True, capture_output=True, env=environment, check=False)
    if first.returncode:
        raise AssertionError("zero-candidate runner failed:\n%s\n%s" % (first.stdout, first.stderr))
    assert "No HPV partial-alignment candidates" in first.stdout
    assert (run_root / "result/NEG/logs/NEG.complete").is_file()
    second = subprocess.run(arguments, text=True, capture_output=True, env=environment, check=False)
    if second.returncode:
        raise AssertionError("completed-sample skip failed:\n%s\n%s" % (second.stdout, second.stderr))
    assert "[SKIP COMPLETE]" in second.stdout

    unicode_data = root / "中文数据"
    unicode_data.mkdir()
    unicode_inputs = []
    for source in (merged, r1, r2):
        destination = unicode_data / source.name
        destination.symlink_to(source)
        unicode_inputs.append(destination)
    unicode_batch = root / "批次日志"
    unicode_batch.mkdir()
    jobs = unicode_batch / "validated_jobs.tsv"
    jobs.write_text(
        "NEG\t%s\t%s\t%s\t.\t.\n" % tuple(unicode_inputs), encoding="utf-8"
    )
    environment["VALIDATED_JOBS_FILE"] = str(jobs)
    dispatched = subprocess.run(
        ["bash", str(HERE / "00_run_one_logged.sh"), "NEG"],
        text=True, capture_output=True, env=environment, check=False,
    )
    if dispatched.returncode:
        raise AssertionError("ASCII-only sample dispatch failed:\n%s\n%s" % (
            dispatched.stdout, dispatched.stderr
        ))
    assert "[SKIP COMPLETE]" in dispatched.stdout
    assert (run_root / "result/NEG/logs/NEG.pipeline.log").is_file()

    stat = merged.stat()
    os.utime(merged, ns=(stat.st_atime_ns, stat.st_mtime_ns + 1_000_000_000))
    stale = subprocess.run(arguments, text=True, capture_output=True, env=environment, check=False)
    assert stale.returncode != 0
    assert "different inputs, code, references, or thresholds" in stale.stdout + stale.stderr


def main() -> None:
    with tempfile.TemporaryDirectory(prefix="hpv_batch_selftest_") as value:
        root = Path(value)
        merged, r1, r2 = test_single_pair_prepare(root)
        test_shared_md5(root, r1, r2)
        test_targeted_umi_guard(root)
        jobs = test_sample_list_and_manifest(root, merged, r1, r2)
        test_zero_event_result(root, jobs)
        test_zero_candidate_runner_and_skip(root, merged, r1, r2)
    print("[PASS] batch plumbing, ASCII dispatch, targeted UMI guard, shared MD5, zero-event run, skip and stale-input guard")


if __name__ == "__main__":
    main()
