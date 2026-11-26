#!/usr/bin/env python3
"""Verify raw FASTQ files against one-file or directory-level MD5 manifests."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path


MD5_RE = re.compile(r"(?i)\b([0-9a-f]{32})\b")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--pair",
        action="append",
        nargs=2,
        metavar=("FASTQ", "MD5_FILE"),
        required=True,
        help="FASTQ and its MD5 sidecar; may be supplied more than once",
    )
    return parser.parse_args()


def file_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        while True:
            chunk = handle.read(8 * 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
    return digest.hexdigest()


def expected_md5(path: Path, fastq: Path) -> str:
    """Return the checksum belonging to *fastq* from a sidecar/manifest.

    A sample directory can contain a shared ``MD5.txt`` with both R1 and R2.
    Taking the first checksum from that file silently assigns the R1 checksum to
    R2, so multi-entry files must be resolved by basename.
    """
    text = path.read_text(encoding="utf-8", errors="replace")
    matches = []
    filename_matches = []
    basename = fastq.name
    for line in text.splitlines():
        match = MD5_RE.search(line)
        if match is None:
            continue
        value = match.group(1).lower()
        matches.append(value)
        tokens = [item.strip("*:") for item in re.split(r"[\s=()]+", line) if item]
        if any(Path(item).name == basename for item in tokens):
            filename_matches.append(value)
    if len(filename_matches) == 1:
        return filename_matches[0]
    if len(filename_matches) > 1:
        if len(set(filename_matches)) == 1:
            return filename_matches[0]
        raise ValueError("Conflicting MD5 values for %s in %s" % (basename, path))
    if len(matches) == 1:
        return matches[0]
    if not matches:
        raise ValueError("No 32-character MD5 value found in %s" % path)
    raise ValueError("Cannot identify the MD5 for %s in multi-entry file %s" % (basename, path))


def main() -> None:
    args = parse_args()
    report = []
    for fastq_value, md5_value in args.pair:
        fastq = Path(fastq_value)
        sidecar = Path(md5_value)
        if not fastq.is_file() or fastq.stat().st_size == 0:
            raise FileNotFoundError("Missing or empty FASTQ: %s" % fastq)
        if not sidecar.is_file() or sidecar.stat().st_size == 0:
            raise FileNotFoundError("Missing or empty MD5 sidecar: %s" % sidecar)
        expected = expected_md5(sidecar, fastq)
        observed = file_md5(fastq)
        status = "PASS" if observed == expected else "FAIL"
        report.append({
            "fastq": str(fastq),
            "md5_file": str(sidecar),
            "expected_md5": expected,
            "observed_md5": observed,
            "status": status,
        })
        if status != "PASS":
            raise ValueError(
                "MD5 mismatch for %s: expected %s, observed %s" %
                (fastq, expected, observed)
            )
    print(json.dumps({"status": "PASS", "files": report}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
