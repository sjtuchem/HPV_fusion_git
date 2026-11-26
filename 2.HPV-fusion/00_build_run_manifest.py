#!/usr/bin/env python3
"""Create a reproducible per-sample input/code/parameter signature."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--input", action="append", default=[])
    parser.add_argument("--reference", action="append", default=[])
    parser.add_argument("--script", action="append", default=[])
    parser.add_argument("--parameter", action="append", default=[])
    parser.add_argument("--output", required=True)
    return parser.parse_args()


def split_named(value: str) -> tuple[str, str]:
    if "=" not in value:
        raise ValueError("Expected NAME=VALUE, received: %s" % value)
    name, item = value.split("=", 1)
    if not name or not item:
        raise ValueError("Expected non-empty NAME=VALUE, received: %s" % value)
    return name, item


def file_info(values: List[str], hash_content: bool) -> Dict[str, dict]:
    result: Dict[str, dict] = {}
    for value in values:
        name, raw_path = split_named(value)
        path = Path(raw_path)
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError("Missing or empty file for %s: %s" % (name, path))
        stat = path.stat()
        item = {"path": str(path.resolve()), "size": stat.st_size}
        if hash_content:
            digest = hashlib.sha256()
            with path.open("rb") as handle:
                for chunk in iter(lambda: handle.read(1024 * 1024), b""):
                    digest.update(chunk)
            item["sha256"] = digest.hexdigest()
        else:
            item["mtime_ns"] = stat.st_mtime_ns
        result[name] = item
    return result


def named_values(values: List[str]) -> Dict[str, str]:
    result = {}
    for value in values:
        name, item = split_named(value)
        result[name] = item
    return result


def main() -> None:
    args = parse_args()
    signed = {
        "sample": args.sample,
        "pipeline_version": args.pipeline_version,
        "inputs": file_info(args.input, hash_content=False),
        "references": file_info(args.reference, hash_content=False),
        "scripts": file_info(args.script, hash_content=True),
        "parameters": named_values(args.parameter),
    }
    canonical = json.dumps(signed, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
    signature = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    document = dict(signed)
    document.update({
        "signature": signature,
        "generated_at": datetime.now(timezone.utc).astimezone().isoformat(),
    })
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + ".tmp.%d" % os.getpid())
    temporary.write_text(
        json.dumps(document, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    temporary.replace(output)
    print(signature)


if __name__ == "__main__":
    main()
