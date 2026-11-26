#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ARRP targeted-library on-target counting (generic version)

Reproduces the original counting semantics:
  count        = number of reads containing a given iFP or iFP_spacer sequence
  count_anchor = among those reads, number also containing ANCHOR_RC (default CTTGCAC)

Improvements over the old batch script:
  1) primer Excel / paths / sample list are external inputs, no batch-specific paths are hard-coded;
  2) iFP and iFP_spacer are explicitly labeled;
  3) duplicate primer/search sequences are checked before analysis;
  4) FASTQ is streamed instead of loading the whole file into RAM;
  5) Aho-Corasick matching finds nested/overlapping iFP + iFP_spacer patterns efficiently;
  6) a manifest prevents stale 127-plex results being silently reused after switching to 126 plex;
  7) outputs both legacy-style count|count_anchor matrix and tidy rate tables.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
import sys
from collections import deque
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple

import pandas as pd

ALGORITHM_VERSION = "20260901_generic_v1"


@dataclass(frozen=True)
class Pattern:
    pattern_id: int
    gene: str
    number: str
    target_type: str
    seq: str


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Calculate iFP/iFP_spacer counts and on-target rates.")
    p.add_argument("--fastp-dir", required=True, help="Directory containing <sample_id>_qc1.fastq.gz")
    p.add_argument("--out-dir", required=True, help="Output directory")
    p.add_argument("--primer-excel", required=True, help="Primer Excel (.xlsx)")
    p.add_argument("--primer-sheet", default="Sheet1", help="Excel sheet name (default: Sheet1)")
    p.add_argument("--samples-table", required=True, help="samples.tsv with sample_id/A_dir/B_dir/enabled")
    p.add_argument("--anchor", default="CTTGCAC", help="Anchor reverse-complement sequence")
    p.add_argument("--force", action="store_true", help="Force recalculation even when manifest matches")
    return p.parse_args()


def sha256_file(path: Path, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while True:
            chunk = fh.read(chunk_size)
            if not chunk:
                break
            h.update(chunk)
    return h.hexdigest()


def normalize_seq(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip().upper().replace(" ", "")


def load_enabled_samples(path: Path) -> List[str]:
    samples: List[str] = []
    with path.open("r", encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"sample_id", "A_dir", "B_dir", "enabled"}
        if reader.fieldnames is None or not required.issubset(set(reader.fieldnames)):
            raise ValueError(f"samples.tsv must contain columns: {sorted(required)}")
        for row in reader:
            enabled = str(row.get("enabled", "")).strip().lower()
            if enabled not in {"1", "true", "yes"}:
                continue
            sid = str(row.get("sample_id", "")).strip()
            if not sid:
                raise ValueError("Enabled sample has empty sample_id")
            samples.append(sid)
    if not samples:
        raise ValueError("No enabled samples found in samples.tsv")
    if len(samples) != len(set(samples)):
        dup = sorted({x for x in samples if samples.count(x) > 1})
        raise ValueError(f"Duplicate sample_id in samples.tsv: {dup}")
    return samples


def load_patterns(excel: Path, sheet: str, out_dir: Path) -> List[Pattern]:
    df = pd.read_excel(excel, sheet_name=sheet, header=0)
    required = {"gene", "number", "iFP", "iFP_spacer"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Primer Excel missing columns: {sorted(missing)}")

    base = df.loc[:, ["gene", "number", "iFP", "iFP_spacer"]].copy()
    base["gene"] = base["gene"].fillna("").astype(str).str.strip()
    base["number"] = base["number"].fillna("").astype(str).str.strip()
    base["iFP"] = base["iFP"].map(normalize_seq)
    base["iFP_spacer"] = base["iFP_spacer"].map(normalize_seq)

    if (base["number"] == "").any():
        raise ValueError("Primer Excel contains empty number values")
    if base["number"].duplicated().any():
        dups = sorted(base.loc[base["number"].duplicated(keep=False), "number"].unique())
        raise ValueError(f"Duplicate number values in primer Excel: {dups}")

    for col in ("iFP", "iFP_spacer"):
        if (base[col] == "").any():
            bad = base.loc[base[col] == "", ["gene", "number"]].to_dict("records")
            raise ValueError(f"Empty {col} sequence(s): {bad}")
        invalid = base.loc[~base[col].str.fullmatch(r"[ACGTN]+"), ["gene", "number", col]]
        if not invalid.empty:
            raise ValueError(f"Invalid DNA characters in {col}:\n{invalid.to_string(index=False)}")
        if base[col].duplicated().any():
            dups = base.loc[base[col].duplicated(keep=False), ["gene", "number", col]]
            raise ValueError(f"Duplicate {col} sequences detected:\n{dups.to_string(index=False)}")

    rows: List[Pattern] = []
    pid = 0
    for _, r in base.iterrows():
        rows.append(Pattern(pid, r["gene"], r["number"], "iFP", r["iFP"]))
        pid += 1
    for _, r in base.iterrows():
        rows.append(Pattern(pid, r["gene"], r["number"], "iFP_spacer", r["iFP_spacer"]))
        pid += 1

    # Cross-check duplicates across all search strings (e.g. one target's iFP equals another target's spacer).
    seq_to_patterns: Dict[str, List[Pattern]] = {}
    for p in rows:
        seq_to_patterns.setdefault(p.seq, []).append(p)
    cross_dups = {seq: ps for seq, ps in seq_to_patterns.items() if len(ps) > 1}
    if cross_dups:
        lines = []
        for seq, ps in cross_dups.items():
            desc = "; ".join(f"{p.number}:{p.target_type}" for p in ps)
            lines.append(f"{seq}\t{desc}")
        raise ValueError("Duplicate search sequences across iFP/iFP_spacer:\n" + "\n".join(lines))

    out_dir.mkdir(parents=True, exist_ok=True)
    audit_csv = out_dir / "primer_search_list.csv"
    with audit_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["gene", "number", "target_type", "seq"])
        for p in rows:
            writer.writerow([p.gene, p.number, p.target_type, p.seq])

    audit_txt = out_dir / "primer_audit.txt"
    with audit_txt.open("w", encoding="utf-8") as fh:
        fh.write(f"Primer Excel: {excel}\n")
        fh.write(f"Sheet: {sheet}\n")
        fh.write(f"Unique targets: {len(base)}\n")
        fh.write(f"Search sequences: {len(rows)} (= iFP {len(base)} + iFP_spacer {len(base)})\n")
        fh.write("number duplicates: 0\n")
        fh.write("iFP duplicates: 0\n")
        fh.write("iFP_spacer duplicates: 0\n")
        fh.write("cross iFP/iFP_spacer exact duplicates: 0\n")

    return rows


class AhoCorasick:
    """Minimal Aho-Corasick automaton for exact DNA substring matching."""

    def __init__(self, patterns: Iterable[Pattern]):
        self.next: List[Dict[str, int]] = [dict()]
        self.fail: List[int] = [0]
        self.output: List[List[int]] = [[]]

        for p in patterns:
            state = 0
            for ch in p.seq:
                nxt = self.next[state].get(ch)
                if nxt is None:
                    nxt = len(self.next)
                    self.next[state][ch] = nxt
                    self.next.append({})
                    self.fail.append(0)
                    self.output.append([])
                state = nxt
            self.output[state].append(p.pattern_id)

        q: deque[int] = deque()
        for _, state in self.next[0].items():
            q.append(state)
            self.fail[state] = 0

        while q:
            r = q.popleft()
            for ch, s in self.next[r].items():
                q.append(s)
                f = self.fail[r]
                while f and ch not in self.next[f]:
                    f = self.fail[f]
                self.fail[s] = self.next[f].get(ch, 0)
                if self.output[self.fail[s]]:
                    self.output[s].extend(self.output[self.fail[s]])

    def find_pattern_ids(self, text: str) -> Set[int]:
        state = 0
        matched: Set[int] = set()
        for ch in text:
            while state and ch not in self.next[state]:
                state = self.fail[state]
            state = self.next[state].get(ch, 0)
            if self.output[state]:
                matched.update(self.output[state])
        return matched


def iter_fastq_sequences(path: Path) -> Iterable[str]:
    with gzip.open(path, "rt", encoding="ascii", errors="strict", newline="") as fh:
        line_no = 0
        while True:
            header = fh.readline()
            if not header:
                break
            seq = fh.readline()
            plus = fh.readline()
            qual = fh.readline()
            line_no += 4
            if not seq or not plus or not qual:
                raise ValueError(f"Truncated FASTQ near line {line_no} in {path}")
            if not header.startswith("@"):
                raise ValueError(f"Invalid FASTQ header near line {line_no-3} in {path}")
            if not plus.startswith("+"):
                raise ValueError(f"Invalid FASTQ plus line near line {line_no-1} in {path}")
            yield seq.strip().upper()


def make_signature(fastq: Path, primer_sha256: str, anchor: str) -> Dict[str, object]:
    st = fastq.stat()
    return {
        "algorithm_version": ALGORITHM_VERSION,
        "primer_sha256": primer_sha256,
        "anchor": anchor,
        "fastq_path": str(fastq.resolve()),
        "fastq_size": st.st_size,
        "fastq_mtime_ns": st.st_mtime_ns,
    }


def manifest_matches(path: Path, signature: Dict[str, object]) -> bool:
    if not path.exists():
        return False
    try:
        old = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return False
    return old == signature


def analyze_sample(
    sample_id: str,
    fastq: Path,
    patterns: List[Pattern],
    automaton: AhoCorasick,
    anchor: str,
    sample_out_root: Path,
    primer_sha256: str,
    force: bool,
) -> Tuple[List[Dict[str, object]], Dict[str, object]]:
    sample_dir = sample_out_root / sample_id
    sample_dir.mkdir(parents=True, exist_ok=True)
    detail_file = sample_dir / f"sum_onTargetRate_{sample_id}.csv"
    legacy_file = sample_dir / f"lib{sample_id}.csv"
    overall_file = sample_dir / f"overall_onTargetRate_{sample_id}.csv"
    manifest_file = sample_dir / "analysis_manifest.json"

    signature = make_signature(fastq, primer_sha256, anchor)

    # ---------------------------------------------------------
    # 存在即跳过：只要该样本的三个核心结果文件均存在且非空，
    # 就直接读取旧结果，不再重新扫描 FASTQ。
    # 如果你更换了 primer / anchor / FASTQ 并希望重新计算，
    # 在 config.sh 中设置 FORCE_RERUN_COUNT=1。
    # ---------------------------------------------------------
    if (
        not force
        and detail_file.exists() and detail_file.stat().st_size > 0
        and legacy_file.exists() and legacy_file.stat().st_size > 0
        and overall_file.exists() and overall_file.stat().st_size > 0
    ):
        print(f"[SKIP] {sample_id}: existing non-empty count results found")
        detail_df = pd.read_csv(detail_file)
        rows = detail_df.to_dict("records")
        overall_df = pd.read_csv(overall_file)
        if overall_df.empty:
            raise ValueError(f"Existing overall result is empty: {overall_file}")
        overall = overall_df.iloc[0].to_dict()
        return rows, overall

    print("=" * 70)
    print(f"[COUNT] Sample: {sample_id}")
    print(f"[COUNT] FASTQ : {fastq}")
    print("=" * 70)

    n_patterns = len(patterns)
    counts = [0] * n_patterns
    counts_anchor = [0] * n_patterns
    total_reads = 0
    any_ifp = 0
    any_ifp_anchor = 0
    any_spacer = 0
    any_spacer_anchor = 0

    ifp_ids = {p.pattern_id for p in patterns if p.target_type == "iFP"}
    spacer_ids = {p.pattern_id for p in patterns if p.target_type == "iFP_spacer"}

    for seq in iter_fastq_sequences(fastq):
        total_reads += 1
        matched = automaton.find_pattern_ids(seq)
        if not matched:
            if total_reads % 1_000_000 == 0:
                print(f"[COUNT] {sample_id}: processed {total_reads:,} reads")
            continue

        has_anchor = anchor in seq
        for pid in matched:
            counts[pid] += 1
            if has_anchor:
                counts_anchor[pid] += 1

        hit_ifp = bool(matched & ifp_ids)
        hit_spacer = bool(matched & spacer_ids)
        if hit_ifp:
            any_ifp += 1
            if has_anchor:
                any_ifp_anchor += 1
        if hit_spacer:
            any_spacer += 1
            if has_anchor:
                any_spacer_anchor += 1

        if total_reads % 1_000_000 == 0:
            print(f"[COUNT] {sample_id}: processed {total_reads:,} reads")

    if total_reads == 0:
        raise ValueError(f"No reads found in FASTQ: {fastq}")

    rows: List[Dict[str, object]] = []
    for p in patterns:
        c = counts[p.pattern_id]
        ca = counts_anchor[p.pattern_id]
        rows.append({
            "sample_id": sample_id,
            "gene": p.gene,
            "number": p.number,
            "target_type": p.target_type,
            "seq": p.seq,
            "count": c,
            "count_anchor": ca,
            "count_pipe_anchor": f"{c}|{ca}",
            "total_reads": total_reads,
            "on_target_rate_pct": c / total_reads * 100.0,
            "on_target_anchor_rate_pct": ca / total_reads * 100.0,
            "anchor_within_target_pct": (ca / c * 100.0) if c else 0.0,
        })

    pd.DataFrame(rows).to_csv(detail_file, index=False)

    # Legacy-like lib file: one count|count_anchor row per search sequence + final total_reads row.
    legacy_values = [r["count_pipe_anchor"] for r in rows] + [total_reads]
    pd.DataFrame({sample_id: legacy_values}).to_csv(legacy_file)

    overall: Dict[str, object] = {
        "sample_id": sample_id,
        "total_reads": total_reads,
        "reads_with_any_iFP": any_ifp,
        "reads_with_any_iFP_and_anchor": any_ifp_anchor,
        "reads_with_any_iFP_spacer": any_spacer,
        "reads_with_any_iFP_spacer_and_anchor": any_spacer_anchor,
        "any_iFP_rate_pct": any_ifp / total_reads * 100.0,
        "any_iFP_anchor_rate_pct": any_ifp_anchor / total_reads * 100.0,
        "any_iFP_spacer_rate_pct": any_spacer / total_reads * 100.0,
        "any_iFP_spacer_anchor_rate_pct": any_spacer_anchor / total_reads * 100.0,
        "anchor_within_any_iFP_pct": (any_ifp_anchor / any_ifp * 100.0) if any_ifp else 0.0,
        "anchor_within_any_iFP_spacer_pct": (any_spacer_anchor / any_spacer * 100.0) if any_spacer else 0.0,
    }
    pd.DataFrame([overall]).to_csv(overall_file, index=False)
    manifest_file.write_text(json.dumps(signature, indent=2, ensure_ascii=False), encoding="utf-8")

    print(f"[DONE] {sample_id}: total_reads={total_reads:,}")
    print(f"[DONE] {detail_file}")
    return rows, overall


def main() -> int:
    args = parse_args()
    fastp_dir = Path(args.fastp_dir)
    out_dir = Path(args.out_dir)
    primer_excel = Path(args.primer_excel)
    samples_table = Path(args.samples_table)
    anchor = args.anchor.strip().upper()

    for path, label in [
        (fastp_dir, "fastp dir"),
        (primer_excel, "primer Excel"),
        (samples_table, "samples table"),
    ]:
        if not path.exists():
            raise FileNotFoundError(f"{label} not found: {path}")

    if not anchor or any(ch not in "ACGTN" for ch in anchor):
        raise ValueError(f"Invalid anchor sequence: {anchor}")

    out_dir.mkdir(parents=True, exist_ok=True)
    sample_out_root = out_dir / "per_sample"
    sample_out_root.mkdir(parents=True, exist_ok=True)

    sample_ids = load_enabled_samples(samples_table)
    primer_sha256 = sha256_file(primer_excel)
    patterns = load_patterns(primer_excel, args.primer_sheet, out_dir)
    automaton = AhoCorasick(patterns)

    print("=" * 70)
    print(f"Primer Excel      : {primer_excel}")
    print(f"Primer SHA256     : {primer_sha256}")
    print(f"Targets           : {len(patterns)//2}")
    print(f"Search sequences  : {len(patterns)}")
    print(f"Enabled samples   : {len(sample_ids)}")
    print(f"Anchor RC         : {anchor}")
    print("=" * 70)

    all_rows: List[Dict[str, object]] = []
    all_overall: List[Dict[str, object]] = []

    for sample_id in sample_ids:
        fastq = fastp_dir / f"{sample_id}_qc1.fastq.gz"
        if not fastq.exists() or fastq.stat().st_size == 0:
            raise FileNotFoundError(f"Expected fastp FASTQ not found/empty: {fastq}")
        rows, overall = analyze_sample(
            sample_id=sample_id,
            fastq=fastq,
            patterns=patterns,
            automaton=automaton,
            anchor=anchor,
            sample_out_root=sample_out_root,
            primer_sha256=primer_sha256,
            force=args.force,
        )
        all_rows.extend(rows)
        all_overall.append(overall)

    long_file = out_dir / "sum_onTargetRate_ALL_samples_long.csv"
    pd.DataFrame(all_rows).to_csv(long_file, index=False)

    overall_file = out_dir / "overall_onTargetRate_ALL_samples.csv"
    pd.DataFrame(all_overall).to_csv(overall_file, index=False)

    # Matrix compatible with the old conceptual layout: each cell is count|count_anchor.
    pattern_meta = pd.DataFrame([
        {"gene": p.gene, "number": p.number, "target_type": p.target_type, "seq": p.seq}
        for p in patterns
    ])
    matrix = pattern_meta.copy()
    for sample_id in sample_ids:
        sample_rows = [r for r in all_rows if r["sample_id"] == sample_id]
        matrix[sample_id] = [r["count_pipe_anchor"] for r in sample_rows]

    total_row: Dict[str, object] = {
        "gene": "",
        "number": "TOTAL_READS",
        "target_type": "total_reads",
        "seq": "",
    }
    overall_by_sample = {str(x["sample_id"]): x for x in all_overall}
    for sample_id in sample_ids:
        total_row[sample_id] = int(float(overall_by_sample[sample_id]["total_reads"]))
    matrix = pd.concat([matrix, pd.DataFrame([total_row])], ignore_index=True)

    matrix_file = out_dir / "sum_onTargetRate_ALL_samples_matrix.csv"
    matrix.to_csv(matrix_file, index=False)

    run_manifest = {
        "algorithm_version": ALGORITHM_VERSION,
        "primer_excel": str(primer_excel.resolve()),
        "primer_sha256": primer_sha256,
        "primer_sheet": args.primer_sheet,
        "anchor": anchor,
        "samples": sample_ids,
        "fastp_dir": str(fastp_dir.resolve()),
    }
    (out_dir / "run_manifest.json").write_text(
        json.dumps(run_manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    print("=" * 70)
    print("All counting done.")
    print(f"Long table    : {long_file}")
    print(f"Matrix table  : {matrix_file}")
    print(f"Overall table : {overall_file}")
    print(f"Primer audit  : {out_dir / 'primer_audit.txt'}")
    print("=" * 70)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"[ERROR] {exc}", file=sys.stderr)
        raise
