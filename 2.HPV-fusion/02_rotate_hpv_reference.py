#!/usr/bin/env python3
"""Move the linear cut of each HPV reference away from its original origin."""

import argparse

from fusion_common import iter_fasta, wrap_fasta


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("reference")
    args = parser.parse_args()
    for accession, sequence in iter_fasta(args.reference):
        midpoint = len(sequence) // 2
        print(">%s\n%s" % (accession, wrap_fasta(sequence[midpoint:] + sequence[:midpoint])))


if __name__ == "__main__":
    main()
