#!/usr/bin/env python3
"""Merge per-sample LRAA .quant.expr files into a cohort expression matrix.

Each input file is a tab-separated table produced by LRAA with columns:
    gene_id, transcript_id, uniq_reads, all_reads, isoform_fraction,
    unique_gene_read_fraction, TPM, exons, introns, splice_hash_code

The output is a wide-format TSV with one row per transcript and columns:
    gene_id, transcript_id, <sample1>_TPM, <sample2>_TPM, ...

Sample names are derived from the file stem (everything before '.quant.expr').
"""

import argparse
import os
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        nargs="+",
        required=True,
        help="One or more LRAA .quant.expr files",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the merged expression matrix TSV",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    frames = []
    for fpath in args.input:
        sample = Path(fpath).name.replace(".quant.expr", "")
        df = pd.read_csv(fpath, sep="\t")
        df = df[["gene_id", "transcript_id", "TPM"]].rename(
            columns={"TPM": f"{sample}_TPM"}
        )
        frames.append(df)

    merged = frames[0]
    for df in frames[1:]:
        merged = merged.merge(df, on=["gene_id", "transcript_id"], how="outer")

    merged = merged.fillna(0).sort_values(["gene_id", "transcript_id"])
    merged.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
