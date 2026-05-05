#!/usr/bin/env python3
"""Apply ID mapping from reannotate_gtf.py to an LRAA .quant.expr file.

Replaces the gene_id and transcript_id columns in the quant table
using the mapping TSV produced by reannotate_gtf.py --mapping.

Columns in quant.expr (tab-separated, no header in some LRAA versions):
    gene_id, transcript_id, uniq_reads, all_reads, isoform_fraction,
    unique_gene_read_fraction, TPM, exons, introns, splice_hash_code
"""

import argparse
import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--quant",
        required=True,
        help="LRAA .quant.expr file to reannotate",
    )
    parser.add_argument(
        "--mapping",
        required=True,
        help="ID mapping TSV from reannotate_gtf.py --mapping",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Path to write the reannotated .quant.expr file",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Read the quant table — LRAA writes a header line starting with gene_id
    quant = pd.read_csv(args.quant, sep="\t")

    # Read the ID mapping
    mapping = pd.read_csv(args.mapping, sep="\t")

    # Build lookup dicts
    tx_map = dict(zip(mapping["old_transcript_id"], mapping["new_transcript_id"]))
    gene_map = dict(zip(mapping["old_gene_id"], mapping["new_gene_id"]))

    # Apply mapping (keep original ID if no mapping exists)
    quant["transcript_id"] = (
        quant["transcript_id"].map(tx_map).fillna(quant["transcript_id"])
    )
    quant["gene_id"] = quant["gene_id"].map(gene_map).fillna(quant["gene_id"])

    quant.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    main()
