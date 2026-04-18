#!/usr/bin/env python
"""
Re-annotate a gffcompare-annotated GTF so that:
  - exact-match transcripts (class_code '=') recover their reference ENST / ENSG IDs
  - novel transcripts get a tool-specific prefix (e.g. StrgTx / LraaTx)

Supports --tool stringtie (default) and --tool lraa.
"""
from gtfparse import read_gtf
import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import logging

logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)
import polars
from pathlib import Path
import typing as t
from typing import Union
import argparse

COMMONS_COL = [
    "seqname",
    "source",
    "feature",
    "start",
    "end",
    "score",
    "strand",
    "frame",
]

# Tool → (novel transcript prefix, novel gene prefix)
TOOL_PREFIXES = {
    "stringtie": ("StrgTx", "StrgGene"),
    "lraa": ("LraaTx", "LraaGene"),
}


def write_gtf(
    df: polars.DataFrame,
    export_path: Union[str, Path],
    tool: str,
    headers: t.List[str] = None,
):
    headers = headers or []
    with open(export_path, "w") as f:
        f.write(f"# re-annotated gtf of merged transcripts from {tool}\n")
        f.write(f"###\n")
        for header in headers:
            f.write(f"{header}\n")
        for _, row in df.iterrows():
            f.write(f"{commons_cols(row)}\t{custom_fields(row)}\n")


def commons_cols(row) -> str:
    return "\t".join([str(row[field] or ".") for field in COMMONS_COL])


def custom_fields(row) -> str:
    return "; ".join(
        [
            f'{field} "{row[field]}"'
            for field in row.keys()
            if (field not in COMMONS_COL) and (row[field])
        ]
    )


p = argparse.ArgumentParser(description="Re-annotate gffcompare GTF with reference IDs")
p.add_argument("gffcmp_results", help="gffcompare-annotated GTF")
p.add_argument(
    "--reference_fai",
    help="reference FASTA index (.fai) for filtering out transcripts with exon boundaries that exceed contig span",
)
p.add_argument("output_file", help="output re-annotated GTF")
p.add_argument(
    "--tool",
    choices=list(TOOL_PREFIXES.keys()),
    default="stringtie",
    help="assembler tool (determines novel transcript/gene prefix)",
)
p.add_argument(
    "--mapping",
    default=None,
    help="optional output TSV mapping old IDs to new IDs",
)
args = p.parse_args()

tx_prefix, gene_prefix = TOOL_PREFIXES[args.tool]

gffcmp = read_gtf(args.gffcmp_results)

annotatedat = pd.DataFrame()
id_mapping_rows = []
df = gffcmp
## nans give an error when doing ballgown estimates
df["strand"] = df["strand"].apply(lambda x: x if x in ["+", "-"] else ".")
## filter out any transcript structures with exon boundaries that EXCEED the span of the contig (e.g. due to a gffcompare bug)
reference_fai = pd.read_csv(
    args.reference_fai,
    sep="\t",
    header=None,
    names=["seqname", "length", "offset", "linebases", "linewidth"],
)
reference_lengths = dict(zip(reference_fai["seqname"], reference_fai["length"]))
df = df[
    df.apply(lambda row: row["end"] <= reference_lengths.get(row["seqname"], 0), axis=1)
]
grouped = df.groupby("transcript_id")

# Iterate through each group
i = 1
# track transcript ids to check for duplicates
tx_ids = set()
for transcript_id, group_df in grouped:
    old_gene_id = list(group_df["gene_id"])[0]
    ## rename gene
    ref_gene = list(group_df["ref_gene_id"])[0]
    if not ref_gene == "":
        group_df["gene_id"] = ref_gene
    else:
        group_df["gene_id"] = "%s%d" % (gene_prefix, i)
    ### rename transcripts so exact matches are given known ids and those which are not
    ### are given a tool-specific novel id
    class_code = list(group_df["class_code"])[0]
    tx_id = list(group_df["cmp_ref"])[0]
    if class_code == "=":
        group_df["transcript_id"] = tx_id
    else:
        tx_id = "%s%d" % (tx_prefix, i)
        group_df["transcript_id"] = tx_id

    new_gene_id = list(group_df["gene_id"])[0]
    new_tx_id = list(group_df["transcript_id"])[0]
    id_mapping_rows.append(
        {
            "old_transcript_id": transcript_id,
            "new_transcript_id": new_tx_id,
            "old_gene_id": old_gene_id,
            "new_gene_id": new_gene_id,
            "class_code": class_code,
        }
    )
    # check for duplicate transcript IDs, which can occur due to multiple
    # transcripts having exact intron matches to the same reference transcript
    # in this case, we keep the first one and skip the rest, but print a warning
    if tx_id in tx_ids:
        print(
            f"Warning: transcript ID {tx_id} already seen, likely to multiple exact intron matches to ref transcript.\n"
            f"Check gffcompare results for transcript {transcript_id}."
        )
    else:
        tx_ids.add(tx_id)
        annotatedat = pd.concat([annotatedat, group_df])
    i = i + 1

write_gtf(annotatedat, args.output_file, args.tool)

if args.mapping:
    mapping_df = pd.DataFrame(id_mapping_rows)
    mapping_df.to_csv(args.mapping, sep="\t", index=False)
