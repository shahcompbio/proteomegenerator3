#!/usr/bin/env python
from gtfparse import read_gtf
import warnings

warnings.filterwarnings("ignore")
import pandas as pd
import logging

logger = logging.getLogger()
logger.setLevel(logging.CRITICAL)
from tqdm import tqdm
import polars
from pathlib import Path
import typing as t
from typing import Union
import sys

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


def write_gtf(
    df: polars.DataFrame, export_path: Union[str, Path], headers: t.List[str] = None
):
    headers = headers or []
    with open(export_path, "w") as f:
        f.write(f"# re-annotated gtf of merged transcripts from StringTie\n")
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


## inputs
gffcmp_results = sys.argv[1]
output_file = sys.argv[2]

gffcmp = read_gtf(gffcmp_results)

annotatedat = pd.DataFrame()
df = gffcmp
## nans give an error when doing ballgown estimates
df["strand"] = df["strand"].apply(lambda x: x if x in ["+", "-"] else ".")
## a test
print(set(df["strand"]))
grouped = df.groupby("transcript_id")

# Get the total number of groups for tqdm
total_groups = len(grouped)
# Iterate through each group
i = 1
for transcript_id, group_df in tqdm(grouped, total=total_groups):
    ## rename gene
    ref_gene = list(group_df["ref_gene_id"])[0]
    if not ref_gene == "":
        group_df["gene_id"] = ref_gene
    else:
        group_df["gene_id"] = "StrgGene%d" % i
    ### rename transcripts so exact matches are given known ids and those which are not
    ### are given a stringtie id
    class_code = list(group_df["class_code"])[0]
    tx_id = list(group_df["cmp_ref"])[0]
    if class_code == "=":
        group_df["transcript_id"] = tx_id
    else:
        group_df["transcript_id"] = "StrgTx%d" % i
    annotatedat = pd.concat([annotatedat, group_df])
    i = i + 1

write_gtf(annotatedat, output_file)
