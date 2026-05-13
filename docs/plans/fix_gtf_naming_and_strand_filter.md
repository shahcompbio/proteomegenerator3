# Fix GTF Naming and Unstranded Transcript Filtering

## Problem

1. **GTF naming**: The union reannotated GTF is named `union.union.reannotated.gtf` because `meta.id` is hardcoded to `"union"` in `gtf_merge_annotate/main.nf` and `ext.prefix` in `modules.config` appends `.union` again.
2. **SQANTI failure**: LRAA SQANTI crashes with `AssertionError: orient must be set to + or -` because some transcripts have strand `.` (unstranded) in the reannotated GTF.

## Solution

### 1. Dynamic `meta.id` in `subworkflows/local/gtf_merge_annotate/main.nf`

Derive `meta.id` from the incoming `assembly_ch`:
- If incoming `meta.id == "merge"` (multi-sample path) → use `"cohort"`
- Otherwise keep the original `meta.id` (which is the `sample_id` in single-sample mode)

Result:
- Multi-sample: `cohort.union.reannotated.gtf`
- Single-sample: `{sample_id}.union.reannotated.gtf`

### 2. Filter unstranded transcripts in `bin/reannotate_gtf.py`

After strand normalization (which sets non-standard strands to `.`), filter out any `transcript_id` that has at least one feature with strand `.`. Print a warning with the count of removed transcripts.

## Files Changed

- `subworkflows/local/gtf_merge_annotate/main.nf`
- `bin/reannotate_gtf.py`
