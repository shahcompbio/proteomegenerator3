# Multi-Assembler Merge Plan

## Overview

Add support for running multiple long-read assemblers (bambu, lraa, stringtie) and merging their GTF outputs into a single consensus assembly. This replaces the current single-assembler-at-a-time approach.

## Execution Order

```
per-tool GTFs (assembly_ch)
    │
    ├──► STRINGTIE_MERGE (no ref GTF) ──► merged GTF
    │         │
    │         ▼
    │    REANNOTATEGTF (merged GTF) ──► reannotated merged GTF
    │
    ├──► GFFCOMPARE (queries=per-tool GTFs, ref=reannotated merged GTF) ──► .tracking provenance
    │
    └──► LRAA_SQANTI (reannotated merged GTF vs ref annotation) ──► classification
              │
              ▼
         GFFREAD ──► PREDICT_ORFS ──► FASTA_MERGE_ANNOTATE (simplified)
```

## Steps

### Step 1: Tool name in GTF filenames

Each assembler's output GTF should be named `{meta.id}.{meta.tool}.gtf` — e.g., `A.stringtie_lr.gtf`, `A.bambu.gtf`, `A.lraa.gtf`. This ensures GFFcompare's `.tracking` file references are human-readable by tool.

**Changes:**

- Update `ext.prefix` in `conf/modules.config` for `REANNOTATEGTF` to `"${meta.id}.${meta.tool}"`
- Ensure bambu's output GTF follows the same `{id}.{tool}` naming convention (bambu doesn't go through `REANNOTATEGTF` — it emits GTFs from `BAMBU_FILTER`)

### Step 2: STRINGTIE_MERGE + REANNOTATEGTF

Group `assembly_ch` by `meta.id`, collect the per-tool GTFs, and run `STRINGTIE_MERGE` **without** the reference annotation (to avoid pulling in undetected transcripts).

Immediately after merging, run `REANNOTATEGTF` on the merged GTF to assign:

- Reference ENST/ENSG IDs for canonical transcripts (exact matches)
- Novel prefixes for new transcripts

This reannotated merged GTF is the single GTF used by all downstream steps.

### Step 3: GFFcompare — single run

Run `GFFCOMPARE` once with:

- **Query**: all per-tool GTFs collected
- **Reference** (`-r`): the reannotated merged GTF

This produces a single `.tracking` file mapping each merged transcript back to which tool(s) contributed it. Builds a provenance table: `merged_transcript_id → [bambu, lraa, stringtie]`.

### Step 4: SQANTI classification

Reuse the existing `LRAA_SQANTI` module (`modules/local/lraa/sqanti/main.nf`) — no new module needed. Alias it as `SQANTI_CLASSIFY` for clarity. Run on the reannotated merged GTF against the reference annotation to classify all assembled isoforms (FSM, ISM, NIC, NNC, etc.).

### Step 5: Simplify downstream (GFFREAD → PREDICT_ORFS → FASTA_MERGE_ANNOTATE)

With a single merged+reannotated GTF, the downstream simplifies:

- **GFFREAD**: run once on the single merged GTF (instead of per-tool)
- **PREDICT_ORFS**: run once on the single merged cDNA FASTA
- **FASTA_MERGE_ANNOTATE**: simplifies significantly:
  - Remove tool-branching logic (no more `meta.tool` checks)
  - Remove `CAT_CAT_SAMPLES` for short-read + long-read merging (already merged at GTF level)
  - Single FASTA path: `TRANSDECODER2FASTA` → optional fusion concat → `SEQKIT_RMDUP` → `SEQKIT_STATS`

### Step 6: Refactor NDR to single value

Simplify bambu to run at **one** NDR value only:

- If user provides `--NDR`, use that value
- If user provides `--recommended_NDR`, use Bambu's recommended NDR (`"DEFAULT"`)
- Remove the current logic that creates a channel of both when both flags are set
- `ch_NDR` is always a single value, not a multi-value channel

## File Changes Summary

| File                                              | Change                                                                                                                                                                        |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `workflows/proteomegenerator3.nf`                 | Simplify NDR to single value (L74-82), complete merge block (L151) with STRINGTIE_MERGE → REANNOTATEGTF → GFFCOMPARE → LRAA_SQANTI, simplify downstream to single merged path |
| `conf/modules.config`                             | Update REANNOTATEGTF prefix to `${meta.id}.${meta.tool}`, add process configs for STRINGTIE_MERGE_ASSEMBLERS and SQANTI_CLASSIFY                                              |
| `subworkflows/local/bam_assembly_bambu/main.nf`   | Ensure output GTF naming matches `{id}.{tool}` convention                                                                                                                     |
| `subworkflows/local/fasta_merge_annotate/main.nf` | Remove tool-branching logic, simplify to single-input path                                                                                                                    |
| `bin/transdecoder2fasta.py`                       | Handle merged transcript ID prefixes from reannotated merged GTF                                                                                                              |
| `nextflow_schema.json`                            | Update NDR/assembler parameter descriptions                                                                                                                                   |

## Key Modules Already Available

- `STRINGTIE_MERGE`: `modules/nf-core/stringtie/merge/main.nf`
- `GFFCOMPARE`: `modules/nf-core/gffcompare/main.nf`
- `LRAA_SQANTI`: `modules/local/lraa/sqanti/main.nf`
- `REANNOTATEGTF`: `modules/local/reannotategtf/main.nf`
- `reannotate_gtf.py`: supports `--tool` flag with prefix mappings (currently stringtie/lraa)
