# Plan: Clean up output directory structure

## TL;DR

Fix empty NDR directories in `proteome/`, reorganize `assembly/union/` into top-level `union_transcriptome/` with proper subdirectories, and label the union GFFREAD output clearly. All other intermediate output directories kept as-is.

## Steps

### Step 1: Consolidate TRANSDECODER2FASTA publishDir (fixes empty NDR dirs)

Replace the 4-entry publishDir array (~L49-90 of modules.config) with a single entry using a dynamic path closure that branches on `meta.tool`. Produces exactly one directory per task — no phantom `proteome/NDR_0.1/` or `proteome/NDR_DEFAULT/`.

### Step 2: Consolidate MERGEFUSIONS/FUSIONFASTA publishDir (parallel with step 1)

Replace the 3-entry array (~L27-45) with a single dynamic path branching on `params.NDR` / `params.recommended_NDR`.

### Step 3: Reorganize `assembly/union/` → top-level `union_transcriptome/`

- `GTF_MERGE_ANNOTATE:STRINGTIE_MERGE` → `publishDir = [enabled: false]`
- `GTF_MERGE_ANNOTATE:GFFCOMPARE` → `publishDir = [enabled: false]`
- `GTF_MERGE_ANNOTATE:REANNOTATEGTF` → `union_transcriptome/`
- `GFFCOMPARE_PROVENANCE` → `union_transcriptome/gffcompare/`
- `GTF_MERGE_ANNOTATE:LRAA_SQANTI` → `union_transcriptome/sqanti/`

### Step 4: Tag union meta and fix GFFREAD naming

- In `gtf_merge_annotate/main.nf`: add `tool: 'union'` to the emitted meta
- In modules.config GFFREAD: set `ext.prefix = { meta.tool ? "${meta.id}.${meta.tool}" : meta.id }`
- Result: `{subject_id}.union.fasta` for union path, `{sample_id}.fasta` for single-assembler

## Relevant files

- `conf/modules.config` — publishDir changes
- `subworkflows/local/gtf_merge_annotate/main.nf` — add tool to meta

## Verification

1. `nf-test test` passes
2. Single-assembler: no empty `proteome/NDR_*` dirs; GFFREAD → `{sample_id}.fasta`
3. Multi-assembler: `union_transcriptome/` with `gffcompare/` and `sqanti/` subdirs; GFFREAD → `{subject_id}.union.fasta`
