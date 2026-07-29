# SQANTI3 Integration Plan

## Summary

Replace the lightweight `LRAA_SQANTI` classification with real SQANTI3 v6.0.1 (QC → ML-based Filter → Rescue) to produce curated transcriptomes before ORF prediction. Route all GTFs (single- and multi-assembler) through the subworkflow, renamed to `gtf_merge_sqanti`. Keep LRAA_MERGE for cross-assembler merging. Skip SQANTI3's built-in ORF prediction; continue using Transdecoder+BLAST for proteogenomics FASTA.

## Architecture

```
                     Multi-assembler                Single-assembler
                     ┌────────────┐                ┌────────────┐
                     │ Bambu GTF  │                │  Bambu GTF │
                     │ LRAA GTF   │                └─────┬──────┘
                     │ StringTie  │                      │
                     └─────┬──────┘                      │
                           │                             │
                    ┌──────▼──────┐                      │
                    │ LRAA_MERGE  │                      │
                    │(splice graph│                      │
                    │  collapse)  │                      │
                    └──────┬──────┘                      │
                           │                             │
                    ┌──────▼──────────────────────────────▼──────┐
                    │              GFFCOMPARE                     │
                    │         (annotate vs ref)                   │
                    ├────────────────────────────────────────────-┤
                    │            REANNOTATEGTF                    │
                    │   (assign ENST/ENSG or novel IDs)           │
                    ├────────────────────────────────────────────-┤
                    │         GFFCOMPARE_PROVENANCE               │
                    │   (track which tool contributed each tx)    │
                    ├────────────────────────────────────────────-┤
                    │            SQANTI3_QC                       │
                    │   (classify isoforms, --skipORF)            │
                    ├────────────────────────────────────────────-┤
                    │          SQANTI3_FILTER                     │
                    │   (ML-based artifact removal)               │
                    ├────────────────────────────────────────────-┤
                    │          SQANTI3_RESCUE                     │
                    │   (recover reference transcripts)           │
                    └──────────────────┬────────────────────────-─┘
                                       │
                                       ▼
                              Curated GTF (emit)
                                       │
                                       ▼
                                   GFFREAD → Transdecoder → FASTA
```

## Key Decisions

| Decision         | Choice              | Rationale                                                      |
| ---------------- | ------------------- | -------------------------------------------------------------- |
| SQANTI3 version  | v6.0.1              | Latest stable; uses TD2 for ORF, improved rescue logic         |
| Merge tool       | Keep LRAA_MERGE     | Already does splice graph reconstruction on GTFs directly      |
| SQANTI3 ORF      | Skip (`--skipORF`)  | Transdecoder+BLAST is more rigorous for proteogenomics         |
| Filter type      | ML-based            | More accurate than rules-based; uses SQANTI3 built-in pipeline |
| Final QC re-run  | No                  | Rescue output already includes classification                  |
| ID handling      | `--force_id_ignore` | Transcript IDs are not PacBio-formatted                        |
| Single-assembler | Same subworkflow    | Skip LRAA_MERGE, run annotation + SQANTI3                      |
| Optional inputs  | None for now        | CAGE, polyA, SJ can be added later                             |

## New Parameters

| Parameter               | Default | Description                                                              |
| ----------------------- | ------- | ------------------------------------------------------------------------ |
| `--skip_sqanti3`        | `false` | Skip SQANTI3 QC/filter/rescue (fast mode, uses reannotated GTF directly) |
| `--sqanti3_filter_type` | `'ml'`  | Filter strategy: `'ml'` or `'rules'`                                     |

## New Modules

### `modules/local/sqanti3/qc/main.nf`

- **Container**: `conesalab/sqanti3:v6.0.1` (or custom build)
- **Input**: GTF, reference GTF, reference FASTA
- **Output**: classification.txt, corrected.gtf, junctions.txt
- **Command**: `sqanti3 qc --isoforms <gtf> --refGTF <ref> --refFasta <fa> --skipORF --force_id_ignore`

### `modules/local/sqanti3/filter/main.nf`

- **Input**: classification.txt, corrected GTF
- **Output**: filtered classification, filtered GTF
- **Command**: `sqanti3 filter ml`

### `modules/local/sqanti3/rescue/main.nf`

- **Input**: filtered classification, filtered GTF, reference GTF, reference FASTA
- **Output**: rescued GTF, rescued classification
- **Command**: `sqanti3 rescue`

## Implementation Phases

### Phase 1: Module Creation

1. Create SQANTI3 QC, Filter, Rescue Nextflow modules
2. Build/reference SQANTI3 v6.0.1 Docker container

### Phase 2: Subworkflow Refactoring

3. Rename `gtf_merge_annotate` → `gtf_merge_sqanti`
4. Add conditional merge logic (>1 GTF → LRAA_MERGE; else skip)
5. Wire SQANTI3 modules into the subworkflow
6. Update main workflow to route ALL GTFs through subworkflow

### Phase 3: Configuration

7. Add `--skip_sqanti3` and `--sqanti3_filter_type` parameters
8. Add process configs (publishDir, resources)
9. Update test configs with `--skip_sqanti3 true`

### Phase 4: Cleanup

10. Delete `modules/local/lraa/sqanti/`
11. Remove LRAA_SQANTI config sections
12. Write nf-tests
13. Update documentation

## Future Enhancements

- CAGE peak BED input (`--sqanti3_cage_peak`)
- polyA motif list (`--sqanti3_polya_motif`)
- Short-read SJ.out.tab input (`--sqanti3_short_read_sj`)
- Parallelization with `-n` chunks for large transcriptomes
