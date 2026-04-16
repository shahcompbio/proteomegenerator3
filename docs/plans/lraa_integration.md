# Plan: Integrate LRAA as an alternative long-read assembler

Target branch: `LRAA_stringtie` (branched from `dev`).

## 1. Goal & scope

Add [LRAA](https://github.com/MethodsDev/LongReadAlignmentAssembler) as a
user-selectable alternative to **Bambu** for long-read transcript assembly and
quantification. Downstream stages (`GFFREAD` → `PREDICT_ORFS` →
`FASTA_MERGE_ANNOTATE`) are unchanged; LRAA only needs to emit a GTF tagged
with `meta.tool = 'lraa'`, matching Bambu's output contract.

Selection is **either/or** via a new `--long_read_assembler` enum
(`bambu` | `lraa`). Running both assemblers in a single invocation is out of
scope for v1.

## 2. Background: how LRAA differs from Bambu

Key differences that drive the design:

| Aspect | Bambu | LRAA |
|---|---|---|
| Multi-sample call | Single R call ingests all `rc_file`s, emits unified GTF + counts | No. One BAM per invocation. |
| Merge utility | Built-in (`mode = "multi-sample"`) | External: `util/merge_LRAA_GTFs.py` |
| Cohort quant | Built into merge call | Per-sample re-run with `--quant_only` against merged GTF, then external matrix assembly |
| Input | Read-class RDS (pre-computed) or BAM | BAM directly (no read-class step) |
| ONT-specific flag | N/A | **Omit `--HiFi`** on merge for ONT/LowFi mode |
| NDR analog | `--NDR` | None; uses own isoform-discovery filters |
| Shipped multi-sample orchestration | Bambu R package handles it | Only single-sample WDL; caller scripts the loop |

**Implication:** the LRAA subworkflow needs three stages — per-sample
discovery → single merge → per-sample re-quant — plus a matrix-assembly step,
whereas Bambu fits in one or two.

## 3. User-facing parameters

Add to `nextflow.config` and `nextflow_schema.json`:

| Param | Type | Default | Notes |
|---|---|---|---|
| `long_read_assembler` | enum(`bambu`,`lraa`) | `bambu` | Top-level switch |
| `skip_lraa_discovery` | bool | `false` | Resume at merge using pre-computed per-sample GTFs from the samplesheet |
| `lraa_min_isoform_fraction` | float | TBD from `LRAA --help` | Expose key LRAA filter(s); finalize list during implementation |

Bambu-specific params (`NDR`, `recommended_NDR`, `yieldsize`, `min_lr_cts`)
stay as-is and are only consumed when `long_read_assembler == 'bambu'`.

## 4. Samplesheet extension

Add a new `filetype` value: `lraa_gtf`. Valid only with
`sequence_type: long_read`.

```csv
subject_id,sample_id,sequence_type,filetype,filepath
PATIENT1,SAMPLE1,long_read,bam,/path/to/s1.bam
PATIENT1,SAMPLE1,long_read,lraa_gtf,/path/to/s1.LRAA.gtf
```

Validation rules:

- `lraa_gtf` filetype is accepted only when `long_read_assembler == 'lraa'`.
- When `skip_lraa_discovery == true`, every sample must have an `lraa_gtf` row.
- BAM is still required even when discovery is skipped, because the
  re-quant step runs against the pipeline-produced merged GTF (not the user's).
- `lraa_quant` (`.quant.expr`) passthrough is **out of scope for v1** — the
  pipeline's merged GTF won't match a user's pre-computed one, so their
  per-sample quant would be inconsistent with cohort merge.

Update `assets/schema_input.json` accordingly.

## 5. New subworkflow: `subworkflows/local/bam_assembly_lraa/`

Mirrors the shape and emit contract of `bam_assembly_bambu`.

### Signature

```groovy
workflow BAM_ASSEMBLY_LRAA {
    take:
    bam_ch            // channel: [ val(meta), path(bam) ]    — post-PREPROCESS_READS
    lraa_gtf_ch       // channel: [ val(meta), path(gtf) ]    — pre-computed, empty when skip_lraa_discovery=false
    skip_multisample  // val
    skip_discovery    // val (params.skip_lraa_discovery)
    sample_count      // val
    ref_gtf_ch        // channel (guided mode reference GTF)
    ref_fasta         // val   (params.fasta)

    main:
    ...

    emit:
    gtf         = ...  // channel: [ val(meta), path(gtf) ]   — matches BAM_ASSEMBLY_BAMBU.out.gtf
    quant_matrix = ... // channel: [ val(meta), path(matrix_tsv) ]  optional, multi-sample only
    versions    = ch_versions
}
```

### Control flow

```
BAM_ASSEMBLY_LRAA
│
├── per-sample GTF source:
│     if skip_discovery: use lraa_gtf_ch from samplesheet
│     else:              LRAA_ASSEMBLY(bam_ch, ref_gtf_ch, ref_fasta)  [guided mode]
│
├── single-sample / skip_multisample branch:
│     emit per-sample GTF directly (meta stays per-sample)
│
└── multi-sample branch:
      LRAA_MERGE_GTFS(collect(per_sample_gtfs), ref_fasta)      // util/merge_LRAA_GTFs.py, NO --HiFi
        → cohort GTF with meta.id = "merge"
      LRAA_QUANT(bam_ch.combine(cohort_gtf), ref_fasta)         // LRAA --quant_only per sample
        → per-sample .quant.expr
      LRAA_QUANT_MERGE(collect(per_sample_expr))                 // Python glue → cohort matrix TSV
        → cohort expression matrix (emitted for downstream consumers)
      emit cohort GTF
```

Notes:
- `PREPROCESS_READS` still runs (MAPQ / read-length filtering). Its `rc_ch`
  output is unused in the LRAA path; `bam_ch` is what LRAA consumes.
- No NDR channel plumbing in the LRAA branch.
- Matrix assembly is a thin Python helper in `bin/` (e.g.
  `merge_lraa_quant.py`) — joins per-sample `.quant.expr` on transcript_id,
  emits a single TSV (transcript_id × sample_id counts/TPM).

## 6. New modules under `modules/local/lraa/`

Each module follows PG3 convention: `main.nf`, `environment.yml`,
`meta.yml`, `tests/main.nf.test`, plus a `stub` block.

### 6.1 `lraa/assembly/main.nf` — `LRAA_ASSEMBLY`

```groovy
process LRAA_ASSEMBLY {
    tag     "${meta.id}"
    label   'process_high_memory'
    container "<official LRAA image — confirm tag during implementation>"

    input:
    tuple val(meta), path(bam), path(ref_gtf, stageAs: 'reference.gtf')
    path ref_genome

    output:
    tuple val(meta), path("${meta.id}.LRAA.gtf"),        emit: gtf
    tuple val(meta), path("${meta.id}.quant.expr"),      optional: true, emit: quant
    tuple val(meta), path("${meta.id}.LRAA.log"),        emit: log
    path "versions.yml",                                 emit: versions

    script:
    def args = task.ext.args ?: ''
    """
    LRAA \\
        --genome ${ref_genome} \\
        --gtf    reference.gtf \\
        --bam    ${bam} \\
        --output_prefix ${meta.id} \\
        --CPU ${task.cpus} \\
        ${args} 2>&1 | tee ${meta.id}.LRAA.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """
}
```

(Guided mode: both `--gtf` and `--genome` supplied. Confirm exact flag names
and output file-naming during implementation by inspecting the container's
`LRAA --help`.)

### 6.2 `lraa/merge/main.nf` — `LRAA_MERGE_GTFS`

```groovy
process LRAA_MERGE_GTFS {
    tag     "merge"
    label   'process_medium'
    container "<same official LRAA image>"

    input:
    path(gtfs, stageAs: 'sample_gtfs/*')
    path ref_genome

    output:
    tuple val([id: 'merge']), path("cohort.LRAA.merged.gtf"),         emit: gtf
    path("cohort.LRAA.merged.gtf.tracking.tsv"),                       emit: tracking
    path "versions.yml",                                               emit: versions

    script:
    """
    merge_LRAA_GTFs.py \\
        --genome ${ref_genome} \\
        --gtf sample_gtfs/*.gtf \\
        --output_gtf cohort.LRAA.merged.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """
}
```

Critical: **do not pass `--HiFi`** — omitted for ONT/LowFi mode per LRAA
docs.

### 6.3 `lraa/quant/main.nf` — `LRAA_QUANT`

```groovy
process LRAA_QUANT {
    tag     "${meta.id}"
    label   'process_high'
    container "<same official LRAA image>"

    input:
    tuple val(meta), path(bam), path(cohort_gtf, stageAs: 'cohort.gtf')
    path ref_genome

    output:
    tuple val(meta), path("${meta.id}.quant.expr"),      emit: quant
    tuple val(meta), path("${meta.id}.quant.tracking"),  optional: true, emit: tracking
    path "versions.yml",                                 emit: versions

    script:
    """
    LRAA \\
        --quant_only \\
        --genome ${ref_genome} \\
        --gtf    cohort.gtf \\
        --bam    ${bam} \\
        --output_prefix ${meta.id} \\
        --CPU ${task.cpus}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """
}
```

### 6.4 `lraa/quant_merge/main.nf` — `LRAA_QUANT_MERGE`

Thin wrapper calling a new `bin/merge_lraa_quant.py`. Pure-Python pandas
join over per-sample `.quant.expr` files on `transcript_id`; emits
`cohort.expr_matrix.tsv`. Runs in the base conda env (no LRAA container
needed).

### 6.5 Container

Use the **official LRAA image**. Pin a specific tag (not `latest`). The
exact image path/tag is TBD — check the project's DockerHub or quay.io
namespace during implementation. Hardcoded in module `container` directives
as with `quay.io/shahlab_singularity/bambu:3.10.0beta`.

## 7. Samplesheet parsing helpers

In `subworkflows/local/utils_nfcore_proteomegenerator3_pipeline/main.nf`,
add a helper parallel to `getLongReadRcFiles`:

```groovy
def getLongReadLraaGtfs(ch_samplesheet) {
    ch_samplesheet
        .filter { row -> row.sequence_type == 'long_read' && row.filetype == 'lraa_gtf' }
        .map    { row -> [[id: row.sample_id, subject_id: row.subject_id], file(row.filepath)] }
}
```

## 8. Wiring in `workflows/proteomegenerator3.nf`

Replace the current unconditional Bambu call with a branch:

```groovy
ch_long_read_lraa_gtfs = getLongReadLraaGtfs(ch_samplesheet)

if (params.long_read_assembler == 'bambu') {
    BAM_ASSEMBLY_BAMBU(
        rc_ch, params.skip_multisample, sample_count,
        ch_NDR, ref_gtf_ch, bam_ch,
    )
    ch_versions = ch_versions.mix(BAM_ASSEMBLY_BAMBU.out.versions)
    assembly_ch = BAM_ASSEMBLY_BAMBU.out.gtf
        .map { meta, gtf -> [meta + [tool: 'bambu'], gtf] }
}
else if (params.long_read_assembler == 'lraa') {
    BAM_ASSEMBLY_LRAA(
        bam_ch,
        ch_long_read_lraa_gtfs,
        params.skip_multisample,
        params.skip_lraa_discovery,
        sample_count,
        ref_gtf_ch,
        params.fasta,
    )
    ch_versions = ch_versions.mix(BAM_ASSEMBLY_LRAA.out.versions)
    assembly_ch = BAM_ASSEMBLY_LRAA.out.gtf
        .map { meta, gtf -> [meta + [tool: 'lraa'], gtf] }
}
```

Everything downstream of `assembly_ch` is untouched — `GFFREAD`,
`PREDICT_ORFS`, `FASTA_MERGE_ANNOTATE`, fusion handling, MultiQC.

**Fusion compatibility:** LRAA has no fusion-aware mode, but PG3 handles
fusions in a separate branch that joins at `FASTA_MERGE_ANNOTATE`, so
`--fusions` works transparently with `--long_read_assembler=lraa`.

## 9. Config & profiles

- `nextflow.config`: add new params with defaults.
- `conf/base.config`: verify `process_high_memory` is adequate for LRAA; add
  a dedicated resource block if profiling shows different footprint.
- New `conf/test_lraa.config` mirroring `conf/test_single_sample.config`
  but with `long_read_assembler = 'lraa'`.
- Add a `test_lraa` profile to `nextflow.config`.

## 10. Tests (nf-test)

1. `modules/local/lraa/*/tests/main.nf.test` — one per module, stub + real.
2. `subworkflows/local/bam_assembly_lraa/tests/main.nf.test`:
   - single-sample (discovery on)
   - multi-sample (discovery on, full pipeline)
   - multi-sample with `skip_lraa_discovery = true` using pre-computed GTFs
3. Top-level `tests/main.nf.test` variant running
   `-profile test_lraa,docker`.
4. Update `tests/nextflow.config` with any LRAA test-data URLs.
5. Confirm snapshot outputs stay deterministic (LRAA may need a random
   seed flag — check).

## 11. Documentation

- `docs/usage.md`: new section on `--long_read_assembler`; document
  `skip_lraa_discovery` resume path; call out Bambu-only vs LRAA-only
  samplesheet rows (`rc_file` vs `lraa_gtf`).
- `docs/output.md`: document LRAA output dirs (`lraa/`, `lraa_merged/`,
  `lraa_quant/`, `cohort.expr_matrix.tsv`).
- `CLAUDE.md`: update Architecture section to list LRAA alongside Bambu;
  update Configuration → Important Parameters.
- `CITATIONS.md`: add LRAA citation.
- `README.md`: one-line mention in the features list.

## 12. Suggested PR breakdown

Keep PRs small and independently reviewable:

1. **PR 1 — LRAA modules.** `modules/local/lraa/{assembly,merge,quant,quant_merge}` with stubs, unit tests, and the `bin/merge_lraa_quant.py` helper. No wiring yet. Mergeable standalone.
2. **PR 2 — `bam_assembly_lraa` subworkflow.** Plus its nf-test suite. Consumes PR 1 modules. Still not wired into main workflow.
3. **PR 3 — Wire into main workflow + schema + samplesheet validation + `test_lraa` profile.** End-to-end runnable.
4. **PR 4 — Docs.**

All PRs target `dev`.

## 13. Out of scope / future work

- De novo (non-guided) LRAA mode.
- `lraa_quant` samplesheet filetype for matrix-only resume.
- Running Bambu and LRAA in parallel and consuming both into
  `assembly_ch` (the current code shape tolerates this, but the UX is
  either/or for v1).
- LRAA-specific MultiQC module (none exists upstream as of writing).
- Single-cell / cluster-guided LRAA modes — PG3 is bulk-focused.

## 14. Open items to resolve during implementation

- [ ] Confirm official LRAA Docker image URI and a stable tag.
- [ ] Confirm exact CLI flag names by running `LRAA --help` inside the container (the sketch above uses plausible names — validate before coding).
- [ ] Decide which LRAA isoform-discovery filters to expose as pipeline params (keep the surface minimal; prefer `--args` passthrough via `ext.args` for power users).
- [ ] Confirm `merge_LRAA_GTFs.py` is on `$PATH` in the official image or needs a full path.
- [ ] Profile memory/CPU on a representative ONT BAM to size resource labels.
