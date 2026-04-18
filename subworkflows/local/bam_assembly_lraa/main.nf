// Transcript assembly and quantification with LRAA
// Three-stage pipeline: per-sample discovery → cohort merge → per-sample re-quant

include { LRAA_ASSEMBLY                      } from '../../../modules/local/lraa/assembly/main'
include { LRAA_MERGE                         } from '../../../modules/local/lraa/merge/main'
include { LRAA_QUANT ; LRAA_QUANT as LRAA_REANNOTATEQUANT } from '../../../modules/local/lraa/quant/main'
include { LRAA_QUANTMERGE                    } from '../../../modules/local/lraa/quantmerge/main'
include { LRAA_SQANTI                        } from '../../../modules/local/lraa/sqanti/main'
include { GFFCOMPARE                         } from '../../../modules/nf-core/gffcompare/main'
include { REANNOTATEGTF                      } from '../../../modules/local/reannotategtf/main'

workflow BAM_ASSEMBLY_LRAA {
    take:
    bam_ch // channel: [ val(meta), path(bam) ]    — post-PREPROCESS_READS
    lraa_gtf_ch // channel: [ val(meta), path(gtf) ]    — pre-computed, empty when skip_lraa_discovery=false
    skip_multisample // val
    skip_discovery // val (params.skip_lraa_discovery)
    sample_count // val
    ref_gtf // val: path(ref_gtf)
    ref_fasta // val: path(genome fasta)
    ref_fai // val: path(genome fai)

    main:
    ch_versions = channel.empty()
    lraa_out_ch = channel.empty()
    quant_input = channel.empty()

    //
    // Step 1: Per-sample GTF source
    //
    if (skip_discovery) {
        // Use pre-computed GTFs from the samplesheet
        per_sample_gtfs = lraa_gtf_ch
    }
    else {
        // Run guided assembly per sample
        assembly_input = bam_ch.combine(Channel.of(file(ref_gtf)))
        LRAA_ASSEMBLY(assembly_input, ref_fasta)
        ch_versions = ch_versions.mix(LRAA_ASSEMBLY.out.versions.first())
        per_sample_gtfs = LRAA_ASSEMBLY.out.gtf
    }

    //
    // Step 2: Single-sample or multi-sample path
    //
    if (skip_multisample || sample_count == 1) {
        // Single-sample: emit per-sample GTF directly
        // Run gffcompare + reannotate to match Bambu/StringTie output contract
        GFFCOMPARE(
            per_sample_gtfs,
            [[], [], []],
            [[id: "ref"], ref_gtf],
        )
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)
        REANNOTATEGTF(GFFCOMPARE.out.annotated_gtf, ref_fai)
        ch_versions = ch_versions.mix(REANNOTATEGTF.out.versions)
        lraa_out_ch = REANNOTATEGTF.out.gtf

        // Classify reannotated isoforms with SQANTI-like categories
        LRAA_SQANTI(
            REANNOTATEGTF.out.gtf,
            ref_gtf,
        )
        ch_versions = ch_versions.mix(LRAA_SQANTI.out.versions)

        // re-run quant after reannotation
        if (!skip_discovery) {
            // Re-quantify each sample against the reannotated GTF
            quant_input = bam_ch.combine(
                REANNOTATEGTF.out.gtf.map { _meta, gtf -> gtf }
            )
            LRAA_REANNOTATEQUANT(quant_input, ref_fasta)
            ch_versions = ch_versions.mix(LRAA_REANNOTATEQUANT.out.versions.first())
        }
    }
    else {
        // Multi-sample: merge → reannotate → re-quant → matrix assembly
        merge_input = per_sample_gtfs
            .map { _meta, gtf -> gtf }
            .collect()
            .map { gtfs -> [[id: "merge"], gtfs] }

        LRAA_MERGE(merge_input, ref_fasta)
        ch_versions = ch_versions.mix(LRAA_MERGE.out.versions)

        // Run gffcompare + reannotate on the merged GTF
        GFFCOMPARE(
            LRAA_MERGE.out.gtf,
            [[], [], []],
            [[id: "ref"], ref_gtf],
        )
        ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)
        REANNOTATEGTF(GFFCOMPARE.out.annotated_gtf, ref_fai)
        ch_versions = ch_versions.mix(REANNOTATEGTF.out.versions)
        lraa_out_ch = REANNOTATEGTF.out.gtf

        // Classify reannotated isoforms with SQANTI-like categories
        LRAA_SQANTI(
            REANNOTATEGTF.out.gtf,
            ref_gtf,
        )
        ch_versions = ch_versions.mix(LRAA_SQANTI.out.versions)

        // Re-quantify each sample against the reannotated cohort GTF
        quant_input = bam_ch.combine(
            REANNOTATEGTF.out.gtf.map { _meta, gtf -> gtf }
        )
        LRAA_QUANT(quant_input, ref_fasta)
        ch_versions = ch_versions.mix(LRAA_QUANT.out.versions.first())

        // Assemble cohort expression matrix
        quant_files = LRAA_QUANT.out.quant
            .map { _meta, expr -> expr }
            .collect()
            .map { exprs -> [[id: "merge"], exprs] }
        LRAA_QUANTMERGE(quant_files)
        ch_versions = ch_versions.mix(LRAA_QUANTMERGE.out.versions)
    }

    emit:
    versions       = ch_versions
    gtf            = lraa_out_ch // channel: [ val(meta), path(gtf) ]
    quant_matrix   = skip_multisample || sample_count == 1 ? Channel.empty() : LRAA_QUANTMERGE.out.matrix // optional
    sqanti_summary = LRAA_SQANTI.out.summary // channel: [ val(meta), path(tsv) ]
    sqanti_plot    = LRAA_SQANTI.out.plot // channel: [ val(meta), path(pdf) ]
}
