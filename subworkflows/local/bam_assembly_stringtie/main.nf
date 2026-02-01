// transcript assembly and quant with stringtie (short-read only for now)
include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'
include { GFFCOMPARE          } from '../../../modules/nf-core/gffcompare/main'
include { REANNOTATESTRINGTIE } from '../../../modules/local/reannotatestringtie/main'
workflow BAM_ASSEMBLY_STRINGTIE {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ref_gtf
    skip_multisample // val
    sample_count // val

    main:
    ch_versions = channel.empty()
    stringtie_out_ch = channel.empty()
    STRINGTIE_STRINGTIE(
        ch_bam,
        ref_gtf,
    )
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)
    if (skip_multisample || sample_count == 1) {
        stringtie_out_ch = stringtie_out_ch.mix(STRINGTIE_STRINGTIE.out.transcript_gtf)
    }
    else {

        // STRINGTIE_STRINGTIE.out.transcript_gtf.view()
        merge_ch = STRINGTIE_STRINGTIE.out.transcript_gtf.collect { _meta, gtf -> gtf }
        // merge assembled transcripts across samples
        STRINGTIE_MERGE(
            merge_ch,
            [],
        )
        ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
        stringtie_out_ch = stringtie_out_ch.mix(
            STRINGTIE_MERGE.out.gtf.map { gtf -> [[id: "merge"], gtf] }
        )
    }
    // run gffcompare to determine which transcripts are non-canonical
    GFFCOMPARE(
        stringtie_out_ch,
        [[], [], []],
        [[id: "ref"], ref_gtf],
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)
    // reannotate stringtie output
    REANNOTATESTRINGTIE(GFFCOMPARE.out.annotated_gtf)
    ch_versions = ch_versions.mix(REANNOTATESTRINGTIE.out.versions)

    emit:
    versions = ch_versions
    gtf      = REANNOTATESTRINGTIE.out.gtf
}
