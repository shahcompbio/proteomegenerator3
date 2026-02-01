// transcript assembly and quant with stringtie (short-read only for now)
include { STRINGTIE_STRINGTIE } from '../../../modules/nf-core/stringtie/stringtie/main'
include { STRINGTIE_MERGE     } from '../../../modules/nf-core/stringtie/merge/main'
include { GFFCOMPARE          } from '../../../modules/nf-core/gffcompare/main'
include { REANNOTATESTRINGTIE } from '../../../modules/local/reannotatestringtie/main'
workflow BAM_ASSEMBLY_STRINGTIE {
    take:
    ch_bam // channel: [ val(meta), [ bam ] ]
    ref_gtf

    main:
    ch_versions = channel.empty()
    STRINGTIE_STRINGTIE(
        ch_bam,
        ref_gtf,
    )
    ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions)
    // STRINGTIE_STRINGTIE.out.transcript_gtf.view()
    merge_ch = STRINGTIE_STRINGTIE.out.transcript_gtf.collect { _meta, gtf -> gtf }
    // merge assembled transcripts across samples
    STRINGTIE_MERGE(
        merge_ch,
        ref_gtf,
    )
    ch_versions = ch_versions.mix(STRINGTIE_MERGE.out.versions)
    // run gffcompare to determine which transcripts are non-canonical
    // STRINGTIE_MERGE.out.gtf.view()
    GFFCOMPARE(
        STRINGTIE_MERGE.out.gtf.map { gtf -> [[id: "merge"], gtf] },
        [[], [], []],
        [[id: "ref"], ref_gtf],
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)
    // reannotate stringtie output
    REANNOTATESTRINGTIE(GFFCOMPARE.out.annotated_gtf)
    ch_versions = ch_versions.mix(REANNOTATESTRINGTIE.out.versions)

    emit:
    versions = ch_versions
}
