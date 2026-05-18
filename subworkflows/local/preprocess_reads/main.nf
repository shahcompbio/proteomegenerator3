// preprocess reads via filtering and read class creation with bambu
include { SAMTOOLS_INDEX                } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_VIEW ; SAMTOOLS_VIEW as FILTER_READS } from '../../../modules/nf-core/samtools/view/main'
include { SAMTOOLS_FAIDX                } from '../../../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_CONVERT              } from '../../../modules/nf-core/samtools/convert/main'
include { BAMBU_READCLASSES             } from '../../../modules/local/bambu/readclasses/main'

workflow PREPROCESS_READS {
    take:
    input_bam_ch // channel: [ val(meta), [ bam ] ]
    filter_reads // boolean; filter reads on mapq and read length
    filter_acc_reads // boolean; filter reads on accessory chromosomes
    long_read_assembler // string: bambu, lraa, or stringtie

    main:

    ch_versions = Channel.empty()
    ch_cram = Channel.empty()

    if (filter_reads || filter_acc_reads) {
        SAMTOOLS_INDEX(input_bam_ch)
        ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions.first())
        bam_bai_ch = input_bam_ch.join(SAMTOOLS_INDEX.out.bai, by: 0)

        SAMTOOLS_FAIDX([[id: 'ref'], params.fasta, []], false)
        ref_fasta_fai = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> [[id: 'ref'], params.fasta, fai] }
    }
    if (filter_reads) {
        ch_filtered = FILTER_READS(
            bam_bai_ch,
            [[], []],
            [],
            "bai",
        )
        ch_bam = ch_filtered.bam
        ch_versions = ch_versions.mix(FILTER_READS.out.versions.first())
    }
    else if (filter_acc_reads) {
        ch_filtered = SAMTOOLS_VIEW(
            bam_bai_ch,
            [[], []],
            [],
            "bai",
        )
        ch_bam = ch_filtered.bam
        ch_versions = ch_versions.mix(SAMTOOLS_VIEW.out.versions.first())
    }
    else {
        ch_bam = input_bam_ch
    }
    // convert to cram if bam was filtered
    if (filter_reads || filter_acc_reads) {
        SAMTOOLS_CONVERT(
            ch_bam.map { meta, bam -> [meta, bam, []] },
            ref_fasta_fai,
        )
        ch_cram = SAMTOOLS_CONVERT.out.cram
    }
    // create read classes with bambu (only when using bambu assembler)
    if (long_read_assembler.split(',').contains('bambu')) {
        BAMBU_READCLASSES(
            ch_bam,
            params.yieldsize,
            params.fasta,
            params.gtf,
        )
        ch_versions = ch_versions.mix(BAMBU_READCLASSES.out.versions)
    }

    emit:
    bam      = ch_bam // channel: [ val(meta), path(bam) ]
    cram     = ch_cram // channel: [ val(meta), path(cram) ]
    reads    = long_read_assembler.split(',').contains('bambu') ? BAMBU_READCLASSES.out.rds : Channel.empty() // channel: [ val(meta), [ rcFile ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
