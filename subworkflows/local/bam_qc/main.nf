// perform QC on filtered bams
include { SAMTOOLS_SORT               } from '../../../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX              } from '../../../modules/nf-core/samtools/index/main'
include { SAMTOOLS_STATS              } from '../../../modules/nf-core/samtools/stats/main'
include { NANOPLOT                    } from '../../../modules/nf-core/nanoplot/main'
include { RSEQC_BAMSTAT               } from '../../../modules/nf-core/rseqc/bamstat/main'
include { UCSC_GTFTOGENEPRED          } from '../../../modules/nf-core/ucsc/gtftogenepred/main'
include { PICARD_COLLECTRNASEQMETRICS } from '../../../modules/nf-core/picard/collectrnaseqmetrics/main'

workflow BAM_QC {
    take:
    ch_rc // channel: [ val(meta), path(rc_file) ] — unused, kept for interface compatibility
    ch_bam // channel: [ val(meta), path(bam) ]

    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Sort and index BAMs
    SAMTOOLS_SORT(ch_bam, [[], [], []], "")
    SAMTOOLS_INDEX(SAMTOOLS_SORT.out.bam)

    // Prepare BAM + BAI channel for tools that need it
    ch_bam_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai, by: 0)

    // SAMTOOLS_STATS: BAM statistics
    SAMTOOLS_STATS(ch_bam_bai, [[], [], []])
    ch_multiqc_files = ch_multiqc_files.mix(SAMTOOLS_STATS.out.stats.collect { it[1] })

    // NANOPLOT: long-read QC plots
    NANOPLOT(SAMTOOLS_SORT.out.bam)
    ch_versions = ch_versions.mix(NANOPLOT.out.versions.first())
    ch_multiqc_files = ch_multiqc_files.mix(NANOPLOT.out.txt.collect { it[1] })

    // RSEQC_BAMSTAT: BAM statistics
    RSEQC_BAMSTAT(ch_bam_bai)
    ch_multiqc_files = ch_multiqc_files.mix(RSEQC_BAMSTAT.out.txt.collect { it[1] })

    // UCSC_GTFTOGENEPRED: convert GTF to refFlat for Picard
    UCSC_GTFTOGENEPRED([[id: 'genome'], params.gtf])
    ref_flat = UCSC_GTFTOGENEPRED.out.refflat.map { _meta, refflat -> refflat }

    // PICARD_COLLECTRNASEQMETRICS: RNA-seq metrics
    PICARD_COLLECTRNASEQMETRICS(SAMTOOLS_SORT.out.bam, ref_flat, params.fasta, [])
    ch_multiqc_files = ch_multiqc_files.mix(PICARD_COLLECTRNASEQMETRICS.out.metrics.collect { it[1] })

    emit:
    bam      = SAMTOOLS_SORT.out.bam // channel: [ val(meta), path(bam) ]
    bai      = SAMTOOLS_INDEX.out.bai // channel: [ val(meta), path(bai) ]
    multiqc  = ch_multiqc_files // channel: collected QC files for MultiQC
    versions = ch_versions // channel: [ path(versions.yml) ]
}
