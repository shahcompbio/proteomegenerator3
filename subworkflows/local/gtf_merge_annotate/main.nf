// Merge per-tool assembler GTFs into a single union assembly,
// reannotate with reference IDs, track provenance, and classify isoforms.

include { LRAA_MERGE                          } from '../../../modules/local/lraa/merge/main'
include { GFFCOMPARE ; GFFCOMPARE as GFFCOMPARE_PROVENANCE } from '../../../modules/nf-core/gffcompare/main'
include { REANNOTATEGTF                       } from '../../../modules/local/reannotategtf/main'
include { LRAA_SQANTI                         } from '../../../modules/local/lraa/sqanti/main'

workflow GTF_MERGE_ANNOTATE {
    take:
    assembly_ch // channel: [ val(meta), path(gtf) ] — per-tool GTFs with meta.tool set
    ref_gtf // val: path to reference GTF
    ref_fai // channel: path to reference FASTA index
    ref_fasta // val: path to reference genome FASTA

    main:
    ch_versions = Channel.empty()
    // Collect all per-tool GTFs and merge with LRAA splice graph reconstruction
    merge_gtfs = assembly_ch
        .map { meta, gtf ->
            [[id: meta.subject_id], gtf]
        }
        .groupTuple()
    LRAA_MERGE(merge_gtfs, ref_fasta)
    ch_versions = ch_versions.mix(LRAA_MERGE.out.versions)

    // Annotate merged GTF against reference annotation
    GFFCOMPARE(
        LRAA_MERGE.out.gtf,
        [[], [], []],
        [[id: "ref"], ref_gtf],
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    // Reannotate with reference IDs (ENST/ENSG for canonical, novel prefixes for new)
    REANNOTATEGTF(GFFCOMPARE.out.annotated_gtf, ref_fai)
    ch_versions = ch_versions.mix(REANNOTATEGTF.out.versions)

    // Provenance tracking: which tool(s) contributed each merged transcript
    GFFCOMPARE_PROVENANCE(
        assembly_ch.map { _meta, gtf -> gtf }.collect().map { gtfs -> [[id: "provenance"], gtfs] },
        [[], [], []],
        REANNOTATEGTF.out.gtf,
    )
    ch_versions = ch_versions.mix(GFFCOMPARE_PROVENANCE.out.versions)

    // SQANTI classification of union isoforms
    LRAA_SQANTI(REANNOTATEGTF.out.gtf, ref_gtf)
    ch_versions = ch_versions.mix(LRAA_SQANTI.out.versions)

    emit:
    gtf      = REANNOTATEGTF.out.gtf.map { meta, gtf -> [meta + [tool: 'union'], gtf] } // channel: [ val(meta), path(gtf) ] — reannotated union GTF
    mapping  = REANNOTATEGTF.out.mapping // channel: [ val(meta), path(tsv) ] — ID mapping table
    tracking = GFFCOMPARE_PROVENANCE.out.tracking // channel: [ val(meta), path(tracking) ] — provenance
    sqanti   = LRAA_SQANTI.out.tsv // channel: [ val(meta), path(tsv) ] — isoform classification
    versions = ch_versions // channel: [ versions.yml ]
}
