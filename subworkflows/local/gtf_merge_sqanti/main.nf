// Merge per-tool assembler GTFs into a single union assembly,
// reannotate with reference IDs, track provenance, and curate with SQANTI3.
//
// Supports both multi-assembler (LRAA_MERGE → annotate → SQANTI3) and
// single-assembler (skip merge, annotate → SQANTI3) modes.

include { LRAA_MERGE                          } from '../../../modules/local/lraa/merge/main'
include { GFFCOMPARE ; GFFCOMPARE as GFFCOMPARE_PROVENANCE } from '../../../modules/nf-core/gffcompare/main'
include { REANNOTATEGTF                       } from '../../../modules/local/reannotategtf/main'
include { SQANTI3_QC                          } from '../../../modules/local/sqanti3/qc/main'
include { SQANTI3_FILTER                      } from '../../../modules/local/sqanti3/filter/main'
include { SQANTI3_RESCUE                      } from '../../../modules/local/sqanti3/rescue/main'

workflow GTF_MERGE_SQANTI {
    take:
    assembly_ch       // channel: [ val(meta), path(gtf) ] — per-tool GTFs with meta.tool set
    ref_gtf           // val: path to reference GTF
    ref_fai           // channel: path to reference FASTA index
    ref_fasta         // val: path to reference genome FASTA
    skip_sqanti3      // val: boolean — skip SQANTI3 QC/filter/rescue
    skip_union_assembly // val: boolean — skip LRAA_MERGE across assemblers

    main:
    ch_versions = Channel.empty()

    //
    // Step 1: Conditional merge (only when multiple assemblers and not skipped)
    //
    if (!skip_union_assembly) {
        // Group GTFs by subject_id to detect multi-assembler vs single-assembler
        grouped_gtfs = assembly_ch
            .map { meta, gtf ->
                [[id: meta.subject_id], gtf]
            }
            .groupTuple()

        // Split into groups that need merging (>1 GTF) and those that don't (1 GTF)
        grouped_gtfs
            .branch {
                merge: it[1].size() > 1
                single: true
            }
            .set { branch_gtfs }

        // Merge multi-assembler GTFs with LRAA splice graph reconstruction
        LRAA_MERGE(branch_gtfs.merge, ref_fasta)
        ch_versions = ch_versions.mix(LRAA_MERGE.out.versions)

        // Combine: merged GTFs + single-assembler GTFs (passed through)
        merged_or_single_gtf = LRAA_MERGE.out.gtf
            .mix(branch_gtfs.single.map { meta, gtfs -> [meta, gtfs[0]] })
    }
    else {
        // Skip union assembly: pass each GTF through individually
        merged_or_single_gtf = assembly_ch
            .map { meta, gtf -> [[id: meta.subject_id], gtf] }
    }

    //
    // Step 2: Annotate against reference
    //
    GFFCOMPARE(
        merged_or_single_gtf,
        [[], [], []],
        [[id: "ref"], ref_gtf],
    )
    ch_versions = ch_versions.mix(GFFCOMPARE.out.versions)

    // Reannotate with reference IDs (ENST/ENSG for canonical, novel prefixes for new)
    REANNOTATEGTF(GFFCOMPARE.out.annotated_gtf, ref_fai)
    ch_versions = ch_versions.mix(REANNOTATEGTF.out.versions)

    //
    // Step 3: Provenance tracking (which tool(s) contributed each transcript)
    //
    GFFCOMPARE_PROVENANCE(
        assembly_ch.map { _meta, gtf -> gtf }.collect().map { gtfs -> [[id: "provenance"], gtfs] },
        [[], [], []],
        REANNOTATEGTF.out.gtf,
    )
    ch_versions = ch_versions.mix(GFFCOMPARE_PROVENANCE.out.versions)

    //
    // Step 4: SQANTI3 QC → Filter → Rescue (unless skipped)
    //
    if (!skip_sqanti3) {
        // QC: classify isoforms against reference (skipORF, force_id_ignore)
        SQANTI3_QC(REANNOTATEGTF.out.gtf, ref_gtf, ref_fasta)
        ch_versions = ch_versions.mix(SQANTI3_QC.out.versions)

        // Filter: ML-based artifact removal
        sqanti3_filter_input = SQANTI3_QC.out.classification
            .join(SQANTI3_QC.out.corrected_gtf, by: 0)
        SQANTI3_FILTER(sqanti3_filter_input)
        ch_versions = ch_versions.mix(SQANTI3_FILTER.out.versions)

        // Rescue: recover reference transcripts for discarded artifacts
        sqanti3_rescue_input = SQANTI3_FILTER.out.classification
            .join(SQANTI3_FILTER.out.filtered_gtf, by: 0)
            .join(SQANTI3_QC.out.corrected_fasta, by: 0)
            .join(SQANTI3_FILTER.out.random_forest, by: 0, remainder: true)
            .map { meta, classif, gtf, fasta, rf -> [meta, classif, gtf, fasta, rf ?: []] }
        SQANTI3_RESCUE(sqanti3_rescue_input, ref_gtf, ref_fasta)
        ch_versions = ch_versions.mix(SQANTI3_RESCUE.out.versions)

        // Final curated GTF comes from rescue
        curated_gtf = SQANTI3_RESCUE.out.rescued_gtf
            .map { meta, gtf -> [meta + [tool: 'union'], gtf] }
        sqanti_classification = SQANTI3_RESCUE.out.classification
    }
    else {
        // Skip SQANTI3: use reannotated GTF directly
        curated_gtf = REANNOTATEGTF.out.gtf
            .map { meta, gtf -> [meta + [tool: 'union'], gtf] }
        sqanti_classification = Channel.empty()
    }

    emit:
    gtf            = curated_gtf              // channel: [ val(meta), path(gtf) ] — curated GTF (SQANTI3 rescued or reannotated)
    mapping        = REANNOTATEGTF.out.mapping // channel: [ val(meta), path(tsv) ] — ID mapping table
    tracking       = GFFCOMPARE_PROVENANCE.out.tracking // channel: [ val(meta), path(tracking) ] — provenance
    classification = sqanti_classification    // channel: [ val(meta), path(txt) ] — SQANTI3 classification (empty if skipped)
    versions       = ch_versions              // channel: [ versions.yml ]
}
