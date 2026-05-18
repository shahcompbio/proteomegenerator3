/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_READS                                    } from '../subworkflows/local/preprocess_reads/main'
include { BAM_QC                                              } from '../subworkflows/local/bam_qc/main'
include { BAM_ASSEMBLY_BAMBU                                  } from '../subworkflows/local/bam_assembly_bambu/main'
include { BAM_ASSEMBLY_LRAA                                   } from '../subworkflows/local/bam_assembly_lraa/main'
include { GFFREAD                                             } from '../modules/nf-core/gffread/main'
include { SAMTOOLS_FAIDX                                      } from '../modules/nf-core/samtools/faidx/main'
include { SAMTOOLS_CONVERT as CRAM_TO_BAM                     } from '../modules/nf-core/samtools/convert/main'
include { PREDICT_ORFS                                        } from '../subworkflows/local/predict_orfs/main'
include { FASTA_MERGE_ANNOTATE                                } from '../subworkflows/local/fasta_merge_annotate/main'
include { GTF_MERGE_ANNOTATE                                  } from '../subworkflows/local/gtf_merge_annotate/main'
include { GTF_MERGE_ANNOTATE as GTF_MERGE_ORFS_ONLY           } from '../subworkflows/local/gtf_merge_annotate/main'
include { BAM_ASSEMBLY_STRINGTIE as BAM_ASSEMBLY_STRINGTIE_LR } from '../subworkflows/local/bam_assembly_stringtie/main'
include { BAM_ASSEMBLY_STRINGTIE as BAM_ASSEMBLY_STRINGTIE_SR } from '../subworkflows/local/bam_assembly_stringtie/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getLongReadBams                                     } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getLongReadRcFiles                                  } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getGtfs                                             } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getShortReadBams                                    } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getFusionTsvs                                       } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getLongReadCrams                                    } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PROTEOMEGENERATOR3 {
    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    // begin workflow
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()
    //
    // Extract typed channels from long-format samplesheet
    //
    ch_long_read_bams = getLongReadBams(ch_samplesheet)
    ch_long_read_crams = getLongReadCrams(ch_samplesheet)
    ch_long_read_rc = getLongReadRcFiles(ch_samplesheet)
    ch_gtfs = getGtfs(ch_samplesheet)
    ch_short_read_bams = getShortReadBams(ch_samplesheet)
    ch_fusion_tsvs = getFusionTsvs(ch_samplesheet)

    //
    // Index reference genome (needed for reannotation and CRAM conversion)
    //
    SAMTOOLS_FAIDX([[id: 'ref'], params.fasta, []], false)
    ref_fai = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> fai }
    ref_fasta_fai = SAMTOOLS_FAIDX.out.fai.map { _meta, fai -> [[id: 'ref'], file(params.fasta), fai] }

    if (!params.orfs_only) {
        //
        // Convert CRAM inputs to BAM (if any)
        //
        CRAM_TO_BAM(
            ch_long_read_crams.map { meta, cram -> [meta, cram, []] },
            ref_fasta_fai,
        )
        ch_long_read_bams = ch_long_read_bams.mix(CRAM_TO_BAM.out.bam)

        //
        // process long-read rnaseq data
        //
        if (!params.skip_preprocessing) {
            PREPROCESS_READS(ch_long_read_bams, params.filter_reads, params.filter_acc_reads, params.long_read_assembler)
            rc_ch = PREPROCESS_READS.out.reads
            bam_ch = PREPROCESS_READS.out.bam
            ch_versions = ch_versions.mix(PREPROCESS_READS.out.versions)
        }
        else {
            // Use provided rc_files when skipping preprocessing
            rc_ch = ch_long_read_rc
            bam_ch = ch_long_read_bams
        }
        // perform qc on filtered bams
        if (!params.skip_qc) {
            BAM_QC(rc_ch, bam_ch)
            ch_versions = ch_versions.mix(BAM_QC.out.versions)
            ch_multiqc_files = ch_multiqc_files.mix(BAM_QC.out.multiqc)
        }
    }
    // end if (!params.orfs_only)

    if (!params.qc_only) {
        if (params.orfs_only) {
            //
            // ORFs-only mode: skip assembly, use pre-computed GTFs
            //
            sample_count = countSamples(params.input)

            if (!params.skip_multisample && sample_count > 1) {
                // Merge all GTFs into a cohort-level GTF via GTF_MERGE_ANNOTATE
                merge_input_ch = ch_gtfs.map { meta, gtf ->
                    [meta + [tool: 'user_gtf'], gtf]
                }
                GTF_MERGE_ORFS_ONLY(merge_input_ch, params.gtf, ref_fai)
                ch_versions = ch_versions.mix(GTF_MERGE_ORFS_ONLY.out.versions)
                downstream_gtf_ch = GTF_MERGE_ORFS_ONLY.out.gtf
            }
            else {
                // Single GTF or skip_multisample: pass GTFs directly (no reannotation)
                downstream_gtf_ch = ch_gtfs
            }
        }
        else {
            //
            // Standard assembly mode
            //
            // make an NDR channel (single value only)
            if (params.recommended_NDR) {
                ch_NDR = channel.of("DEFAULT")
            }
            else {
                ch_NDR = channel.of(params.NDR)
            }
            ref_gtf_ch = channel.of(params.gtf)
            // run sample assembly & quant with read classes
            // count samples to make sure multisample isn't run on single samples
            sample_count = countSamples(params.input)

            //
            // Long-read assembly: select assembler
            //
            assembly_ch = Channel.empty()

            if (params.long_read_assembler.split(',').contains('bambu')) {
                BAM_ASSEMBLY_BAMBU(
                    rc_ch,
                    params.skip_multisample,
                    sample_count,
                    ch_NDR,
                    ref_gtf_ch,
                    bam_ch,
                )
                ch_versions = ch_versions.mix(BAM_ASSEMBLY_BAMBU.out.versions)
                assembly_ch = assembly_ch.mix(BAM_ASSEMBLY_BAMBU.out.gtf.map { meta, gtf -> [meta + [tool: 'bambu'], gtf] })
            }
            if (params.long_read_assembler.split(',').contains('lraa')) {
                BAM_ASSEMBLY_LRAA(
                    bam_ch,
                    params.skip_multisample,
                    sample_count,
                    params.gtf,
                    params.fasta,
                    ref_fai,
                )
                ch_versions = ch_versions.mix(BAM_ASSEMBLY_LRAA.out.versions)
                assembly_ch = assembly_ch.mix(BAM_ASSEMBLY_LRAA.out.gtf.map { meta, gtf -> [meta + [tool: 'lraa'], gtf] })
            }
            if (params.long_read_assembler.split(',').contains('stringtie')) {
                BAM_ASSEMBLY_STRINGTIE_LR(
                    bam_ch,
                    params.gtf,
                    params.skip_multisample,
                    sample_count,
                    ref_fai,
                )
                ch_versions = ch_versions.mix(BAM_ASSEMBLY_STRINGTIE_LR.out.versions)
                assembly_ch = assembly_ch.mix(BAM_ASSEMBLY_STRINGTIE_LR.out.gtf.map { meta, gtf -> [meta + [tool: 'stringtie_lr'], gtf] })
            }

            //
            // process short-read rnaseq data (if provided)
            //
            if (params.short_reads) {
                BAM_ASSEMBLY_STRINGTIE_SR(
                    ch_short_read_bams,
                    params.gtf,
                    params.skip_multisample,
                    sample_count,
                    ref_fai,
                )
                ch_versions = ch_versions.mix(BAM_ASSEMBLY_STRINGTIE_SR.out.versions)
                // combine LR and SR assemblies
                stringtie_ch = BAM_ASSEMBLY_STRINGTIE_SR.out.gtf.map { meta, gtf -> [meta + [tool: 'stringtie_sr'], gtf] }
                assembly_ch = assembly_ch.mix(stringtie_ch)
            }
            //
            // Merge assembler GTFs into consensus assembly (only when multiple assemblers)
            //
            def assembler_count = params.long_read_assembler.split(',').size() + (params.short_reads ? 1 : 0)
            if (assembler_count > 1) {
                GTF_MERGE_ANNOTATE(assembly_ch, params.gtf, ref_fai)
                ch_versions = ch_versions.mix(GTF_MERGE_ANNOTATE.out.versions)
                downstream_gtf_ch = GTF_MERGE_ANNOTATE.out.gtf
            }
            else {
                // Single assembler: use its GTF directly (already annotated)
                downstream_gtf_ch = assembly_ch
            }
        }
        // end assembly mode selection

        //
        // Downstream: single path
        //
        // Extract cDNA
        GFFREAD(downstream_gtf_ch, params.fasta)
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
        // Predict ORFs with transdecoder
        PREDICT_ORFS(GFFREAD.out.gffread_fasta, params.uniprot_proteome)
        ch_versions = ch_versions.mix(PREDICT_ORFS.out.versions)
        // Make uniprot-style fasta for msfragger and create index tables
        ch_orfs = PREDICT_ORFS.out.ORFs
            .join(downstream_gtf_ch, by: 0)
            .combine(PREDICT_ORFS.out.swissprot.map { _meta, fasta -> fasta })
        FASTA_MERGE_ANNOTATE(
            ch_orfs,
            params.input,
            params.skip_multisample,
            PREDICT_ORFS.out.swissprot,
            ch_fusion_tsvs,
            params.fusions,
        )
        ch_versions = ch_versions.mix(FASTA_MERGE_ANNOTATE.out.versions)
    }
    // end if (!params.qc_only)
    // collect versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'proteomegenerator3_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()

    summary_params = paramsSummaryMap(
        workflow,
        parameters_schema: "nextflow_schema.json"
    )
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml')
    )
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def countSamples(input) {
    def lines = file(input).readLines()
    def header = lines[0].split(',')
    def sampleIdx = header.findIndexOf { it == 'sample_id' }

    // Get unique sample names (excluding header) for long-format samplesheet
    def samples = lines[1..-1]
        .collect { it.split(',')[sampleIdx] }
        .unique()

    def sample_count = samples.size()
    if (sample_count == 1) {
        println("1 sample detected; switching to single sample mode")
    }
    return sample_count
}
