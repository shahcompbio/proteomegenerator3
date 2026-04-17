/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { PREPROCESS_READS                                    } from '../subworkflows/local/preprocess_reads/main'
include { BAM_ASSEMBLY_BAMBU                                  } from '../subworkflows/local/bam_assembly_bambu/main'
include { BAM_ASSEMBLY_LRAA                                   } from '../subworkflows/local/bam_assembly_lraa/main'
include { GFFREAD                                             } from '../modules/nf-core/gffread/main'
include { CAT_CAT                                             } from '../modules/nf-core/cat/cat/main'
include { PREDICT_ORFS                                        } from '../subworkflows/local/predict_orfs/main'
include { FASTA_MERGE_ANNOTATE                                } from '../subworkflows/local/fasta_merge_annotate/main'
include { BAM_ASSEMBLY_STRINGTIE as BAM_ASSEMBLY_STRINGTIE_LR } from '../subworkflows/local/bam_assembly_stringtie/main'
include { BAM_ASSEMBLY_STRINGTIE as BAM_ASSEMBLY_STRINGTIE_SR } from '../subworkflows/local/bam_assembly_stringtie/main'
include { MULTIQC                                             } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap                                    } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                                } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                              } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                              } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getLongReadBams                                     } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getLongReadRcFiles                                  } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getLongReadLraaGtfs                                 } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getShortReadBams                                    } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
include { getFusionTsvs                                       } from '../subworkflows/local/utils_nfcore_proteomegenerator3_pipeline'
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
    ch_long_read_rc = getLongReadRcFiles(ch_samplesheet)
    ch_long_read_lraa_gtfs = getLongReadLraaGtfs(ch_samplesheet)
    ch_short_read_bams = getShortReadBams(ch_samplesheet)
    ch_fusion_tsvs = getFusionTsvs(ch_samplesheet)

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
    // perform assembly & quantification with bambu
    // make an NDR channel
    if (params.recommended_NDR && params.NDR != null) {
        ch_NDR = channel.of("DEFAULT", params.NDR)
    }
    else if (params.recommended_NDR) {
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
    if (params.long_read_assembler == 'bambu') {
        BAM_ASSEMBLY_BAMBU(
            rc_ch,
            params.skip_multisample,
            sample_count,
            ch_NDR,
            ref_gtf_ch,
            bam_ch,
        )
        ch_versions = ch_versions.mix(BAM_ASSEMBLY_BAMBU.out.versions)
        assembly_ch = BAM_ASSEMBLY_BAMBU.out.gtf.map { meta, gtf -> [meta + [tool: 'bambu'], gtf] }
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
        assembly_ch = BAM_ASSEMBLY_LRAA.out.gtf.map { meta, gtf -> [meta + [tool: 'lraa'], gtf] }
    }
    else if (params.long_read_assembler == 'stringtie') {
        BAM_ASSEMBLY_STRINGTIE_LR(
            bam_ch,
            params.gtf,
            params.skip_multisample,
            sample_count,
        )
        ch_versions = ch_versions.mix(BAM_ASSEMBLY_STRINGTIE_LR.out.versions)
        assembly_ch = BAM_ASSEMBLY_STRINGTIE_LR.out.gtf.map { meta, gtf -> [meta + [tool: 'stringtie_lr'], gtf] }
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
        )
        ch_versions = ch_versions.mix(BAM_ASSEMBLY_STRINGTIE_SR.out.versions)
        // combine LR and SR assemblies
        stringtie_ch = BAM_ASSEMBLY_STRINGTIE_SR.out.gtf.map { meta, gtf -> [meta + [tool: 'stringtie'], gtf] }
        assembly_ch = assembly_ch.mix(stringtie_ch)
    }
    // extract cDNA
    GFFREAD(assembly_ch, params.fasta)
    ch_versions = ch_versions.mix(GFFREAD.out.versions)
    // predict ORFs with transdecoder and output fasta for msfragger
    PREDICT_ORFS(GFFREAD.out.gffread_fasta, params.uniprot_proteome)
    ch_versions = ch_versions.mix(PREDICT_ORFS.out.versions)
    // make uniprot-style fasta for msfragger and create index tables
    ch_orfs = PREDICT_ORFS.out.ORFs
        .join(assembly_ch, by: 0)
        .combine(PREDICT_ORFS.out.swissprot.map { _meta, fasta -> fasta })
    // ch_orfs.view { v -> "ch_orfs: ${v}" }
    FASTA_MERGE_ANNOTATE(
        ch_orfs,
        params.input,
        params.skip_multisample,
        PREDICT_ORFS.out.swissprot,
        ch_fusion_tsvs,
        params.fusions,
        params.short_reads,
    )
    ch_versions = ch_versions.mix(FASTA_MERGE_ANNOTATE.out.versions)
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
