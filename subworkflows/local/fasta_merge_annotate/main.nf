// make uniprot-style fasta for msfragger and create index tables
include { TRANSDECODER2FASTA         } from '../../../modules/local/transdecoder2fasta/main.nf'
include { MERGEFUSIONS               } from '../../../modules/local/mergefusions/main.nf'
include { CAT_CAT as CAT_CAT_SAMPLES } from '../../../modules/nf-core/cat/cat/main.nf'
include { CAT_CAT                    } from '../../../modules/nf-core/cat/cat/main.nf'
include { SEQKIT_RMDUP               } from '../../../modules/local/seqkit_rmdup/main.nf'
include { SEQKIT_STATS               } from '../../../modules/nf-core/seqkit/stats/main'
include { FUSIONFASTA                } from '../../../modules/local/fusionfasta/main.nf'

workflow FASTA_MERGE_ANNOTATE {
    take:
    ch_orfs // channel: [ val(meta), [ transdecoder_peps, gtf, swissprot_fasta ] ]
    samplesheet // samplesheet path for multi-sample fusion merging
    skip_multisample // boolean determining multi or single sample mode
    swissprot_fasta // swissprot fasta
    ch_fusion_tsvs // channel [ val(meta), path(fusion_tsv) ] - pre-filtered fusion entries
    run_fusions // flag to include fusions
    short_reads // flag to include short read data

    main:

    ch_versions = Channel.empty()

    // uniprot-style fasta and index table for transdecoder orfs
    TRANSDECODER2FASTA(ch_orfs)
    // TRANSDECODER2FASTA.out.fasta.view { v -> "Transdecoder fasta channel: ${v}" }
    if (short_reads) {
        TRANSDECODER2FASTA.out.fasta
            .branch { meta, _fasta ->
                bambu: meta.tool == 'bambu'
                return tuple(meta.subject_id ?: meta.id, meta, _fasta)
                stringtie: meta.tool == 'stringtie'
                return tuple(meta.subject_id ?: meta.id, meta, _fasta)
            }
            .set { branch_ch }
        fasta_ch = branch_ch.bambu
            .combine(branch_ch.stringtie, by: 0)
            .map { id, _meta1, bambu_fasta, _meta2, stringtie_fasta ->
                [[id: id], [bambu_fasta, stringtie_fasta]]
            }
        CAT_CAT_SAMPLES(fasta_ch)
        ch_versions = ch_versions.mix(CAT_CAT_SAMPLES.out.versions.first())
        fasta_ch = CAT_CAT_SAMPLES.out.file_out
    }
    else {
        fasta_ch = TRANSDECODER2FASTA.out.fasta
    }
    ch_versions = ch_versions.mix(TRANSDECODER2FASTA.out.versions.first())
    // multisample workflow with merging
    if (!skip_multisample & run_fusions) {
        // merge fusions if multisample mode has been enabled
        MERGEFUSIONS(samplesheet)
        ch_versions = ch_versions.mix(MERGEFUSIONS.out.versions)
        // concatenate fusions, non-canonical proteins, and swissprot
        cat_ch = fasta_ch
            .combine(MERGEFUSIONS.out.fasta)
            .combine(swissprot_fasta)
            .map { meta1, novel_proteins, fusions, _meta2, sp_fasta ->
                [meta1, [sp_fasta, novel_proteins, fusions]]
            }
    }
    else if (skip_multisample & run_fusions) {
        // ch_fusion_tsvs is already in the format [meta, fusion_tsv]
        FUSIONFASTA(ch_fusion_tsvs)
        // concatenate fusions, non-canonical proteins, and swissprot
        // this time one item per sample
        // TRANSDECODER2FASTA.out.fasta.view()
        cat_ch = fasta_ch
            .map { meta, fasta ->
                [[id: meta.id], meta, fasta]
            }
            .combine(
                FUSIONFASTA.out.fasta.map { meta, fasta -> [[id: meta.subject_id ?: meta.id], fasta] },
                by: 0
            )
            .map { _meta1, meta2, novel_proteins, fusions ->
                tuple(meta2, novel_proteins, fusions)
            }
            .combine(swissprot_fasta)
            .map { meta1, novel_proteins, fusions, _meta2, sp_fasta ->
                [meta1, [sp_fasta, novel_proteins, fusions]]
            }
    }
    else {
        // concatenate transdecoder orfs and swissprot only
        cat_ch = fasta_ch
            .combine(swissprot_fasta)
            .map { meta1, novel_proteins, _meta2, sp_fasta ->
                [meta1, [sp_fasta, novel_proteins]]
            }
    }
    // concat fasta files
    // cat_ch.view()
    CAT_CAT(cat_ch)
    ch_versions = ch_versions.mix(CAT_CAT.out.versions.first())
    // remove duplicates from fasta
    SEQKIT_RMDUP(CAT_CAT.out.file_out)
    ch_versions = ch_versions.mix(SEQKIT_RMDUP.out.versions.first())
    // compute some basic stats on the final proteome
    SEQKIT_STATS(SEQKIT_RMDUP.out.fastx)
    ch_versions = ch_versions.mix(SEQKIT_STATS.out.versions.first())

    emit:
    predicted_orfs = TRANSDECODER2FASTA.out.fasta // channel: [ val(meta), [ fasta ] ]
    orf_tsv        = TRANSDECODER2FASTA.out.tsv // channel: [ val(meta), [ tsv ] ]
    proteome_fasta = SEQKIT_RMDUP.out.fastx // channel: [ val(meta), [fasta]]
    versions       = ch_versions // channel: [ versions.yml ]
}
