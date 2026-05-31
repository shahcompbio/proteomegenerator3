// SQANTI3 Rescue — recover reference transcripts for discarded artifacts
// Expands filtered transcriptome with matching reference models

process SQANTI3_RESCUE {
    tag "${meta.id}"
    label 'process_medium'

    container "quay.io/shahlab_singularity/sqanti3:6.0.1"

    input:
    tuple val(meta), path(classification), path(filtered_gtf), path(corrected_fasta), path(random_forest)
    path ref_gtf
    path ref_fasta

    output:
    tuple val(meta), path("${prefix}*rescued*.gtf"), emit: rescued_gtf
    tuple val(meta), path("${prefix}*rescued*classification*"), emit: classification
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def rescue_type = task.ext.rescue_type ?: 'rules'
    def rf_arg = (rescue_type == 'ml' && random_forest) ? "-r ${random_forest}" : ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TZ=UTC

    sqanti3_rescue.py \\
        --filter_class ${classification} \\
        --filtered_isoforms_gtf ${filtered_gtf} \\
        -rg ${ref_gtf} \\
        -rf ${ref_fasta} \\
        -s ${rescue_type} \\
        ${rf_arg} \\
        --corrected_isoforms_fasta ${corrected_fasta} \\
        -o ${prefix} \\
        -d . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py -v 2>&1 | grep 'SQANTI3' | sed 's/.*SQANTI3 //' || echo "6.0.1")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_rescued.gtf
    touch ${prefix}_rescued_classification.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: "6.0.1"
    END_VERSIONS
    """
}
