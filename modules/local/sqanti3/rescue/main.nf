// SQANTI3 Rescue — recover reference transcripts for discarded artifacts
// Expands filtered transcriptome with matching reference models

process SQANTI3_RESCUE {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sqanti3:6.0.1--hdfd78af_0'
        : 'biocontainers/sqanti3:6.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(classification), path(filtered_gtf)
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
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sqanti3 rescue \\
        --isoforms ${classification} \\
        --gtf ${filtered_gtf} \\
        --refGTF ${ref_gtf} \\
        --refFasta ${ref_fasta} \\
        -o ${prefix} \\
        -d . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3 --version 2>&1 | sed -n 's/.*version //p' || echo "6.0.1")
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
