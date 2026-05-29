// SQANTI3 Quality Control — classify isoforms against reference annotation
// Produces classification table, corrected GTF, and junctions file

process SQANTI3_QC {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/sqanti3:6.0.1--hdfd78af_0'
        : 'biocontainers/sqanti3:6.0.1--hdfd78af_0'}"

    input:
    tuple val(meta), path(gtf)
    path ref_gtf
    path ref_fasta

    output:
    tuple val(meta), path("${prefix}_classification.txt"), emit: classification
    tuple val(meta), path("${prefix}_corrected.gtf"), emit: corrected_gtf
    tuple val(meta), path("${prefix}_junctions.txt"), emit: junctions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    sqanti3 qc \\
        --isoforms ${gtf} \\
        --refGTF ${ref_gtf} \\
        --refFasta ${ref_fasta} \\
        --skipORF \\
        --force_id_ignore \\
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
    touch ${prefix}_classification.txt
    touch ${prefix}_corrected.gtf
    touch ${prefix}_junctions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: "6.0.1"
    END_VERSIONS
    """
}
