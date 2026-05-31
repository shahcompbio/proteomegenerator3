// SQANTI3 Filter — remove artifact isoforms using ML or rules-based approach
// Takes QC output and produces a curated transcriptome

process SQANTI3_FILTER {
    tag "${meta.id}"
    label 'process_medium'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://anaconesalab/sqanti3:v6.0.1'
        : 'docker.io/anaconesalab/sqanti3:v6.0.1'}"

    input:
    tuple val(meta), path(classification), path(corrected_gtf)

    output:
    tuple val(meta), path("*result_classification.txt"), emit: classification
    tuple val(meta), path("*inclusion-list.txt"), emit: inclusion_list
    tuple val(meta), path("*.filtered.gtf"), emit: filtered_gtf
    tuple val(meta), path("*randomforest.RData"), emit: random_forest, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def filter_type = task.ext.filter_type ?: 'ml'
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    export TZ=UTC

    sqanti3_filter.py \\
        ${filter_type} \\
        --sqanti_class ${classification} \\
        --filter_gtf ${corrected_gtf} \\
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
    touch ${prefix}_MLresult_classification.txt
    touch ${prefix}_inclusion-list.txt
    touch ${prefix}.filtered.gtf
    touch ${prefix}_randomforest.RData

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: "6.0.1"
    END_VERSIONS
    """
}
