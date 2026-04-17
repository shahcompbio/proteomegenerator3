// Per-sample re-quantification against a merged cohort GTF
process LRAA_QUANT {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:v0.15.0"

    input:
    tuple val(meta), path(bam), path(cohort_gtf, stageAs: 'cohort.gtf')
    path ref_genome

    output:
    tuple val(meta), path("${prefix}.quant.expr"), emit: quant
    tuple val(meta), path("${prefix}.quant.tracking"), emit: tracking
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    LRAA \\
        --quant_only \\
        --genome ${ref_genome} \\
        --gtf    cohort.gtf \\
        --bam    ${bam} \\
        --output_prefix ${prefix} \\
        --num_parallel_contigs ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.quant.expr
    touch ${prefix}.quant.tracking

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """
}
