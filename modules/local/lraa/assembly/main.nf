// Per-sample guided transcript assembly with LRAA
process LRAA_ASSEMBLY {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:0.15.0"

    input:
    tuple val(meta), path(bam), path(ref_gtf, stageAs: 'reference.gtf')
    path ref_genome

    output:
    tuple val(meta), path("${prefix}.gtf"), emit: gtf
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
        --genome ${ref_genome} \\
        --gtf    reference.gtf \\
        --bam    ${bam} \\
        --output_prefix ${prefix} \\
        --num_parallel_contigs ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: "\$(LRAA --version 2>&1 | head -1)"
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.gtf
    touch ${prefix}.quant.expr
    touch ${prefix}.quant.tracking

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: "\$(LRAA --version 2>&1 | head -1)"
    END_VERSIONS
    """
}
