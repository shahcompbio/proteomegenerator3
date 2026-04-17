// SQANTI-like isoform structure classification using LRAA's built-in utility
// Classifies assembled isoforms against a reference annotation

process LRAA_SQANTI {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:0.15.0"

    input:
    tuple val(meta), path(input_gtf)
    path ref_gtf

    output:
    tuple val(meta), path("${prefix}.iso_cats.tsv"), emit: tsv
    tuple val(meta), path("${prefix}.iso_cats.summary_counts.tsv"), emit: summary
    tuple val(meta), path("${prefix}.iso_cats.summary_counts.pdf"), emit: plot
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    SQANTI-like_cats_for_reads_or_isoforms.py \\
        --ref_gtf ${ref_gtf} \\
        --input_gtf ${input_gtf} \\
        --output_prefix ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: "\$(LRAA --version 2>&1 | head -1)"
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.iso_cats.tsv
    touch ${prefix}.iso_cats.summary_counts.tsv
    touch ${prefix}.iso_cats.summary_counts.pdf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: "\$(LRAA --version 2>&1 | head -1)"
    END_VERSIONS
    """
}
