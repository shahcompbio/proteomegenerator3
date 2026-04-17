// Merge per-sample LRAA GTFs into a cohort-level GTF
process LRAA_MERGE {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:0.15.0"

    input:
    tuple val(meta), path(gtfs, stageAs: 'sample_gtfs/*')
    path ref_genome

    output:
    tuple val(meta), path("*.LRAA.merged.gtf"), emit: gtf
    tuple val(meta), path("*.LRAA.merged.gtf.tracking.tsv"), emit: tracking
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_LRAA_GTFs.py \\
        --genome ${ref_genome} \\
        --gtf sample_gtfs/*.gtf \\
        --output_gtf ${prefix}.LRAA.merged.gtf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.LRAA.merged.gtf"
    touch "${prefix}.LRAA.merged.gtf.tracking.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """
}
