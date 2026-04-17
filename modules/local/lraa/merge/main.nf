// Merge per-sample LRAA GTFs into a cohort-level GTF
process LRAA_MERGE {
    tag "merge"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "us-central1-docker.pkg.dev/methods-dev-lab/lraa/lraa:v0.15.0"

    input:
    path gtfs, stageAs: 'sample_gtfs/*'
    path ref_genome

    output:
    tuple val([id: 'merge']), path("cohort.LRAA.merged.gtf"), emit: gtf
    path ("cohort.LRAA.merged.gtf.tracking.tsv"), emit: tracking
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    merge_LRAA_GTFs.py \\
        --genome ${ref_genome} \\
        --gtf sample_gtfs/*.gtf \\
        --output_gtf cohort.LRAA.merged.gtf \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """

    stub:
    """
    touch cohort.LRAA.merged.gtf
    touch cohort.LRAA.merged.gtf.tracking.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        lraa: \$(LRAA --version 2>&1 | head -1)
    END_VERSIONS
    """
}
