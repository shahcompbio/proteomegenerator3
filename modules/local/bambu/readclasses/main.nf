// create read classes with bambu
process BAMBU_READCLASSES {
    tag "${meta.id}"
    // temporary label for testing
    label 'process_high_memory'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/bambu:HongYhong_fix"

    input:
    tuple val(meta), path(bam)
    val yieldsize
    path ref_genome
    path ref_gtf

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.rds"), emit: rds
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    create_readclasses.R \\
        --bam=${bam} \\
        --yieldsize=${yieldsize} \\
        --ref_genome=${ref_genome} \\
        --ref_gtf=${ref_gtf}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bambu/HongYhong_fix: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.rds

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        r-base: \$(echo \$(R --version 2>&1) | sed 's/^.*R version //; s/ .*\$//')
        bambu/HongYhong_fix: \$(Rscript -e "library(bambu); cat(as.character(packageVersion('bambu')))")
    END_VERSIONS
    """
}
