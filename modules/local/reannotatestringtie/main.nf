// re-annotate stringtie output in a sane way; otherwise we don't know
// which transcripts belong to which genes, or which are canonical vs non-canonical
process REANNOTATESTRINGTIE {
    tag "${meta.id}"
    label 'process_single'

    // TODO nf-core: See section in main README for further information regarding finding and adding container addresses to the section below.
    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/biopython:v250501"

    input:
    tuple val(meta), path(annotated_gtf)

    output:
    tuple val(meta), path("*.reannotated.gtf"), emit: gtf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reannotate_gtf.py \\
        ${annotated_gtf} \\
        ${prefix}.reannotated.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        polars: \$(python -c "import polars; print(polars.__version__)")
        gtfparse: \$(python -c "import gtfparse; print(gtfparse.__version__)")
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reannotated.gtf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        polars: \$(python -c "import polars; print(polars.__version__)")
        gtfparse: \$(python -c "import gtfparse; print(gtfparse.__version__)")
    END_VERSIONS
    """
}
