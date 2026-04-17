// re-annotate a gffcompare-annotated GTF so that exact-match transcripts
// recover their reference IDs and novel transcripts get a tool-specific prefix
process REANNOTATEGTF {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/biopython:v250501"

    input:
    tuple val(meta), path(annotated_gtf)

    output:
    tuple val(meta), path("*.reannotated.gtf"),    emit: gtf
    tuple val(meta), path("*.id_mapping.tsv"),     emit: mapping
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    reannotate_gtf.py \\
        $args \\
        --mapping ${prefix}.id_mapping.tsv \\
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
    touch ${prefix}.id_mapping.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        polars: \$(python -c "import polars; print(polars.__version__)")
        gtfparse: \$(python -c "import gtfparse; print(gtfparse.__version__)")
    END_VERSIONS
    """
}
