// Merge per-sample LRAA .quant.expr files into a cohort expression matrix
process LRAA_QUANTMERGE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/biopython:v250501"

    input:
    tuple val(meta), path(quant_files, stageAs: 'quant_files/*')

    output:
    tuple val(meta), path("*.expr_matrix.tsv"), emit: matrix
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    merge_lraa_quant.py \\
        --input quant_files/*.quant.expr \\
        --output ${prefix}.expr_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.expr_matrix.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
