// Apply ID mapping from reannotate_gtf.py to LRAA .quant.expr files
// so that gene_id and transcript_id columns match the reannotated GTF

process LRAA_REANNOTATEQUANT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/biopython:v250501"

    input:
    tuple val(meta), path(quant_expr), path(id_mapping)

    output:
    tuple val(meta), path("*.reannotated.quant.expr"), emit: quant
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    reannotate_quant.py \\
        --quant ${quant_expr} \\
        --mapping ${id_mapping} \\
        --output ${prefix}.reannotated.quant.expr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.reannotated.quant.expr

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
