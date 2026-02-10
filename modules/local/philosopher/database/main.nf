// download fasta files from uniprot with philosopher
process PHILOSOPHER_DATABASE {
    tag "${meta.id}_download"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "quay.io/shahlab_singularity/fragpipe:23.1"

    input:
    val meta
    val reviewed
    val isoforms

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.fas"), emit: fasta
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reviewed_arg = reviewed == true ? '--reviewed' : ''
    def isoforms_arg = isoforms == true ? '--isoform' : ''
    def philosopher = '/fragpipe_bin/fragpipe-23.1/fragpipe-23.1/tools/Philosopher/philosopher-v5.1.2'
    """
    ${philosopher} workspace --init --nocheck
    ${philosopher} database --id ${meta.id} --nodecoys ${reviewed_arg} ${isoforms_arg} ${args}
    ${philosopher} workspace --clean --nocheck
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        philosopher: 5.1.2
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.fas

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        philosopher: 5.1.2
    END_VERSIONS
    """
}
