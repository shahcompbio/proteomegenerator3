// SQANTI3 Quality Control — classify isoforms against reference annotation
// Produces classification table, corrected GTF, and junctions file

process SQANTI3_QC {
    tag "${meta.id}"
    label 'process_high'

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'docker://anaconesalab/sqanti3:v6.0.1'
        : 'docker.io/anaconesalab/sqanti3:v6.0.1'}"

    input:
    tuple val(meta), path(gtf)
    path ref_gtf
    path ref_fasta

    output:
    tuple val(meta), path("${prefix}_classification.txt"), emit: classification
    tuple val(meta), path("${prefix}_corrected.gtf"), emit: corrected_gtf
    tuple val(meta), path("${prefix}_junctions.txt"), emit: junctions
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Fix missing libbz2.so.1 symlink in container (needed by gtfToGenePred)
    mkdir -p lib_fix
    ln -sf /usr/local/lib/libbz2.so.1.0.8 lib_fix/libbz2.so.1
    export LD_LIBRARY_PATH="\$PWD/lib_fix:\$LD_LIBRARY_PATH"

    sqanti3_qc.py \\
        --isoforms ${gtf} \\
        --refGTF ${ref_gtf} \\
        --refFasta ${ref_fasta} \\
        -o ${prefix} \\
        -d . \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: \$(sqanti3_qc.py -v 2>&1 | grep 'SQANTI3' | sed 's/.*SQANTI3 //' || echo "6.0.1")
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_classification.txt
    touch ${prefix}_corrected.gtf
    touch ${prefix}_junctions.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sqanti3: "6.0.1"
    END_VERSIONS
    """
}
