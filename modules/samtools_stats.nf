// =======================================
// MODULE: SAMTOOLS_STATS
// =======================================

nextflow.enable.dsl = 2

process SAMTOOLS_STATS {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/alignment_stats/samtools_stats/${meta.id}", mode: 'copy', pattern: "*.stats", saveAs: { fn -> fn }

    input:
    tuple val(meta), path(input_bam), path(input_index)
    path fasta

    output:
    tuple val(meta), path("*.stats"), emit: stats
    path  "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def user_opts = params.samtools_stats_opts ?: "" // Para opciones personalizadas
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = fasta ? "--reference ${fasta}" : ""
    """
    samtools \\
        stats \\
        --threads ${task.cpus} \\
        ${reference} \\
        ${args} \\
        ${user_opts} \\
        ${input_bam} \\
        > ${prefix}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix_stub = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix_stub}.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: stub
    END_VERSIONS
    """
}