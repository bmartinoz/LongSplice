// =======================================
// MODULE: SAMTOOLS_IDXSTATS
// =======================================

nextflow.enable.dsl = 2

process SAMTOOLS_IDXSTATS {
    tag "$meta.id"
    label 'process_low'

    publishDir "${params.outdir}/alignment_stats/samtools_idxstats/${meta.id}", mode: 'copy', pattern: "*.idxstats", saveAs: { fn -> fn }

    // Directivas conda y container eliminadas

    input:
    tuple val(meta), path(bam), path(bai)

    output:
    tuple val(meta), path("*.idxstats"), emit: idxstats
    path  "versions.yml"               , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    samtools \\
        idxstats \\
        $bam \\
        > ${prefix}.idxstats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}