// =======================================
// MODULE: SAMTOOLS_INDEX
// =======================================

nextflow.enable.dsl = 2

process SAMTOOLS_INDEX {
    tag "${meta.id}"
    //label 'process_low'

    publishDir { "${params.outdir}/alignment/${meta.id}" }, mode: 'copy',
        pattern: "*.{bai,csi,crai}",
        saveAs: { filename -> filename }

    input:
    tuple val(meta), path(bam)

    output:
    // Salida OBLIGATORIA: (meta, bam, bai) para BAMs
    tuple val(meta), path(bam), path("${bam}.bai"), emit: indexed
    path  "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def th = task.cpus ?: 2
    def args = task.ext.args ?: ''

    """
    samtools index -@ ${th} ${args} ${bam}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    """
    # Stub: crea un .bai vacío
    touch ${bam}.bai
    echo "${task.process}:" > versions.yml
    echo "  samtools: stub" >> versions.yml
    """
}
