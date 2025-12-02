// =======================================
// MODULE: SAMTOOLS_SORT
// =======================================

nextflow.enable.dsl = 2

process SAMTOOLS_SORT {
    tag "$meta.id"
    //label 'process_medium'

    publishDir "${params.outdir}/alignment/${meta.id}", mode: 'copy',
        pattern: "${ext.prefix ?: meta.id + '.sorted'}.bam",
        saveAs: { filename -> filename }

    input:
    tuple val(meta), path(bam), val(ext)

    output:
    tuple val(meta), path("${ext.prefix ?: meta.id + '.sorted'}.bam"), emit: bam
    path "versions.yml", emit: versions

    script:
    def args   = ext.args ?: ''
    def prefix = ext.prefix ?: "${meta.id}.sorted"

    if ("${bam.baseName}" == "${prefix}") {
       error "Input BAM base name and output prefix are the same. Use 'ext.prefix'."
    }

    """
    samtools sort $args -@ $task.cpus -o ${prefix}.bam -T ${prefix}.tmp $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    def prefix_stub = ext.prefix ?: "${meta.id}.sorted"
    """
    touch ${prefix_stub}.bam
    echo "${task.process}:" > versions.yml
    echo "  samtools: stub" >> versions.yml
    """
}
