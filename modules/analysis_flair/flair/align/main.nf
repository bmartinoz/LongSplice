/*
* flair align: aligns fastq files to using minmap2 and converts the bam file
* to a bed file
* NOTE: flair align should allow hifi reads and use PacBio's wrapper function of minmap2
*/
process FLAIR_ALIGN {
    tag "$meta.id"
    label 'process_medium'

    input:
    tuple val(meta), path(fastq)
    path ref_fasta
    path ref_index
    path gtf

    output:
    tuple val(meta), path("*.bam"), path("*.bam.bai")        , emit: bam
    tuple val(meta), path("*.bed")                           , emit: bed
    path "versions.yml"                                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args_flair_align ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    flair align \\
        -g ${ref_fasta} \\
        -r ${fastq} \\
        -o ${prefix}.flair.aligned \\
        --quality ${params.min_mapq} \\
        --threads ${task.cpus} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flair: \$( flair --version | sed 's/flair //' )
    END_VERSIONS
    """

}