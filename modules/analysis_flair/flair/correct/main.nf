/*
* flair correctr: corrects alignments to the annotated splice sites
*/ 
process FLAIR_CORRECT {
    tag "${meta.id}"
    publishDir "${params.outdir}/correct", mode: 'copy'
    label 'process_medium'

    input: 
        tuple val(meta), path(bed)
        path ref_fasta
        path gtf

    output: 
        tuple val(meta), path("${meta.id}.flair_all_corrected.bed"),         emit: corrected_bed
        tuple val(meta), path("${meta.id}.flair_all_inconsistent.bed"),      emit: inconsistent_bed
        path "versions.yml",                                                 emit: versions      

    script:
    """
    flair correct  \
        -f ${gtf} \
        -q ${bed} \
        -o ${meta.id}.flair \
        --threads $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flair: \$( flair --version | sed 's/flair //' )
    END_VERSIONS
    """
}