/*
* flair collapse for isoform identification
*
* The important output files are
*  - prefix.isoforms.gtf - your custom transcriptome which you can align to if you want
*  - prefix.isoforms.bed - the easiest way to visualize your isoforms on the UCSC genome browser or IGV, can also be useful for FLAIR-quantify
*  - prefix.combined.isoform.read.map.txt - all detected isoforms associated with the reads that support them
*/

process FLAIR_COLLAPSE {
    tag './collapse/combined_samples.flair.collapse.*'
    publishDir "${params.outdir}/collapse", mode: 'copy'
    label 'process_high'

    input:
        path(ref_fasta)
        path(gtf)
        tuple val(meta), path(combined_fastq)
        tuple val(meta_bed), path(combined_corrected_bed)

    output:
        tuple val(meta), path('*flair.collapse.isoforms.fa'), emit: isoforms_fa
        tuple val(meta), path('*flair.collapse.isoforms.bed'), emit: isoforms_bed
        tuple val(meta), path('*flair.collapse.isoform.read.map.txt'), emit: isoform_read_map
        tuple val(meta), path('*flair.collapse.isoforms.gtf'), emit: isoforms_gtf      

    script: 
    """
    flair collapse \
        -g ${ref_fasta} \
        --gtf ${gtf} \
        -q ${combined_corrected_bed} \
        -r ${combined_fastq} \
        --annotation_reliant generate \
        --generate_map \
        --check_splice \
        --stringent \
        --output ${meta.id}_flair.collapse \
        --threads $task.cpus
    """
}