// modules/flair/diffsplice/main.nf
nextflow.enable.dsl = 2

/*
 * FLAIR diffsplice - Alternative splicing analysis
 * Detects: intron retention (ir), alt 3'/5' splicing, cassette exons (es)
 * Optional: DRIMSeq statistical testing with 3+ replicates per condition
 */
process FLAIR_DIFFSPLICE {
    tag "flair_diffsplice"
    publishDir "${params.outdir}/diffsplice", mode: 'copy'
    label 'process_medium'
    
    input:
    tuple val(meta), path(isoforms_bed)            // combined_samples.flair.collapse.isoforms.bed
    path(counts_matrix)           // flair.quantify.counts.tsv
    
    output:
    path("diffsplice_out"), emit: diffsplice
    
    script:

    """
    mkdir -p diffsplice_out
    
    flair diffsplice \
        --isoforms ${isoforms_bed} \
        --counts_matrix ${counts_matrix} \
        --out_dir diffsplice_out \
        --threads $task.cpus \
        --out_dir_force \
        --test
    """
}