// modules/flair/diffexp/main.nf
nextflow.enable.dsl = 2

/*
 * FLAIR diffexp - Differential expression and usage analysis
 * Runs DESeq2 (genes & isoforms) and DRIMSeq (isoforms only)
 * Requires: R, DESeq2, DRIMSeq packages
 * NOTE: Control condition should be alphabetically lower than test condition
 */
process FLAIR_DIFFEXP {
    tag "flair_diffexp"
    publishDir "${params.outdir}/diffexp", mode: 'copy'
    label 'process_medium'
    
    input:
    path(counts_matrix)           // flair.quantify.counts.tsv
    
    output:
    path("diffexp_out"), emit: diffexp
    
    script:
    def exp_arg = task.ext.exp_thresh ? task.ext.diffexp_thresh : "1"
    """
    mkdir -p diffexp_out
    
    flair diffexp \
        --counts_matrix ${counts_matrix} \
        --out_dir diffexp_out \
        --exp_thresh ${exp_arg} \
        --threads $task.cpus \
        --out_dir_force
    
    """
}