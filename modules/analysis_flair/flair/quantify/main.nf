// modules/flair/quantify/main.nf
nextflow.enable.dsl = 2

/*
* flair quantify process that creates flair.quantify.*.isoform.read.map.txt and
* flair.quanitfy.counts.tsv files in the ./quant foler
* NOTE: the parameters suggested here are from the tutorial and for hg38
*/
process FLAIR_QUANTIFY {
    tag  "flair_quantify"
    publishDir "${params.outdir}/quant", mode: 'copy'
    label 'process_high'
    

    input:
        tuple val(meta), path(collapse_bed)
        tuple val(meta2), path(collapse_fa)
        path(sample_manifest_tsv)


    output:
        path('flair.quantify.*.isoform.read.map.txt'), emit: quantify_isoform_read_map
        path('flair.quantify.counts.tsv'),            emit: quantify_counts

    script:
    """
    flair quantify \
        -r ${sample_manifest_tsv} \
        -i ${collapse_fa} \
        --generate_map \
        --isoform_bed ${collapse_bed} \
        --stringent \
        --check_splice \
        --threads $task.cpus \
        --output flair.quantify    
    """        
}