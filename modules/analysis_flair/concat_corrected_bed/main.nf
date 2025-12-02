/*
* Concatenate all corrected bed from the FLAIR_CORRECT channel
* and deposit to './collapse/combined_samples.all_corrected.bed'
*/
process CONCAT_CORRECTED_BED {
    tag './collapse/combined_samples.all_corrected.bed'
    publishDir "${params.outdir}/collapse", mode: 'copy'

    input: 
        tuple val(meta), path(corrected_files)


    output: 
        tuple val(meta), path("combined_samples.all_corrected.bed"), emit: combined_bed


    script:
    def bed_files = corrected_files.findAll { it.toString().endsWith('flair_all_corrected.bed') }.join(' ')
    """
    cat ${bed_files} > combined_samples.all_corrected.bed
    """

}