nextflow.enable.dsl = 2

include { DESEQ2 as DESEQ2_RUN } from '../modules/deseq2.nf'

/*
 * Subworkflow: ANALYSIS_DESEQ2
 */

workflow ANALYSIS_DESEQ2 {

  take:
    counts_matrix
    samplesheet
    contrast_variable
    reference_level
    target_level
    gtf

  main:
    samplesheet_ch      = samplesheet.ifEmpty( Channel.value('') )
    contrast_var_ch     = contrast_variable.ifEmpty( Channel.value('group') )
    reference_level_ch  = reference_level.ifEmpty( Channel.value('') )
    target_level_ch     = target_level.ifEmpty( Channel.value('') )
    gtf_ch              = gtf.ifEmpty( Channel.value('') )

    DESEQ2_RUN(
      counts_matrix,
      samplesheet_ch,
      contrast_var_ch,
      reference_level_ch,
      target_level_ch,
      gtf_ch
    )

    results_csv       = DESEQ2_RUN.out.results_csv
    deseq2_txt        = DESEQ2_RUN.out.deseq2_txt
    normalized_counts = DESEQ2_RUN.out.normalized_counts
    plots             = DESEQ2_RUN.out.plots
    rds               = DESEQ2_RUN.out.rds
    summary           = DESEQ2_RUN.out.summary
    deseq2_version    = DESEQ2_RUN.out.versions

  emit:
    results_csv
    deseq2_txt
    normalized_counts
    plots
    rds
    summary
    deseq2_version
}
