nextflow.enable.dsl = 2

include { NANOPLOT as NANOPLOT_PRE_TRIM }  from '../modules/nanoplot.nf'
include { NANOPLOT as NANOPLOT_POST_TRIM } from '../modules/nanoplot.nf'
include { PYCHOPPER }                      from '../modules/pychopper.nf'
include { CHOPPER   }                      from '../modules/chopper.nf'
include { RESTRANDER }                     from '../modules/restrander.nf'

workflow PREPROCESS_ONT {
    take:
    // tuple val(meta_map_o_id), path(reads_path)
    samples

    main:
    // -------- PRE NanoPlot (opcional)
    def np_pre_reports = Channel.empty()
    if (params.run_nanoplot) {
        NANOPLOT_PRE_TRIM(samples, "PRE_TRIM_Report")
        np_pre_reports = NANOPLOT_PRE_TRIM.out.report_dir
    }

    // -------- PYCHOPPER (opcional)
    def ch_for_chopper = Channel.empty()
    def ch_pychopper_reports_coll      = Channel.empty().collect()
    def ch_pychopper_stats_coll        = Channel.empty().collect()
    def ch_pychopper_unclassified_coll = Channel.empty().collect()
    def ch_pychopper_rescued_coll      = Channel.empty().collect()

    if (!params.skip_pychopper) {
        def py_args_final = (params.pychopper_args ?: "").trim()
        def py_prefixes   = samples.map { it[0] instanceof Map ? it[0].id : it[0] }

        PYCHOPPER(samples, py_args_final, py_prefixes)

        ch_for_chopper                 = PYCHOPPER.out.fastq
        ch_pychopper_reports_coll      = PYCHOPPER.out.report.collect().ifEmpty([])
        ch_pychopper_stats_coll        = PYCHOPPER.out.stats.collect().ifEmpty([])
        ch_pychopper_unclassified_coll = PYCHOPPER.out.unclassified.collect().ifEmpty([])
        ch_pychopper_rescued_coll      = PYCHOPPER.out.rescued.collect().ifEmpty([])
    } else {
        ch_for_chopper                 = samples
        ch_pychopper_reports_coll      = ch_pychopper_reports_coll.ifEmpty([])
        ch_pychopper_stats_coll        = ch_pychopper_stats_coll.ifEmpty([])
        ch_pychopper_unclassified_coll = ch_pychopper_unclassified_coll.ifEmpty([])
        ch_pychopper_rescued_coll      = ch_pychopper_rescued_coll.ifEmpty([])
    }

    // -------- CHOPPER (opcional)
    def ch_processed_reads = Channel.empty()
    def ch_chopper_logs_collected = Channel.empty().collect()

    if (!params.skip_chopper) {
        def chopper_args_final   = params.chopper_args ?: ""
        def chopper_prefixes_val = ch_for_chopper.map { it[0] instanceof Map ? it[0].id : it[0] }

        CHOPPER(ch_for_chopper, chopper_args_final, chopper_prefixes_val)

        ch_processed_reads        = CHOPPER.out.chopped_reads
        ch_chopper_logs_collected = CHOPPER.out.chopper_logs.collect().ifEmpty([])
    } else {
        ch_processed_reads        = ch_for_chopper
        ch_chopper_logs_collected = ch_chopper_logs_collected.ifEmpty([])
    }

    // -------- RESTRANDER (opcional) — entre CHOPPER y NanoPlot POST
    def ch_restranded_reads = Channel.empty()
    def ch_restrander_logs  = Channel.empty().collect()

    if (!params.skip_restrander) {
        def restrander_prefixes   = ch_processed_reads.map { it[0] instanceof Map ? it[0].id : it[0] }
        def restrander_config_val = params.restrander_config   // pasar como val, no como path

        RESTRANDER(
            ch_processed_reads,       // tuple(meta, fastq)
            restrander_config_val,    // val cfg_path
            restrander_prefixes       // val prefix
        )

        ch_restranded_reads = RESTRANDER.out.restranded_reads
        ch_restrander_logs  = RESTRANDER.out.restrander_log.collect().ifEmpty([])
    } else {
        ch_restranded_reads = ch_processed_reads
        ch_restrander_logs  = ch_restrander_logs.ifEmpty([])
    }

    // -------- POST NanoPlot (opcional) --------
    def np_post_reports = Channel.empty()
    if (params.run_nanoplot) {
        NANOPLOT_POST_TRIM(ch_restranded_reads, "POST_TRIM_Report")
        np_post_reports = NANOPLOT_POST_TRIM.out.report_dir
    }

    emit:
    cleaned_reads               = ch_restranded_reads
    nanoplot_reports_pre        = np_pre_reports
    nanoplot_reports_post       = np_post_reports
    pychopper_reports_out       = ch_pychopper_reports_coll
    pychopper_stats_out         = ch_pychopper_stats_coll
    pychopper_unclassified_out  = ch_pychopper_unclassified_coll
    pychopper_rescued_out       = ch_pychopper_rescued_coll
    chopper_logs                = ch_chopper_logs_collected
    restrander_logs             = ch_restrander_logs
}
