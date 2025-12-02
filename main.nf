#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
========================================================================================
    LongSplice Pipeline
========================================================================================
*/

params.generateevents_boundary_mode = params.generateevents_boundary_mode ?: ''

include { DO_INPUT }                 from './subworkflows/do_input.nf'
include { DO_GENOME }                from './subworkflows/do_genome.nf'
include { PREPROCESS_ONT }           from './subworkflows/preprocess_ont.nf'
include { DO_ALIGNMENT }             from './subworkflows/do_alignment.nf'
include { DO_SORT_INDEX }            from './subworkflows/do_sort_index.nf'
include { DO_COVERAGE }              from './subworkflows/do_coverage.nf'
include { DO_BIGBED }                from './subworkflows/do_bigbed.nf'
include { DO_ALIGNMENT_STATS }       from './subworkflows/do_alignment_stats.nf'
include { QUANTIFICATION_BAMBU }     from './subworkflows/quantification_bambu.nf'
include { ANALYSIS_DTU }             from './subworkflows/analysis_dtu.nf'
include { ANALYSIS_DESEQ2 }          from './subworkflows/analysis_deseq2.nf'
include { MAKE_CONTRASTSHEET }       from './subworkflows/make_contrastsheet.nf'
include { ANALYSIS_FLAIR }           from './subworkflows/analysis_flair.nf'
include { ANALYSIS_SUPPA }           from './subworkflows/analysis_suppa.nf'
include { MAKE_DESIGN_TSV }          from './subworkflows/make_design_tsv.nf'
include { DO_REPORTING }             from './subworkflows/do_reporting.nf'


/* ===========================
   DEFAULTS
   =========================== */

params.samplesheet                  = params.samplesheet ?: ''
params.quantification_method        = params.quantification_method ?: 'bambu'

params.alpha                        = params.alpha ?: 0.05
params.lfc_threshold                = params.lfc_threshold ?: 1
params.lfc_shrink_method            = params.lfc_shrink_method ?: ''

params.deseq2_contrast_variable     = params.deseq2_contrast_variable ?: 'group'
params.deseq2_reference_level       = params.deseq2_reference_level   ?: ''
params.deseq2_target_level          = params.deseq2_target_level      ?: ''
params.deseq2_level                 = params.deseq2_level ?: 'gene'   // 'gene' o 'transcript'

params.run_deseq2                   = params.run_deseq2 ?: true
params.run_dexseq                   = params.run_dexseq ?: true
params.run_suppa                    = params.run_suppa ?: false
params.run_suppa_diff               = params.run_suppa_diff ?: false   // activa diffSplice si hay design.tsv

params.skip_preprocess_ont          = params.skip_preprocess_ont ?: false
params.skip_quantification          = params.skip_quantification ?: false
params.skip_dtu                     = params.skip_dtu ?: false
params.run_flair                    = params.run_flair ?: false
params.skip_chopper                 = params.skip_chopper ?: false
params.skip_coverage                = params.skip_coverage ?: false
params.skip_bigbed                  = params.skip_bigbed ?: false
// Importante: `skip_alignment` se define en nextflow.config; aquí solo lo consumimos.

workflow {

    if ( !file(params.samplesheet).exists() ) {
        exit 1, "Samplesheet no encontrado: ${params.samplesheet}"
    }

    // 1) INPUT
    DO_INPUT( Channel.value( file(params.samplesheet) ) )
    def ch_samples             = DO_INPUT.out.samples             // tuple(metaMap, fastq)
    def ch_raw_validated_sheet = DO_INPUT.out.validated_samplesheet

    // 2) GENOMA (siempre requerido)
    DO_GENOME(ch_samples)
    def ch_fasta_path  = DO_GENOME.out.fasta
    def ch_gtf_path    = DO_GENOME.out.gtf
    def ch_fai_path    = DO_GENOME.out.fai
    def ch_chrom_sizes = DO_GENOME.out.chrom_sizes

    // =========================
    // BRANCH según skip_alignment
    // =========================
    if ( params.skip_alignment as boolean ) {
        // -----------------------------
        // MODO SOLO FLAIR
        // -----------------------------
        if ( !(params.run_flair as boolean) ) {
            exit 1, "[Main] --skip_alignment true requiere --run_flair true para ejecutar algo."
        }

        log.info "[Main] skip_alignment=true → Solo DO_INPUT → DO_GENOME → PREPROCESS_ONT → ANALYSIS_FLAIR"

        // Forzar PREPROCESS_ONT en este modo (ignoramos params.skip_preprocess_ont)
        PREPROCESS_ONT(ch_samples)
        def ch_nanoplot_reports_pre  = PREPROCESS_ONT.out.nanoplot_reports_pre
        def ch_nanoplot_reports_post = PREPROCESS_ONT.out.nanoplot_reports_post
        def ch_preprocess_logs       = PREPROCESS_ONT.out.chopper_logs
        def ch_cleaned_reads      = PREPROCESS_ONT.out.cleaned_reads  // tuple(metaMap, fastq)

        // Junctions opcionales (canal normalizado)
        def ch_junctions_input = params.shortread_junctions
            ? Channel.fromPath(params.shortread_junctions, checkIfExists: true)
            : Channel.value(file('NO_JUNCTIONS'))

        // Llamada a FLAIR (firma usa refs + lecturas + samplesheet + junctions)
        ANALYSIS_FLAIR(
            ch_cleaned_reads,
            ch_fasta_path,
            ch_fai_path,
            ch_gtf_path

        )

    } else {
        // -----------------------------
        // MODO PIPELINE COMPLETO
        // -----------------------------

        // 3) PREPROCESAMIENTO ONT (respetando tu toggle)
        def ch_reads_for_alignment   = Channel.empty()
        def ch_nanoplot_reports_pre  = Channel.empty()
        def ch_nanoplot_reports_post = Channel.empty()
        def ch_preprocess_logs       = Channel.empty()

        if (!params.skip_preprocess_ont) {
            PREPROCESS_ONT(ch_samples)
            ch_reads_for_alignment   = PREPROCESS_ONT.out.cleaned_reads   // tuple(metaMap, fastq)
            ch_nanoplot_reports_pre  = PREPROCESS_ONT.out.nanoplot_reports_pre
            ch_nanoplot_reports_post = PREPROCESS_ONT.out.nanoplot_reports_post
            ch_preprocess_logs       = PREPROCESS_ONT.out.chopper_logs
        } else {
            log.info "[Main] Omitiendo PREPROCESS_ONT."
            ch_reads_for_alignment = ch_samples  // ya es tuple(metaMap, fastq)
        }

        // 4) ALINEAMIENTO
        DO_ALIGNMENT(ch_fasta_path, ch_reads_for_alignment)
        def ch_bam_aligned = DO_ALIGNMENT.out.align_bam   // tuple(metaMap, bam)

        // 5) SORT + INDEX
        DO_SORT_INDEX( ch_bam_aligned )
        def ch_sorted_shared = DO_SORT_INDEX.out.sortbam   // tuple(metaMap, bam, bai)

        // Reusar el mismo canal
        def ch_sorted_for_cov       = ch_sorted_shared
        def ch_sorted_for_stats     = ch_sorted_shared
        def ch_sorted_for_bambu_raw = ch_sorted_shared

        // Normaliza para Bambu (meta,bam,bai) — robusto si falta .bai
        def ch_sorted_for_bambu = ch_sorted_for_bambu_raw.map { t ->
            def meta = t[0]
            def bam  = t[1]
            def bai  = (t.size() > 2) ? t[2] : null
            if (bai == null) {
                def cand1 = file("${bam}.bai")
                def cand2 = file(bam.toString().replaceAll(/\.bam$/, '.bai'))
                if      (cand1.exists()) bai = cand1
                else if (cand2.exists()) bai = cand2
                else error "[BAMBU] No se encontró índice (.bai) para ${bam}"
            }
            tuple(meta, bam, bai)
        }

        // 6) BIGWIG (opcional)
        if ( !params.skip_coverage ) {
            DO_COVERAGE(ch_sorted_for_cov)
            def ch_bigwig_files = DO_COVERAGE.out.bigwig
        } else {
            log.info "[Main] BIGWIG coverage omitido."
        }

        // 6.b) BIGBED (opcional)
        if (!params.skip_bigbed) {
            DO_BIGBED( ch_sorted_for_cov, ch_chrom_sizes )
            def ch_bigbed_files = DO_BIGBED.out.bigbed
        } else {
            log.info "[Main] BIGBED omitido."
        }

        // 7) STATS
        DO_ALIGNMENT_STATS(ch_sorted_for_stats, ch_fasta_path)
        def ch_flagstat_reports = DO_ALIGNMENT_STATS.out.flagstat_reports
        def ch_idxstats_reports = DO_ALIGNMENT_STATS.out.idxstats_reports
        def ch_stats_reports    = DO_ALIGNMENT_STATS.out.stats_reports

        // 8) Bambu (+ DTU + DE + SUPPA)
        if (!params.skip_quantification) {
            log.info "[Main] Ejecutando cuantificación con Bambu."
            QUANTIFICATION_BAMBU(
                ch_sorted_for_bambu,
                ch_fasta_path,
                ch_gtf_path,
                ch_fai_path
            )
            QUANTIFICATION_BAMBU.out.counts_transcript.set { ch_tx_counts }
            QUANTIFICATION_BAMBU.out.counts_gene      .set { ch_gene_counts }
            QUANTIFICATION_BAMBU.out.gtf_extended     .set { ch_gtf_ext }
            def ch_sheet_val = ch_raw_validated_sheet

            if (!params.skip_dtu && params.run_dexseq) {
                ANALYSIS_DTU(ch_tx_counts, ch_sheet_val, ch_gtf_ext)
            } else {
                log.info "[Main] DEXSeq omitido."
            }

            if (params.run_deseq2) {
                def ch_contrast_var = Channel.value( params.deseq2_contrast_variable ?: 'group' )
                def ch_ref_level    = Channel.value( params.deseq2_reference_level   ?: '' )
                def ch_tgt_level    = Channel.value( params.deseq2_target_level      ?: '' )
                def ch_counts_for_deseq2 = (params.deseq2_level?.toLowerCase() == 'transcript') ? ch_tx_counts : ch_gene_counts

                ANALYSIS_DESEQ2(
                    ch_counts_for_deseq2,
                    ch_sheet_val,
                    ch_contrast_var,
                    ch_ref_level,
                    ch_tgt_level,
                    ch_gtf_ext
                )
            } else {
                log.info "[Main] DESeq2 omitido."
            }

            if (params.run_suppa) {
                MAKE_DESIGN_TSV( ch_raw_validated_sheet )
                def ch_design_tsv = MAKE_DESIGN_TSV.out.design
                ANALYSIS_SUPPA(ch_gtf_ext, ch_tx_counts, ch_design_tsv)
            } else {
                log.info "[Main] SUPPA2 omitido."
            }
        } else {
            log.info "[Main] Cuantificación (bambu) omitida."
            if (params.run_suppa) {
                def path_counts = file("${params.outdir ?: 'results'}/bambu/counts_transcript.txt")
                def path_gtf    = file("${params.outdir ?: 'results'}/bambu/extended_annotations.gtf")
                if (!path_counts.exists() || !path_gtf.exists()) {
                    exit 1, "[SUPPA] No se encontraron artefactos de bambu: ${path_counts} / ${path_gtf}"
                }
                MAKE_DESIGN_TSV( ch_raw_validated_sheet )
                def ch_design_tsv2 = MAKE_DESIGN_TSV.out.design
                ANALYSIS_SUPPA(Channel.value(path_gtf), Channel.value(path_counts), ch_design_tsv2)
            } else {
                log.info "[Main] SUPPA2 omitido."
            }
        }

        // 9) FLAIR (independiente)
        if (params.run_flair) {
            log.info "[Main] Ejecutando análisis FLAIR2."
            def ch_junctions_input = params.shortread_junctions
                ? Channel.fromPath(params.shortread_junctions, checkIfExists: true)
                : Channel.value(file('NO_JUNCTIONS'))


            ANALYSIS_FLAIR(
                ch_samples,
                ch_fasta_path,
                ch_fai_path,
                ch_gtf_path


            )
        } else {
            log.info "[Main] FLAIR2 omitido."
        }
    }

    // (Opcional) DO_REPORTING aquí si quieres consolidar reportes globales
}
