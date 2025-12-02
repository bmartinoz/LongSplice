// =======================================
// SUBWORKFLOW: DO_ALIGNMENT_STATS
// =======================================

nextflow.enable.dsl = 2

include { SAMTOOLS_FLAGSTAT } from '../modules/samtools_flagstat.nf'
include { SAMTOOLS_IDXSTATS } from '../modules/samtools_idxstats.nf'
include { SAMTOOLS_STATS    } from '../modules/samtools_stats.nf'

workflow DO_ALIGNMENT_STATS {
    take:
        ch_sorted_bam_indexed // Canal: tuple(meta, sorted_bam_path, index_path)
        ch_fasta_path         // Canal: path al FASTA (único elemento)

    main:
        //ch_versions = Channel.empty()

        // 1. Samtools Flagstat
        SAMTOOLS_FLAGSTAT(ch_sorted_bam_indexed)
        ch_flagstat_reports = SAMTOOLS_FLAGSTAT.out.flagstat
        //ch_versions = ch_versions.mix(SAMTOOLS_FLAGSTAT.out.versions.ifEmpty(Channel.empty()))

        // 2. Samtools Idxstats
        SAMTOOLS_IDXSTATS(ch_sorted_bam_indexed)
        ch_idxstats_reports = SAMTOOLS_IDXSTATS.out.idxstats
        //ch_versions = ch_versions.mix(SAMTOOLS_IDXSTATS.out.versions.ifEmpty(Channel.empty()))

        // 3. Samtools Stats
        SAMTOOLS_STATS(ch_sorted_bam_indexed, ch_fasta_path)
        ch_stats_reports = SAMTOOLS_STATS.out.stats
        //ch_versions = ch_versions.mix(SAMTOOLS_STATS.out.versions.ifEmpty(Channel.empty()))

    emit:
        flagstat_reports = ch_flagstat_reports
        idxstats_reports = ch_idxstats_reports
        stats_reports    = ch_stats_reports
        //versions         = ch_versions.unique()
}