// subworkflows/analysis_dtu.nf
nextflow.enable.dsl=2

include { DEXSEQ } from '../modules/dexseq.nf'

workflow ANALYSIS_DTU {
    take:
    counts_transcript
    validated_samplesheet
    extended_gtf

    main:
    DEXSEQ(
        counts_transcript,
        validated_samplesheet,
        extended_gtf
    )

    //ch_versions = DEXSEQ.out.versions.ifEmpty(Channel.empty())

    emit:
    dexseq_results = DEXSEQ.out.results_txt
    log_file       = DEXSEQ.out.log_file
    //versions       = ch_versions
}
