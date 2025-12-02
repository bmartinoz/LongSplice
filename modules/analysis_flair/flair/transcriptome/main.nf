// modules/flair/transcriptome/main.nf
nextflow.enable.dsl = 2

/*
 * FLAIR_TRANSCRIPTOME: genera transcriptoma de alta confianza directamente desde BAM
 * (más rápido y con menos memoria que correct + collapse).
 * Recibe meta como Map {id: ...} y usa meta.id para prefijos seguros.
 */

process FLAIR_TRANSCRIPTOME {
    tag { meta.id }
    label 'process_high'

    cpus   { params.flair_transcriptome_cpus ?: 8 }
    memory { params.flair_transcriptome_mem  ?: '24 GB' }
    time   { params.flair_transcriptome_time ?: '36h' }

    publishDir (params.flair_transcriptome_publishdir ?: "${params.outdir ?: 'results'}/flair/transcriptome"),
               mode: 'copy',
               overwrite: true

    input:
    // Firma alineada a ANALYSIS_FLAIR: (metaMap, bam, bai) + refs + junctions
    tuple val(meta), path(sample_bam), path(sample_bai)
    path genome_fasta
    path gtf
    path shortread_junctions   // puede venir como marcador 'NO_JUNCTIONS'

    output:
    // Usa meta.id en nombres de salida (no interpolar el Map completo).
    tuple val(meta), path("${meta.id}.flair.isoforms.bed"),         emit: bed
    tuple val(meta), path("${meta.id}.flair.isoforms.gtf"),         emit: gtf
    tuple val(meta), path("${meta.id}.flair.isoforms.fa"),          emit: fa
    tuple val(meta), path("${meta.id}.flair.isoforms.cds.fa"),      emit: cds
    tuple val(meta), path("${meta.id}.flair.read.map.txt"),         emit: readmap, optional: true

    script:
    // Prefijo seguro
    def sid = (meta instanceof Map ? meta.id : meta.toString())

    // Flags opcionales
    def gtf_arg        = (gtf && gtf.name != 'NO_GTF') ? "-f ${gtf}" : ''
    def junctions_arg  = (shortread_junctions && shortread_junctions.name != 'NO_JUNCTIONS') ? "-j ${shortread_junctions}" : ''

    // Parámetros tunables (con defaults sensatos)
    def support        = params.flair_transcriptome_support          ?: 3
    def ss_window      = params.flair_transcriptome_ss_window        ?: 15
    def end_window     = params.flair_transcriptome_end_window       ?: 100
    def max_ends       = params.flair_transcriptome_max_ends         ?: 2
    def stringent_flag = params.flair_transcriptome_stringent        ?  '--stringent'   : ''
    def check_splice   = params.flair_transcriptome_check_splice     ?  '--check_splice': ''
    def no_align_annot = params.flair_transcriptome_noaligntoannot   ?  '--noaligntoannot' : ''
    def no_redundant   = params.flair_transcriptome_no_redundant     ?: 'none'
    def filter_mode    = params.flair_transcriptome_filter           ?: 'default'
    def predict_cds    = params.flair_transcriptome_predictCDS       ?  '--predictCDS'  : ''

    """
    set -euo pipefail

    flair transcriptome \\
        -b ${sample_bam} \\
        -g ${genome_fasta} \\
        -o ${sid}.flair \\
        -t ${task.cpus} \\
        ${gtf_arg} \\
        ${junctions_arg} \\
        --ss_window ${ss_window} \\
        -s ${support} \\
        ${stringent_flag} \\
        ${check_splice} \\
        -w ${end_window} \\
        ${no_align_annot} \\
        -n ${no_redundant} \\
        --max_ends ${max_ends} \\
        --filter ${filter_mode} \\
        ${predict_cds}

    """
}
