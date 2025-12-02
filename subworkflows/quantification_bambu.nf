// =======================================
// SUBWORKFLOW: QUANTIFICATION_BAMBU
// =======================================

nextflow.enable.dsl = 2

include { BAMBU } from '../modules/bambu.nf'

workflow QUANTIFICATION_BAMBU {
    take:
      ch_indexed_bam  // (meta, bam) o (meta, bam, bai)
      ch_fasta
      ch_gtf
      ch_fai          // se toma pero NO se usa aquí

    main:
      //ch_versions = Channel.empty()

      // 1) Normaliza a obtener SOLO el BAM (ignora meta y bai para Bambu)
      ch_bams_only = ch_indexed_bam.map { t ->
        // t = [meta, bam] o [meta, bam, bai]
        t[1]
      }

      // 2) Colecta todos los BAMs en una única lista (List<Path>)
      ch_bams_list = ch_bams_only.collect()

      // 3) Construye el canal requerido por BAMBU: (fasta, gtf)
      ch_ref = ch_fasta.combine(ch_gtf).map { f, g -> tuple(f, g) }

      // 4) Ejecuta BAMBU
      BAMBU(
        ch_ref,
        ch_bams_list,
        params.bambu_discovery,
        params.bambu_ndr,
        params.bambu_opts
      )

      //ch_versions = ch_versions.mix( BAMBU.out.versions.ifEmpty(Channel.empty()) )

    emit:
      counts_gene       = BAMBU.out.ch_gene_counts
      counts_transcript = BAMBU.out.ch_transcript_counts
      gtf_extended      = BAMBU.out.extended_gtf
      //versions          = ch_versions.unique()
}
