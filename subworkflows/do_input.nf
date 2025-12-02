// =======================================
// SUBWORKFLOW: DO_INPUT  (corregido con fallback seguro)
// =======================================

nextflow.enable.dsl = 2

include { SAMPLESHEET_INPUT } from '../modules/samplesheet_input.nf'

workflow DO_INPUT {
    take:
      samplesheet_path // Value-channel con el Path al samplesheet (viene de main.nf)

    main:
      //ch_versions = Channel.empty()
      ch_samplesheet = samplesheet_path.map { p ->
          def f = file(p)
          if( !f.exists() )   error "Samplesheet not found: ${f}"
          if( f.size() == 0 ) error "Samplesheet file is empty: ${f}"

          def header = null
          try { header = f.readLines().first() }
          catch (Throwable e) { error "Error reading header from samplesheet ${f}: ${e.message}" }

          if( !header || header.trim().isEmpty() )
              error "Samplesheet file '${f}' seems empty or header could not be read."

          tuple(f, header)
      }

      // Validación
      SAMPLESHEET_INPUT(ch_samplesheet)
      ch_cleaned = SAMPLESHEET_INPUT.out.cleaned
      //ch_versions = ch_versions.mix( SAMPLESHEET_INPUT.out.versions.ifEmpty(Channel.empty()) )

      def exp_fastq = ['group','replicate','input_file','fasta','gtf','technology'] as Set

      def ch_original_only = ch_samplesheet.map { f, hdr -> file(f) }
      def ch_cleaned_only  = ch_cleaned

      def ch_safe_sheet = ch_original_only.combine(ch_cleaned_only).map { orig, cleaned ->
          def pick = file(cleaned)
          def ok = false
          try {
              def head = pick.readLines().first()
              def cols = head.split(',')*.trim()*.toLowerCase() as Set
              ok = exp_fastq.every { cols.contains(it) }
          } catch (Throwable ignore) { ok = false }
          ok ? pick : file(orig)
      }

      // Parseo de filas desde la hoja "segura"
      ch_samples = ch_safe_sheet
          .splitCsv(header: true, sep: ',', strip: true)
          .map { row ->
              try {
                  def found_cols = row.keySet()*.toString()
                  def lower = found_cols.collect { it.toLowerCase() } as Set
                  if( !exp_fastq.every{ lower.contains(it) } ) {
                      throw new Exception("CSV row missing one or more required columns: ${exp_fastq}. Found: ${found_cols}")
                  }

                  // Validar existencia de archivos
                  def input_file = file(row.input_file, checkIfExists: true)
                  def fasta_file = file(row.fasta,      checkIfExists: true)
                  def gtf_file   = file(row.gtf,        checkIfExists: true)

                  // Meta
                  def meta = [:]
                  meta.id             = "${row.group}-${row.replicate}".replaceAll("[^a-zA-Z0-9_.-]", "_")
                  meta.condition      = row.group.toString()
                  meta.replicate      = row.replicate.toString()
                  meta.fasta          = fasta_file
                  meta.gtf            = gtf_file
                  meta.technology     = row.technology.toString()
                  meta.is_transcripts = false
                  meta.annot          = row.group 

                  if( !(meta.technology in ['ONT','PacBio']) )
                      log.warn "Sample ${meta.id}: Technology '${meta.technology}' no reconocida explícitamente. Se asume compatible."

                  tuple(meta, input_file)
              }
              catch (Throwable e) {
                  error "Error processing samplesheet row for group '${row.group ?: 'MISSING'}' replicate '${row.replicate ?: 'MISSING'}': ${e.message}\n       Row data: ${row}"
              }
          }

    emit:
      samples                   = ch_samples
      validated_samplesheet     = ch_safe_sheet
      validated_samplesheet_raw = ch_cleaned
      //versions                  = ch_versions.unique()
}
