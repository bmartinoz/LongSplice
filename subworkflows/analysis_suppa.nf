// =======================================
// SUBWORKFLOW: ANALYSIS_SUPPA
// =======================================

nextflow.enable.dsl = 2

// ===== Módulos SUPPA =====
include { SUPPA_GENERATE_IOI      } from '../modules/analysis_suppa/suppa_generate_ioi.nf'
include { SUPPA_GENERATE_EVENTS   } from '../modules/analysis_suppa/suppa_generate_events.nf'
include { SUPPA_COUNTS_TO_TPM     } from '../modules/analysis_suppa/suppa_counts_to_tpm.nf'
include { SUPPA_PSI_PER_EVENT     } from '../modules/analysis_suppa/suppa_psi_per_event.nf'
include { SUPPA_SPLIT_TPM         } from '../modules/analysis_suppa/suppa_split_tpm.nf'
include { SUPPA_DIFFSPLICE        } from '../modules/analysis_suppa/suppa_diffsplice.nf'

// -------------------- SUBWORKFLOW --------------------
workflow ANALYSIS_SUPPA {
  take:
  ch_gtf
  ch_counts_tx
  ch_design

  main:
  //ch_versions = Channel.empty()

  // (1) IOI (desde GTF) + IOE (eventos)
  def gen = SUPPA_GENERATE_IOI(ch_gtf)
  //ch_versions = ch_versions.mix( gen.versions.ifEmpty(Channel.empty()) )
  def ioe = SUPPA_GENERATE_EVENTS(ch_gtf)

  // (2) TPM desde counts + GTF
  def tpm = SUPPA_COUNTS_TO_TPM(ch_counts_tx, ch_gtf)

  // (3) Etiquetar IOE por tipo
  def ch_ioe_typed = ioe.ioe_all
    .flatten()
    .map { Path f ->
      def bn = f.baseName
      def m  = (bn =~ /^events_(SE|MX|RI|FL|A5|A3|AF|AL)(?:_strict)?$/)
      def t  = m ? m[0][1] : 'EV'
      tuple(t, f)  // (type, ioe_file)
    }

  // (3b) Broadcast del TPM a cada IOE
  def ch_psi_inputs = ch_ioe_typed
    .combine(tpm.tpm_tsv)
    .map { type, ioe_file, tpmf -> tuple(type, ioe_file, tpmf) }

  // PSI por evento
  def psi_out = SUPPA_PSI_PER_EVENT(ch_psi_inputs)

  // (4) diffSplice opcional
  if (params.run_suppa_diff && ch_design) {

    // PSI etiquetado por tipo
    def ch_psi_typed = psi_out.psi.map { Path f ->
      def m = (f.baseName =~ /^(SE|MX|RI|FL|A5|A3|AF|AL)$/)
      def type = m ? m[0][1] : 'EV'
      tuple(type, f)
    }

    // Empareja
    def ch_evt = ch_psi_typed
      .join(ch_ioe_typed)
      .map { type, psi_file, ioe_file -> tuple(type, ioe_file, psi_file) }

    // TPM por condición
    def tpm_conds = SUPPA_SPLIT_TPM(tpm.tpm_tsv, ch_design)
    def ch_diff_inputs = ch_evt
      .combine(tpm_conds.tpm_conds)
      .combine(ch_design)

    SUPPA_DIFFSPLICE(ch_diff_inputs)
  }

  emit:
  ioe_files  = ioe.ioe_all
  tpm_matrix = tpm.tpm_tsv
  psi_files  = psi_out.psi
  //versions   = ch_versions.unique()
}
