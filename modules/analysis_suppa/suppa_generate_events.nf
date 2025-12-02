// modules/analysis_suppa/suppa_generate_events.nf
nextflow.enable.dsl=2

process SUPPA_GENERATE_EVENTS {
  label 'process_low'
  container params.suppa_container

  publishDir "${params.outdir ?: 'results'}/suppa/annotation",
             mode: 'copy', overwrite: true

  input:
  path gtf

  output:
  path "events_*.ioe", emit: ioe_all

  script:
    def BOUND_OPT = params.generateevents_boundary_mode ? "-b ${params.generateevents_boundary_mode}" : ""
    """
    set -euo pipefail
    # SS genera A5 y A3; FL genera AF y AL
    suppa.py generateEvents -i "${gtf}" -o events -f ioe -e SE SS MX RI FL ${BOUND_OPT}
    """
}