// =======================================
// MODULE: RESTRANDER
// =======================================

nextflow.enable.dsl = 2

process RESTRANDER {
  tag "$meta.id"
  //label 'process_low'

  publishDir "${params.outdir}/restrander/${meta.id}",
             mode: params.publish_dir_mode ?: 'copy',
             pattern: "*"

  input:
  tuple val(meta), path(reads)  //tuple
  val   cfg_path
  val   prefix

  output:
  tuple val(meta), path("${prefix}.restranded.fastq.gz"), emit: restranded_reads
  path  "${prefix}.restrander.log", emit: restrander_log
  path  "${prefix}.stats.json",     emit: restrander_stats
  path  "versions.yml",             emit: versions

  when:
  !params.skip_restrander

  shell:
  '''
  set -euo pipefail

  in_fastq="!{reads}"
  cfg="!{cfg_path}"
  out_fastq="!{prefix}.restranded.fastq.gz"
  log_file="!{prefix}.restrander.log"
  stats_file="!{prefix}.stats.json"

  {
    echo "[\$(date -u +'%F %T')] start"
    echo "Input : ${in_fastq}"
    echo "Config: ${cfg}"
    echo "Output: ${out_fastq}"
    restrander --version || true

    [[ -s "${in_fastq}" ]] || { echo "ERROR: input FASTQ not found or empty: ${in_fastq}"; exit 2; }
    [[ -f "${cfg}"     ]]  || { echo "ERROR: config JSON not found: ${cfg}"; exit 3; }

    if ! restrander "${in_fastq}" "${out_fastq}" "${cfg}" > "${stats_file}"; then
      echo "ERROR: restrander exited non-zero" >&2
      exit 1
    fi

    [[ -s "${out_fastq}" ]] || { echo "ERROR: output FASTQ missing: ${out_fastq}"; exit 4; }
    [[ -s "${stats_file}" ]] || { echo "WARN: stats JSON empty"; }

    echo "[\$(date -u +'%F %T')] done"
  } 2>&1 | tee "${log_file}"

  {
    echo "RESTRANDER:"
    (restrander --version 2>/dev/null | awk 'NF{print "  " $0}') || echo "  version: unknown"
  } > versions.yml
  '''
}
