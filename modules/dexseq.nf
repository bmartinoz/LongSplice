// modules/dexseq.nf
nextflow.enable.dsl=2

process DEXSEQ {
  tag "${samplesheet.baseName}"
  label 'process_medium'
  container 'longsplice:latest'

  publishDir "${params.outdir}/analysis_dtu/dexseq/${samplesheet.baseName}",
             mode: params.publish_dir_mode ?: 'copy',
             overwrite: true,
             saveAs: { it }

  input:
  path counts_transcript, stageAs: "counts_transcript.txt"
  path samplesheet,       stageAs: "samplesheet_validated.csv"
  path extended_gtf,      stageAs: "extended_annotations.gtf"

  output:
  path "${task.ext.prefix ?: 'dexseq_results'}.txt", emit: results_txt, optional: true
  path "versions.yml",                               emit: versions
  path "*.dexseq_script.log",                        emit: log_file, optional: true
  path "dexseq_dropped_samples.txt",                 emit: dropped, optional: true
  path "dexseq_results.str.txt",                     emit: struct_dump, optional: true

  script:
  def prefix     = task.ext.prefix ?: "dexseq_results"
  def args       = task.ext.args ?: (params.dexseq_opts ?: '')
  def dexseq_r   = "${projectDir}/bin/run_dexseq.r"   // usa SIEMPRE tu script local

  """
  #!/bin/bash
  set -euo pipefail

  BASH_SCRIPT_LOG="${prefix}.dexseq_script.log"
  {
    echo "=== Wrapper DEXSeq (llamando a script local stageado) ==="
    echo "Fecha: \$(date)"
    echo "PWD: \$(pwd)"
    echo "Inputs:"
    echo "  counts_transcript.txt     : \$(ls -lh counts_transcript.txt 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado')"
    echo "  samplesheet_validated.csv : \$(ls -lh samplesheet_validated.csv 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado')"
    echo "  extended_annotations.gtf  : \$(ls -lh extended_annotations.gtf 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado')"
    echo "------------------------------------"
    echo "projectDir run_dexseq.r     : ${dexseq_r}"
    echo "which run_dexseq.r (contenedor) -> \$(which run_dexseq.r || echo 'NO ENCONTRADO')"
  } > "\${BASH_SCRIPT_LOG}"

  # Stagear nuestro script (aseguramos que es el que se ejecuta)
  cp "${dexseq_r}" ./run_dexseq.r
  chmod +x ./run_dexseq.r

  # Comando principal
  CMD=( ./run_dexseq.r
        --method "bambu"
        --counts "counts_transcript.txt"
        --samplesheet "samplesheet_validated.csv"
        --annotation "extended_annotations.gtf"
        --outfile "${prefix}.txt"
        --ncore "\${NXF_TASK_CPUS:-4}" )

  # Args extra (opcionales)
  if [ -n "${args}" ]; then
    CMD+=( ${args} )
  fi

  {
    echo "Comando:"
    printf ' %q' "\${CMD[@]}"; echo
    echo "------------------------------------"
  } >> "\${BASH_SCRIPT_LOG}"

  set +e
  "\${CMD[@]}" >> "\${BASH_SCRIPT_LOG}" 2>&1
  RC=\$?
  set -e
  if [ \$RC -ne 0 ]; then
    echo "------------------------------------" >> "\${BASH_SCRIPT_LOG}"
    echo "ERROR: run_dexseq.r falló con código \$RC." >> "\${BASH_SCRIPT_LOG}"
    tail -n 200 "\${BASH_SCRIPT_LOG}" >&2 || true
    exit \$RC
  fi

  # versions.yml
  {
    printf '%s\\n' "\"${task.process}\":" > versions.yml
    echo -n "  r-base: " >> versions.yml
    R --version 2>&1 | awk '/R version/ {print \$3}' >> versions.yml || echo "NA" >> versions.yml
    echo -n "  bioconductor-dexseq: " >> versions.yml
    Rscript -e "cat(as.character(packageVersion('DEXSeq')))" >> versions.yml || echo "NA" >> versions.yml
    echo -n "  bioconductor-drimseq: " >> versions.yml
    Rscript -e "cat(as.character(packageVersion('DRIMSeq')))" >> versions.yml || echo "NA" >> versions.yml
    echo -n "  bioconductor-stager: " >> versions.yml
    Rscript -e "cat(as.character(packageVersion('stageR')))" >> versions.yml || echo "NA" >> versions.yml
    echo >> versions.yml
  }

  # archivo opcional de "dropped"
  if [ ! -s "dexseq_dropped_samples.txt" ]; then
    echo "No se generó listado de muestras descartadas (no requerido)." > dexseq_dropped_samples.txt || true
  fi
  """
}
