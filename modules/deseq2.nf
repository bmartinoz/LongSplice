// modules/deseq2.nf
nextflow.enable.dsl=2

process DESEQ2 {
  tag "${samplesheet.baseName}"
  label 'process_medium'
  container 'longsplice:latest'

  publishDir "${params.outdir}/analysis_deseq2/${samplesheet.baseName}",
             mode: params.publish_dir_mode ?: 'copy',
             overwrite: true,
             saveAs: { it }

  /*
   Entradas:
     - counts_matrix (path)  -> stageAs counts_matrix.txt   (gene o transcript; formato Bambu-like)
     - samplesheet (path)    -> stageAs samplesheet_validated.csv
     - contrast_variable (val)  e.g. "group" o "condition"
     - reference_level  (val)   e.g. "WT"
     - target_level     (val)   e.g. "KO"
     - gtf (val) (opcional, para anotación)
  */
  input:
  path counts_matrix,     stageAs: "counts_matrix.txt"
  path samplesheet,       stageAs: "samplesheet_validated.csv"
  val  contrast_variable
  val  reference_level
  val  target_level
  val  gtf

  /*
   Salidas (lo que escribe run_deseq2.r):
     - deseq2_results.{txt,csv}
     - normalized_counts.csv
     - dds.rds  (+ opcional vsd.rds)
     - *.pdf (MA, PCA, Volcano)
     - analysis_summary.txt
     - versions.yml (lo genera el script; aquí lo reforzamos si faltara)
     - *.deseq2_script.log (log del wrapper bash)
  */
  output:
  path "deseq2_results.csv",      emit: results_csv,       optional: true
  path "deseq2_results.txt",      emit: deseq2_txt,        optional: true
  path "normalized_counts.csv",   emit: normalized_counts, optional: true
  path "*.pdf",                   emit: plots,             optional: true
  path "dds.rds",                 emit: rds,               optional: true
  path "vsd.rds",                 emit: vsd,               optional: true
  path "analysis_summary.txt",    emit: summary,           optional: true
  path "versions.yml",            emit: versions
  path "*.deseq2_script.log",     emit: log_file,          optional: true

 
  script:
  def qmethod  = params.quantification_method ?: 'bambu'
  def alpha    = params.alpha ?: 0.05
  def lfc_thr  = params.lfc_threshold ?: 1
  def shrink   = params.lfc_shrink_method ?: ''
  """
  #!/bin/bash
  set -euo pipefail

  LOGFILE="deseq2_results.deseq2_script.log"
  {
    echo "=== Wrapper DESeq2 (llamando a bin/run_deseq2.r) ==="
    echo "Fecha: \$(date)"
    echo "PWD: \$(pwd)"
    echo "Inputs:"
    echo "  counts_matrix.txt        : \$(ls -lh counts_matrix.txt 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado')"
    echo "  samplesheet_validated.csv: \$(ls -lh samplesheet_validated.csv 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado')"
    echo "  gtf                      : '${gtf}'"
    echo "Parámetros:"
    echo "  quant_method    : '${qmethod}'"
    echo "  contrast_var    : '${contrast_variable}'"
    echo "  ref_level       : '${reference_level}'"
    echo "  tgt_level       : '${target_level}'"
    echo "  alpha           : '${alpha}'"
    echo "  lfc_threshold   : '${lfc_thr}'"
    echo "  lfc_shrink      : '${shrink}'"
    echo "------------------------------------"
    echo "which run_deseq2.r -> \$(which run_deseq2.r || echo 'NO ENCONTRADO')"
  } > "\${LOGFILE}"

  # Construir comando (los opcionales pueden ir vacíos; el script los maneja)
  CMD=( run_deseq2.r
        "${qmethod}"
        "counts_matrix.txt"
        "samplesheet_validated.csv"
        "${contrast_variable}"
        "${reference_level}"
        "${target_level}"
        "${gtf}"
        "${alpha}"
        "${lfc_thr}"
        "${shrink}" )

  {
    echo "Comando:"
    printf ' %q' "\${CMD[@]}"; echo
    echo "------------------------------------"
  } >> "\${LOGFILE}"

  set +e
  "\${CMD[@]}" >> "\${LOGFILE}" 2>&1
  RC=\$?
  set -e
  if [ \$RC -ne 0 ]; then
    echo "------------------------------------" >> "\${LOGFILE}"
    echo "ERROR: run_deseq2.r falló con código \$RC." >> "\${LOGFILE}"
    tail -n 200 "\${LOGFILE}" >&2 || true
    exit \$RC
  fi

  # Reforzar/crear versions.yml si faltara
  if [ ! -s versions.yml ]; then
    {
      printf '%s\\n' "DESEQ2:"
      echo -n "  r-base: " 
      R --version 2>&1 | awk '/R version/ {print \$3}' || echo "NA"
      echo -n "  bioconductor-deseq2: "
      Rscript -e "cat(as.character(packageVersion('DESeq2')))" 2>/dev/null || echo "NA"
      echo
    } > versions.yml
  fi
  """
}
