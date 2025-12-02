// modules/analysis_suppa/suppa_generate_ioi.nf
nextflow.enable.dsl=2

process SUPPA_GENERATE_IOI {
  label 'process_low'
  container params.suppa_container

  publishDir "${params.outdir ?: 'results'}/suppa/annotation",
             mode: 'copy', overwrite: true

  input:
  path gtf

  output:
  path "annotation.ioi", emit: ioi
  path "versions.yml",  emit: versions

  script:
  """
  set -euo pipefail

  GTF_IN="\${gtf-}"
  if [[ -z "\$GTF_IN" || ! -s "\$GTF_IN" ]]; then
      GTF_IN="\$(ls -1 *.gtf 2>/dev/null | head -n1 || true)"
  fi
  if [[ -z "\$GTF_IN" || ! -s "\$GTF_IN" ]]; then
      echo '[SUPPA] ERROR: no se encontró GTF de entrada en SUPPA_GENERATE_IOI.' >&2
      ls -la || true
      exit 1
  fi

  suppa.py generateEvents -i "\$GTF_IN" -o annotation -f ioi

  python3 - <<'PY' > SUPPA_VER.txt 2>/dev/null || true
try:
    import importlib.metadata as im
    print(im.version('SUPPA2'))
except Exception:
    print('NA')
PY

  {
    echo 'ANALYSIS_SUPPA:'
    printf '  suppa2: '
    cat SUPPA_VER.txt
  } > versions.yml || true
  """
}