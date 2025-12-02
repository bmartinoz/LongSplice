// modules/analysis_suppa/suppa_psi_per_event.nf
nextflow.enable.dsl=2

process SUPPA_PSI_PER_EVENT {
  label 'process_low'
  container params.suppa_container

  publishDir "${params.outdir ?: 'results'}/suppa/psi",
             mode: 'copy', overwrite: true

  input:
  tuple val(type), path(ioe), path(tpm)

  output:
  path "${type}.psi", emit: psi

  shell:
  '''
  set -euo pipefail
  export LC_ALL=C
  export LANG=C.UTF-8
  export PYTHONUTF8=1

  type="!{type}"
  ioe="!{ioe}"
  tpm="!{tpm}"

  [[ -n "$type" ]] || { echo "[SUPPA] type vacío" >&2; exit 1; }
  [[ -s "$ioe"  ]] || { echo "[SUPPA] IOE no encontrado o vacío: $ioe" >&2; ls -la || true; exit 1; }
  [[ -s "$tpm"  ]] || { echo "[SUPPA] TPM no encontrado o vacío: $tpm" >&2; ls -la || true; exit 1; }

  echo "[SUPPA] psiPerEvent :: type=$type"
  echo "[SUPPA]   IOE: $ioe"
  echo "[SUPPA]   TPM: $tpm"

  suppa.py psiPerEvent -i "$ioe" -e "$tpm" -o "$type"

  if [[ -f "$type.ioe.psi" ]]; then
      mv "$type.ioe.psi" "$type.psi"
  elif [[ -f "$type.psi" ]]; then
      :
  else
      PSI_FOUND="$(ls -1 *.psi 2>/dev/null | wc -l || true)"
      if [[ "$PSI_FOUND" == "1" ]]; then
          PSI_NAME="$(ls -1 *.psi)"
          echo "[SUPPA] Renombrando $PSI_NAME -> $type.psi"
          mv "$PSI_NAME" "$type.psi"
      else
          echo "[SUPPA] ERROR: No se encontró salida .psi para type=$type" >&2
          ls -la || true
          exit 1
      fi
  fi

  if [[ -s "$type.psi" ]]; then
    first_col="$(head -n1 "$type.psi" | cut -f1)"
    if [[ "$first_col" != "event_id" ]]; then
      awk 'NR==1{print "event_id\t"$0; next}1' "$type.psi" > "$type.psi.tmp" && mv "$type.psi.tmp" "$type.psi"
    fi
  fi

  [[ -s "$type.psi" ]] || { echo "[SUPPA] $type.psi está vacío" >&2; exit 1; }
  echo "[SUPPA] OK -> $type.psi"
  '''
}