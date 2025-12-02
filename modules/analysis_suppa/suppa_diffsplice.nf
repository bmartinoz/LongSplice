// modules/analysis_suppa/suppa_diffsplice.nf
nextflow.enable.dsl=2

process SUPPA_DIFFSPLICE {
  label 'process_low'
  tag { event_type }
  container params.suppa_container

  publishDir "${params.outdir ?: 'results'}/suppa/diffsplice",
             mode: 'copy', overwrite: true

  input:
  tuple val(event_type), path(ioe), path(psi_file), path(cond1_tpm), path(cond2_tpm), path(design_tsv)

  output:
  path "diff_*.dpsi", emit: dpsi

  shell:
  '''
  set -euo pipefail

  EVT_TYPE="!{event_type}"
  IOE_FILE="!{ioe}"
  PSI_FILE="!{psi_file}"
  COND1_TPM="!{cond1_tpm}"
  COND2_TPM="!{cond2_tpm}"
  DESIGN_TSV="!{design_tsv}"



# Split PSI por condición
python3 - "$PSI_FILE" "$DESIGN_TSV" << 'PY'
import sys

TAB = '\t'
NL = chr(10)

# Archivos de entrada
psi_file, design_file = sys.argv[1], sys.argv[2]

# Leer archivo de diseño (sample ↔ condición)
conds = {}
with open(design_file, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        line = line.strip()
        if not line:
            continue
        sample, condition = line.split(TAB)
        conds.setdefault(condition, []).append(sample)

# Validar número de condiciones
if len(conds) != 2:
    raise SystemExit("Error: design.tsv debe contener exactamente 2 condiciones distintas.")

# Leer archivo PSI
with open(psi_file, 'r', encoding='utf-8', errors='ignore') as f:
    header = f.readline().rstrip().split(TAB)
    samples = header[1:]
    rows = [line.rstrip().split(TAB) for line in f if line.strip()]

# Índices de columnas de muestras
idx = {s: i for i, s in enumerate(samples)}

# Condiciones y sus muestras
c1, c2 = list(conds.keys())
c1_samples = [s for s in conds[c1] if s in idx]
c2_samples = [s for s in conds[c2] if s in idx]
c1_idx = [idx[s] for s in c1_samples]
c2_idx = [idx[s] for s in c2_samples]

# Verificar que ambas condiciones tengan muestras en el PSI
if not c1_idx or not c2_idx:
    raise SystemExit("Error: Ninguna muestra de alguna condición está presente en el PSI.")

# Escribir archivos separados por condición
with open('cond1.psi', 'w', encoding='utf-8') as o1, open('cond2.psi', 'w', encoding='utf-8') as o2:
    # Encabezados
    o1.write("event_id" + TAB + TAB.join(c1_samples) + NL)
    o2.write("event_id" + TAB + TAB.join(c2_samples) + NL)

    # Filas de eventos
    for parts in rows:
        event_id = parts[0]
        values = parts[1:]
        o1.write(event_id + TAB + TAB.join(values[i] for i in c1_idx) + NL)
        o2.write(event_id + TAB + TAB.join(values[i] for i in c2_idx) + NL)
PY



  OUTBASE="diff_${EVT_TYPE}"

  # Ejecutar SUPPA diffSplice
  suppa.py diffSplice \
    -i "$IOE_FILE" \
    -p cond1.psi cond2.psi \
    --tpm "$COND1_TPM" "$COND2_TPM" \
    --method !{ params.suppa_method ?: 'empirical' } \
    --tpm-threshold 10 \
    -nan 1 \
    -gc  \
    --area 1000 \
    -o "$OUTBASE"

  # Normalizar el nombre del .dpsi
  set +e
  found=""
  for cand in "${OUTBASE}.diffSplice.events.dpsi" "${OUTBASE}.diffSplice.dpsi" "${OUTBASE}.dpsi"; do
    if [[ -s "$cand" ]]; then
      [[ "$cand" != "${OUTBASE}.dpsi" ]] && mv -f "$cand" "${OUTBASE}.dpsi"
      found="${OUTBASE}.dpsi"
      break
    fi
  done
  if [[ -z "$found" ]]; then
    shopt -s nullglob
    dpsi=( *.dpsi )
    shopt -u nullglob
    if [[ ${#dpsi[@]} -eq 1 ]]; then
      [[ "${dpsi[0]}" != "${OUTBASE}.dpsi" ]] && mv -f "${dpsi[0]}" "${OUTBASE}.dpsi"
      found="${OUTBASE}.dpsi"
    fi
  fi
  set -e

  [[ -n "$found" && -s "$found" ]] || { echo "[SUPPA] ERROR: no se encontró ningún .dpsi tras diffSplice"; ls -la >&2 || true; exit 1; }

  if head -n1 "${OUTBASE}.dpsi" | grep -q '^cond1-cond2_dPSI[[:space:]]'; then
    { echo -e "event_id\t$(head -n1 "${OUTBASE}.dpsi")"; tail -n +2 "${OUTBASE}.dpsi"; } > "${OUTBASE}.dpsi.tmp"
    mv -f "${OUTBASE}.dpsi.tmp" "${OUTBASE}.dpsi"
  fi

  echo "[SUPPA] OK -> ${OUTBASE}.dpsi"
  '''
}
