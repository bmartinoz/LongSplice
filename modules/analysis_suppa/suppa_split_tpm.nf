// modules/analysis_suppa/suppa_split_tpm.nf
nextflow.enable.dsl=2

process SUPPA_SPLIT_TPM {
  label 'process_low'
  container 'longsplice:latest'

  publishDir "${params.outdir ?: 'results'}/suppa/diffsplice",
             mode: 'copy', overwrite: true

  input:
  path tpm_tsv
  path design_tsv

  output:
  tuple path('cond1.tpm.tsv'), path('cond2.tpm.tsv'), emit: tpm_conds

  shell:
  '''
  set -euo pipefail

  IN_TPM="!{tpm_tsv}"
  IN_DESIGN="!{design_tsv}"

  if [[ -z "$IN_TPM" || ! -s "$IN_TPM" ]]; then
      IN_TPM="bambu.tpm.tsv"
  fi
  if [[ -z "$IN_DESIGN" || ! -s "$IN_DESIGN" ]]; then
      IN_DESIGN="design.tsv"
  fi

  [[ -s "$IN_TPM" ]]    || { echo "TPM no encontrado: $IN_TPM" >&2; exit 1; }
  [[ -s "$IN_DESIGN" ]] || { echo "design.tsv no encontrado: $IN_DESIGN" >&2; exit 1; }

  python3 - "$IN_TPM" "$IN_DESIGN" << 'PY'
import sys, os

TAB = chr(9)
NL  = chr(10)

if len(sys.argv) < 3:
    raise SystemExit("Uso interno: python - <tpm.tsv> <design.tsv>")

tpm    = sys.argv[1]
design = sys.argv[2]

if not os.path.exists(tpm):
    raise SystemExit(f"TPM no encontrado: {tpm}")
if not os.path.exists(design):
    raise SystemExit(f"design.tsv no encontrado: {design}")

# design.tsv: sample<TAB>condition
conds = {}
with open(design, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        line = line.rstrip()
        if not line:
            continue
        try:
            s, c = line.split(TAB)
        except ValueError:
            raise SystemExit("design.tsv: cada línea debe ser 'sample<TAB>condition'")
        conds.setdefault(c, []).append(s)

if len(conds) != 2:
    raise SystemExit("design.tsv debe contener exactamente 2 condiciones.")

# tpm: header = sólo muestras (TAB-sep); filas = transcript_id<TAB>valores
rows = []
with open(tpm, 'r', encoding='utf-8', errors='ignore') as f:
    samples = f.readline().rstrip().split(TAB)
    for line in f:
        parts = line.rstrip().split(TAB)
        tid, vals = parts[0], parts[1:]
        rows.append((tid, vals))

keys = list(conds.keys())
c1, c2 = keys[0], keys[1]

idx_map = { s:i for i,s in enumerate(samples) }
c1_idx = [idx_map[s] for s in conds[c1] if s in idx_map]
c2_idx = [idx_map[s] for s in conds[c2] if s in idx_map]

with open('cond1.tpm.tsv','w', encoding='utf-8') as o1, open('cond2.tpm.tsv','w', encoding='utf-8') as o2:
    # Agregar el encabezado con la columna de ID
    o1.write("event_id" + TAB + TAB.join(conds[c1]) + NL)
    o2.write("event_id" + TAB + TAB.join(conds[c2]) + NL)
    for tid, vals in rows:
        v1 = [vals[i] for i in c1_idx]
        v2 = [vals[i] for i in c2_idx]
        o1.write(tid + TAB + TAB.join(v1) + NL)
        o2.write(tid + TAB + TAB.join(v2) + NL)

PY
  '''
}