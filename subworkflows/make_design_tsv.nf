// =======================================
// SUBWORKFLOW: MAKE_DESIGN_TSV
// =======================================

nextflow.enable.dsl = 2

process BUILD_DESIGN_TSV {
  label 'process_low'
  container params.suppa_container

  publishDir "${params.outdir ?: 'results'}/suppa",
             mode: 'copy', overwrite: true

  input:
  path validated_csv

  output:
  path "design.tsv", emit: design

  script:
  """
  set -euo pipefail

  # Pasamos valores al entorno del script Python
  export IN_CSV="${validated_csv}"
  export DESIGN_TEMPLATE='${params.design_sample_template ?: '{group}-{replicate}.sorted'}'

  python - << 'PY'
import csv, os, sys

in_csv  = os.environ['IN_CSV']
tmpl    = os.environ.get('DESIGN_TEMPLATE', '{group}-{replicate}.sorted')

rows = []
with open(in_csv, newline='') as f:
    reader = csv.DictReader(f)
    for r in reader:
        # normaliza claves/valores
        r = { (k or '').strip(): (v or '').strip() for k,v in r.items() }
        if not r.get('group') or not r.get('replicate'):
            continue
        sample = tmpl.format(
            group=r['group'],
            replicate=r['replicate'],
            input_file=r.get('input_file',''),
            technology=r.get('technology','')
        )
        cond = r['group']
        rows.append((sample, cond))

if not rows:
    sys.exit("No se pudieron generar filas para design.tsv (revisa 'group' y 'replicate').")

# Quitar duplicados preservando orden
seen = set()
uniq = []
for s,c in rows:
    if (s,c) not in seen:
        uniq.append((s,c))
        seen.add((s,c))

with open('design.tsv','w') as o:
    for s,c in uniq:
        o.write(f"{s}\\t{c}\\n")
PY
  """
}

workflow MAKE_DESIGN_TSV {
  take:
  ch_validated_sheet   // path al samplesheet_validated.csv

  main:
  def outp = BUILD_DESIGN_TSV(ch_validated_sheet)

  emit:
  design = outp.design
}
