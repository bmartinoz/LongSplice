// subworkflows/make_contrastsheet.nf  (corregido)

nextflow.enable.dsl = 2

//
// (1) Proceso: exponer el samplesheet con nombre estándar sin fallar si ya está con ese nombre
//
process MAKE_VALIDATED_SAMPLESHEET_PROC {
    tag { samplesheet_in.baseName }
    label 'process_low'
    // container 'longsplice:latest'

    input:
      path samplesheet_in, stageAs: "samplesheet_validated.csv"

    output:
      path "samplesheet_validated.csv", emit: validated

    script:
    """
    set -euo pipefail
    # Si el archivo ya está en el workdir con el mismo nombre (por el stageAs), no hacemos nada.
    if [ -e "samplesheet_validated.csv" ]; then
      :
    else
      # En caso raro de que no esté, intentamos copiar o enlazar.
      cp -f "\${samplesheet_in}" "samplesheet_validated.csv" 2>/dev/null || ln -sf "\${samplesheet_in}" "samplesheet_validated.csv"
    fi
    """
}

//
// (2) Proceso: generar contrasts.csv de forma determinística (sin R)
//
process MAKE_CONTRASTSHEET_PROC {
    tag { samplesheet_validated.baseName }
    label 'process_low'
    // container 'longsplice:latest'   // (opcional; usa el de tu config)

    input:
      // Lo leemos con el mismo nombre estándar para simplificar el script
      path samplesheet_validated, stageAs: "samplesheet_validated.csv"

    output:
      path "contrasts.csv", emit: contrasts_csv

    script:
    """
    set -euo pipefail
    python3 - <<'PY'
import csv, itertools, sys, pathlib
p = pathlib.Path("samplesheet_validated.csv")

# Leer CSV (tolerante a encoding)
with p.open(encoding='utf-8', errors='ignore', newline='') as fh:
    rows = list(csv.DictReader(fh))

if not rows:
    print("ERROR: samplesheet vacío", file=sys.stderr)
    sys.exit(1)

hdr = [h.strip().lower() for h in rows[0].keys()]

def pick_col(candidates):
    # match exact primero
    for c in candidates:
        if c in hdr:
            return c
    # si no exact, intenta contains
    for c in candidates:
        for h in hdr:
            if c in h:
                return h
    return None

# Buscar columna de condición
cond_col = pick_col(['group','condition','phenotype','status','class','label'])
if not cond_col:
    print("ERROR: no se encontró columna de condición (group/condition/phenotype/etc.)", file=sys.stderr)
    sys.exit(1)

# Normalizar y recolectar condiciones no vacías
conds = sorted({ (r.get(cond_col) or '').strip() for r in rows if (r.get(cond_col) or '').strip() })
if len(conds) < 2:
    print("ERROR: se requieren ≥2 niveles de condición para contrasts", file=sys.stderr)
    sys.exit(1)

# Todas las combinaciones pares (A vs B, A<B para determinismo)
pairs = list(itertools.combinations(conds, 2))

with open('contrasts.csv','w',encoding='utf-8',newline='') as out:
    w = csv.writer(out)
    w.writerow(['treatment','control'])
    w.writerows(pairs)
PY
    """
}

//
// (3) Subworkflow: encadenar ambos procesos
//
workflow MAKE_CONTRASTSHEET {
    take:
      ch_samplesheet_in  // conectar desde main.nf (p. ej., DO_INPUT.out.validated_samplesheet)

    main:
      MAKE_VALIDATED_SAMPLESHEET_PROC(ch_samplesheet_in)
      def ch_validated = MAKE_VALIDATED_SAMPLESHEET_PROC.out.validated

      MAKE_CONTRASTSHEET_PROC(ch_validated)
      def ch_contrasts = MAKE_CONTRASTSHEET_PROC.out.contrasts_csv

    emit:
      samplesheet_validated = ch_validated
      contrasts             = ch_contrasts
}
