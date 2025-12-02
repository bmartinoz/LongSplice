// =======================================
// MODULE: NANOPLOT
// =======================================

nextflow.enable.dsl = 2

process NANOPLOT {
    tag "$meta.id (${report_subdir_name})"
    publishDir "${params.outdir}/nanoplot/${meta.id}",
               mode: 'copy',
               pattern: "${report_subdir_name}"

    label 'process_low'

    input:
    tuple val(meta), path(reads)
    val report_subdir_name

    output:
    tuple val(meta), path(report_subdir_name), emit: report_dir
    path "versions.yml"                       , emit: versions

    script:
    def default_nanoplot_opts = "--minlength 0 --minqual 0"
    def nanoplot_args = task.ext.args ?: (params.nanoplot_opts ?: default_nanoplot_opts)
    def threads = 1  // evita multiprocessing problemático (Py 3.12)

    """
    set -euo pipefail
    mkdir -p "${report_subdir_name}" .tmp .mpl

    # Backend y rutas seguras para matplotlib/TMP
    export MPLBACKEND=Agg
    export TMPDIR="\$PWD/.tmp"
    export MPLCONFIGDIR="\$PWD/.mpl"

    # Chequeo rápido de integridad (sin regex; pattern matching)
    if [[ "${reads}" == *.gz ]]; then
      gzip -t "${reads}" || { echo "FASTQ corrupto: ${reads}" >&2; exit 2; }
    else
      [[ -s "${reads}" ]] || { echo "FASTQ vacío: ${reads}" >&2; exit 2; }
    fi

    run_nanoplot() {
      NanoPlot \\
        --fastq "${reads}" \\
        --loglength \\
        --N50 \\
        --outdir "${report_subdir_name}" \\
        --threads ${threads} \\
        ${nanoplot_args}
    }

    # Intento 1
    if ! run_nanoplot ; then
      echo "[WARN] NanoPlot falló; reintento con --downsample 200000…" >&2
      rm -rf "${report_subdir_name}" && mkdir -p "${report_subdir_name}"
      run_nanoplot --downsample 200000
    fi

    # Validación (si no hay archivos, no abortar: deja un warning y continúa)
    if [ ! -d "${report_subdir_name}" ] || [ -z "\$(ls -A "${report_subdir_name}")" ]; then
      echo "[WARN] NanoPlot terminó pero no generó archivos (posible problema de export estático / Kaleido)." >&2
      echo "No static plots due to missing chromium deps or headless export problem." > "${report_subdir_name}/EXPORT_WARNING.txt"
    fi

    cat > versions.yml <<EOF
"${task.process} (${report_subdir_name})":
    NanoPlot: \$(NanoPlot --version 2>/dev/null | sed -E 's/.* ([0-9.]+).*/\\1/')
EOF
    """

    stub:
    """
    mkdir -p "${report_subdir_name}"
    echo "EXPORT_WARNING: stub run (no plots)" > "${report_subdir_name}/EXPORT_WARNING.txt"
    cat > versions.yml <<EOF
"${task.process} (${report_subdir_name})":
  NanoPlot: stub
EOF
    """
}
