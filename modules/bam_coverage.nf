// =======================================
// MODULE: BAM_COVERAGE
// =======================================

nextflow.enable.dsl = 2

process BAM_COVERAGE {
    tag   "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir ?: 'results'}/coverage", mode: 'copy',
        saveAs: { filename -> filename }

    conda (params.enable_conda ? "bioconda::deeptools=3.5.5" : null)

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                 'docker://longsplice:latest' : 'longsplice:latest' }"

    input:
    tuple val(meta), path(bam)

    output:
    path "*.bw",         emit: bigwig
    path "versions.yml", emit: versions

    script:
    def norm     = params.bigwig_norm    ?: 'CPM'
    def binsize  = params.bigwig_binsize ?: 50
    def nthreads = task.cpus ?: 4

    """
    set -euo pipefail

    # ---------- Asegurar que exista índice .bai ----------
    if [ ! -e "${bam}.bai" ]; then
        samtools index ${bam}
    fi

    # ---------- Generar bigWig ----------
    bamCoverage \\
        -b ${bam} \\
        -o ${meta.id}.bw \\
        --outFileFormat bigwig \\
        --normalizeUsing ${norm} \\
        --binSize ${binsize} \\
        -p ${nthreads}

    # ---------- Registrar versiones ----------
    BAMCOV_VERSION=\$(bamCoverage --version 2>&1 | awk '{print \$NF}' || echo "unknown")
    DEEPTOOLS_VERSION=\$(python - <<'PY'
    import pkg_resources, sys
    try:
        print(pkg_resources.get_distribution('deeptools').version)
    except Exception:
        sys.stdout.write('unknown')
    PY
    )

    cat <<- 'END_VERSIONS' > versions.yml
    \"${task.process}\":
        bamCoverage:  \${BAMCOV_VERSION}
        deeptools:    \${DEEPTOOLS_VERSION}
    END_VERSIONS
    """
}
