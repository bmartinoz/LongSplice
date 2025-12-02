// modules/bambu.nf
process BAMBU {
    label 'process_medium'

    publishDir "${params.outdir ?: 'results'}/bambu",
        mode: 'copy',
        overwrite: true,
        enabled: { (task?.ext?.publishDir == null) ? true : task.ext.publishDir }

    // Contenedor: Singularity vs Docker
    container ( ((workflow.containerEngine ?: '') == 'singularity')
        ? 'https://depot.galaxyproject.org/singularity/mulled-v2-3979298465aa41b6990a37f126117735517f279f:af84e1451079953504788a38873b931f0e724970-0'
        : 'longsplice:latest' )

    input:
    tuple path(fasta), path(gtf)
    // lista de BAMs (staged por Nextflow)
    path bams
    val discovery
    val ndR
    val bambu_opts

    output:
    path "counts_gene.txt"         , emit: ch_gene_counts
    path "counts_transcript.txt"   , emit: ch_transcript_counts
    path "extended_annotations.gtf", emit: extended_gtf
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    def disco_flag = discovery ? true : false
    def ndr_flag   = ndR != null ? " ${ndR}" : ''
    def extra_opts = bambu_opts ?: ''


    """
    set -euo pipefail

    echo "[BAMBU] FASTA: ${fasta}" >&2
    echo "[BAMBU]   GTF: ${gtf}" >&2
    echo "[BAMBU] nBAMs: ${bams.size()}" >&2

    # Mostrar BAMs
    for f in ${ bams.collect { "'" + it + "'" }.join(' ') }; do
      echo "  - \$f" >&2
    done

    if [ ${bams.size()} -eq 0 ]; then
      echo "[BAMBU] ERROR: no llegaron BAMs al proceso" >&2
      exit 1
    fi

    run_bambu.r \\
      --tag=. \\
      --ncore=${task.cpus} \\
      --annotation="${gtf}" \\
      --fasta="${fasta}" \\
      --discovery="${disco_flag}" \\
      --ndr="${ndr_flag}" \\
      ${ bams.collect { "'" + it + "'" }.join(' ') }

    # ===== versions.yml =====
    Rscript -e 'cat(paste(R.version[["major"]], R.version[["minor"]], sep="."))' > RVER.txt 2>/dev/null || true
    Rscript -e 'suppressPackageStartupMessages(library(bambu)); cat(as.character(packageVersion("bambu")))' > BAMBU_VER.txt 2>/dev/null || true
    Rscript -e 'suppressPackageStartupMessages(library(BSgenome)); cat(as.character(packageVersion("BSgenome")))' > BSG_VER.txt 2>/dev/null || true

    {
      echo '"QUANTIFICATION_BAMBU:BAMBU":'
      echo -n '  r-base: ';                 if [ -s RVER.txt ]; then cat RVER.txt; else echo NA; fi
      echo -n '  bioconductor-bambu: ';     if [ -s BAMBU_VER.txt ]; then cat BAMBU_VER.txt; else echo NA; fi
      echo -n '  bioconductor-bsgenome: ';  if [ -s BSG_VER.txt ]; then cat BSG_VER.txt; else echo NA; fi
    } > versions.yml
    """
}
