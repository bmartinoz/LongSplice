process BAM_MERGE {
    label 'process_medium'

    publishDir "${params.outdir ?: 'results'}/merged_bam", mode: 'copy',
        saveAs: { filename -> filename }

    conda (params.enable_conda ? "bioconda::deeptools=3.5.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                 'docker://longsplice:latest' : 'longsplice:latest' }"

    input:
    path bam_files

    output:
    path "merged.genomealigned.bam", emit: merged_bam
    path "versions.yml", emit: versions

    script:
    def SORT_MEM = (params.samtools_sort_mem ?: '8G')

    """
    set -euo pipefail


    """"

    }