// =======================================
// MODULE: MINIMAP2_INDEX
// =======================================

nextflow.enable.dsl = 2

process MINIMAP2_INDEX {
    tag "${fasta.getSimpleName()}"
    label 'process_medium'

    input:
    path fasta

    output:
    tuple path(fasta), path("${fasta.simpleName}.mmi"), emit: mmi
    path "versions.yml", emit: versions

    script:
    def base = fasta.simpleName

    """
    minimap2 -d ${base}.mmi ${fasta} -t ${task.cpus}

    if [ ! -s "${base}.mmi" ]; then
        echo "Error: Minimap2 index file ${base}.mmi was not generated or is empty." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
    END_VERSIONS
    """

    stub:
    """
    touch ${fasta.simpleName}.mmi
    echo "${task.process}:" > versions.yml
    echo "  minimap2: stub" >> versions.yml
    """
}
