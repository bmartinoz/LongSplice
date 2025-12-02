// =======================================
// MODULE: SAMTOOLS_FAIDX
// =======================================

nextflow.enable.dsl = 2

process SAMTOOLS_FAIDX {

    tag "$id"
    label 'process_low'

    input:
    // Recibe id desde do_genome
    tuple val(id), path(fasta)

    output:
    tuple val(id), path("${id}.fa.fai"), emit: fai
    path "versions.yml"                 , emit: versions


    script:
    """
    samtools faidx ${fasta}

    # Check que el índice se creó
    if [ ! -f "${fasta}.fai" ]; then
        echo "Error: FASTA index file ${fasta}.fai was not created." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
     """
     touch ${fasta}.fai
     echo "${task.process}:" > versions.yml
     echo "  samtools: stub" >> versions.yml
     """
}
