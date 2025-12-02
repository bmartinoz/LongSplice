// =============================
// MODULE: GET_CHROM_SIZES
// =============================

nextflow.enable.dsl = 2

process GET_CHROM_SIZES {

    tag "$id"
    label 'process_low'

    input:
    // Recibe id desde do_genome
    tuple val(id), path(fai)

    output:
    // Emite id y usa id en el nombre de archivo
    tuple val(id), path("${id}.sizes"), emit: sizes
    path "versions.yml"               , emit: versions

    script:
    """
    # Usa id para el nombre del archivo de salida
    cut -f1,2 ${fai} > ${id}.sizes

    # Check de archivo vacío/existencia
    if [ ! -s "${id}.sizes" ]; then
        echo "Error: Output file ${id}.sizes was not generated or is empty." >&2
        exit 1
    fi

    # Captura versión de 'cut' (coreutils) o samtools que generó el fai
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
        cut: \$(cut --version | head -n1)
    END_VERSIONS
    """

    stub:
    """
    touch ${id}.sizes
    echo "${task.process}:" > versions.yml
    echo "  samtools: stub" >> versions.yml
    echo "  cut: stub" >> versions.yml
    """
}
