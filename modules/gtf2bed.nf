// =============================
// MODULE: GTF2BED (using BEDOPS convert2bed)
// =============================

nextflow.enable.dsl = 2

process GTF2BED {
    
    tag "$id"
    label 'process_medium'

    input:
    // Recibe id desde do_genome
    tuple val(id), path(gtf)

    output:
    // Emite id y usa id en el nombre de archivo
    tuple val(id), path("${id}.bed"), emit: gtf_bed
    path "versions.yml"             , emit: versions

    script:
    """
    # Usa id para el nombre del archivo de salida
    convert2bed --input=gtf < ${gtf} > ${id}.bed

    # Check de archivo vacío/existencia
    if [ ! -s "${id}.bed" ]; then
        echo "Error: Output file ${id}.bed was not generated or is empty." >&2
        exit 1
    fi

    # Captura la versión de convert2bed
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        convert2bed: \$(convert2bed --version 2>&1 | grep -oE 'version [0-9.]+' | awk '{print \$2}' || echo 'unknown')
        bedops: \$(bedops --version 2>&1 | head -n 1 || echo 'unknown')
    END_VERSIONS
    """

    stub:
    """
    touch ${id}.bed
    echo "${task.process}:" > versions.yml
    echo "  convert2bed: stub" >> versions.yml
    echo "  bedops: stub" >> versions.yml
    """
}
