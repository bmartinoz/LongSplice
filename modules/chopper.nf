// modules/chopper.nf
nextflow.enable.dsl = 2

process CHOPPER {
    tag "${meta.id}"
    label 'process_medium'
    container 'longsplice:latest'

    input:
    tuple val(meta), path(fastq_input)
    val chopper_args_in
    val output_prefix_in

    output:
    tuple val(meta), path("${output_prefix_in}.chopped.fastq.gz"), emit: chopped_reads
    path("${output_prefix_in}.chopper.log"),                       emit: chopper_logs
    path "versions.yml"                                         , emit: versions

    script:
    def out_fastq_gz = "${output_prefix_in}.chopped.fastq.gz"
    """
    #!/bin/bash
    set -euo pipefail

    LOG_FILE="${output_prefix_in}.chopper.log"
    TEMP_UNCOMPRESSED_FASTQ="${output_prefix_in}.temp_uncompressed.fastq"

    echo "=== Iniciando Chopper para ${output_prefix_in} ===" > \$LOG_FILE
    echo "Input FASTQ (variable fastq_input): ${fastq_input}" >> \$LOG_FILE
    echo "Tamaño del Input FASTQ: \$(ls -lh ${fastq_input} 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado_o_error')" >> \$LOG_FILE
    echo "Output FASTQ (variable out_fastq_gz): ${out_fastq_gz}" >> \$LOG_FILE
    echo "Threads: ${task.cpus}" >> \$LOG_FILE
    echo "Chopper args (variable chopper_args_in): ${chopper_args_in ?: "Ninguno"}" >> \$LOG_FILE
    echo "------------------------------------" >> \$LOG_FILE

    if [ ! -f "${fastq_input}" ]; then
        echo "ERROR CRÍTICO: Archivo de entrada ${fastq_input} no encontrado." >> \$LOG_FILE
        exit 1
    fi
    if [ ! -s "${fastq_input}" ]; then
        echo "ERROR CRÍTICO: Archivo de entrada ${fastq_input} ESTÁ VACÍO." >> \$LOG_FILE
        touch "${out_fastq_gz}" # Crear salida vacía para que NanoPlot falle si este es el caso
        # ... (código para versions.yml y exit 0 si quieres que NanoPlot muestre el error)
        # O simplemente exit 1 aquí para que CHOPPER falle directamente
        exit 1
    fi

    INPUT_TO_PIPE=""
    if [[ "${fastq_input}" == *.gz ]]; then
        echo "Descomprimiendo ${fastq_input} a \${TEMP_UNCOMPRESSED_FASTQ}" >> \$LOG_FILE
        zcat -f "${fastq_input}" > "\${TEMP_UNCOMPRESSED_FASTQ}"
        ZCAT_EXIT_CODE=\$?
        if [ \$ZCAT_EXIT_CODE -ne 0 ]; then
            echo "ERROR: zcat -f \"${fastq_input}\" falló con código \$ZCAT_EXIT_CODE." >> \$LOG_FILE
            ls -lh "${fastq_input}" >> \$LOG_FILE
            exit \$ZCAT_EXIT_CODE
        fi
        INPUT_TO_PIPE="\${TEMP_UNCOMPRESSED_FASTQ}"
        echo "Tamaño del archivo descomprimido temporal: \$(ls -lh \${INPUT_TO_PIPE} 2>/dev/null | awk '{print \$5}' || echo 'no_encontrado_o_error')" >> \$LOG_FILE
        echo "Primeras líneas del archivo descomprimido temporal:" >> \$LOG_FILE
        head -n 8 "\${INPUT_TO_PIPE}" >> \$LOG_FILE || echo "Error leyendo head de \${INPUT_TO_PIPE}" >> \$LOG_FILE
    else
        echo "La entrada ${fastq_input} no es .gz, usando directamente." >> \$LOG_FILE
        INPUT_TO_PIPE="${fastq_input}"
    fi

    if [ ! -s "\${INPUT_TO_PIPE}" ]; then
        echo "ERROR CRÍTICO: El archivo de entrada para chopper (\${INPUT_TO_PIPE}) está vacío después de la descompresión/asignación." >> \$LOG_FILE
        touch "${out_fastq_gz}" # Crear salida vacía
        exit 1 # Fallar aquí
    fi
    
    echo "------------------------------------" >> \$LOG_FILE
    echo "Comando Chopper a ejecutar:" >> \$LOG_FILE
    echo "cat \"\${INPUT_TO_PIPE}\" \\
        | chopper \\
            --threads ${task.cpus} \\
            ${chopper_args_in} \\
        2>> \$LOG_FILE \\
        | gzip -c > \"${out_fastq_gz}\"" >> \$LOG_FILE
    echo "------------------------------------" >> \$LOG_FILE

    cat "\${INPUT_TO_PIPE}" \\
        | chopper \\
            --threads ${task.cpus} \\
            ${chopper_args_in} \\
        2>> \$LOG_FILE \\
        | gzip -c > "${out_fastq_gz}"
    
    CHOPPER_PIPES_EXIT_CODES=( "\${PIPESTATUS[@]}" )
    CAT_EXIT_CODE=\${CHOPPER_PIPES_EXIT_CODES[0]}
    CHOPPER_TOOL_EXIT_CODE=\${CHOPPER_PIPES_EXIT_CODES[1]}
    GZIP_EXIT_CODE=\${CHOPPER_PIPES_EXIT_CODES[2]}

    rm -f "\${TEMP_UNCOMPRESSED_FASTQ}" # Limpiar archivo temporal

    if [ \$CAT_EXIT_CODE -ne 0 ]; then
        echo "ERROR: cat del archivo temporal falló con código \$CAT_EXIT_CODE" >> \$LOG_FILE
        exit \$CAT_EXIT_CODE
    fi
    if [ \$CHOPPER_TOOL_EXIT_CODE -ne 0 ]; then
        echo "ERROR: chopper (herramienta) falló con código \$CHOPPER_TOOL_EXIT_CODE" >> \$LOG_FILE
        exit \$CHOPPER_TOOL_EXIT_CODE
    fi
    if [ \$GZIP_EXIT_CODE -ne 0 ]; then
        echo "ERROR: gzip falló con código \$GZIP_EXIT_CODE" >> \$LOG_FILE
        exit \$GZIP_EXIT_CODE
    fi

    if [ ! -s "${out_fastq_gz}" ]; then
        echo "ERROR FINAL: El archivo de salida ${out_fastq_gz} está vacío después de ejecutar chopper y gzip." >> \$LOG_FILE
        exit 1
    fi

    echo "------------------------------------" >> \$LOG_FILE
    echo "Chopper completado exitosamente para ${output_prefix_in}." >> \$LOG_FILE

    CHOPPER_TOOL_VERSION=\$(chopper --version 2>&1 | sed -n 's/.*chopper //p' | sed -n 's/ .*//p' || \
                           chopper --version 2>&1 | grep -oP '(\\d+\\.)+\\d+' || \
                           echo "no_disponible")

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: \${CHOPPER_TOOL_VERSION}
    END_VERSIONS
    """

    stub:
    """
    echo "STUB: Creando archivos dummy para CHOPPER ${output_prefix_in}"
    touch "${output_prefix_in}.chopped.fastq.gz"
    touch "${output_prefix_in}.chopper.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        chopper: stub_version
    END_VERSIONS
    """
}