// modules/pychopper.nf
nextflow.enable.dsl = 2

process PYCHOPPER {
    tag "${meta.id}"
    label 'process_high'
    
    container 'longsplice:latest'

    input:
    tuple val(meta), path(fastq_input)
    val pychopper_args_in
    val output_prefix_in

    output:
    tuple val(meta), path("${output_prefix_in}.chopped.fastq.gz"), emit: fastq
    path "${output_prefix_in}.chopper_stats.tsv",                  emit: stats, optional: true
    path "${output_prefix_in}.chopper_report.pdf",                emit: report, optional: true
    path "${output_prefix_in}.unclassified.fastq.gz",            emit: unclassified, optional: true
    path "${output_prefix_in}.rescued.fastq.gz",                  emit: rescued, optional: true
    path "${output_prefix_in}.pychopper_script.log",              emit: log_file, optional: true
    path "versions.yml",                                         emit: versions

    script:
    def user_kit = params.cdna_kit 

    def pychopper_cmd_parts = ['pychopper']
    pychopper_cmd_parts << "-S ${output_prefix_in}.chopper_stats.tsv"
    
    if (user_kit && user_kit.toLowerCase() != "auto" && user_kit.toLowerCase() != "none" && user_kit.trim() != "") {
        pychopper_cmd_parts << "-k ${user_kit}"
    }
    
    pychopper_cmd_parts << "-t ${task.cpus}"
    pychopper_cmd_parts << "-r ${output_prefix_in}.chopper_report.pdf"
    pychopper_cmd_parts << "-u ${output_prefix_in}.unclassified.fastq"
    pychopper_cmd_parts << "-w ${output_prefix_in}.rescued.fastq"
    
    if (pychopper_args_in && pychopper_args_in.trim() != "") {
        pychopper_cmd_parts << pychopper_args_in
    }
    
    pychopper_cmd_parts << "\"${fastq_input}\""
    pychopper_cmd_parts << "\"${output_prefix_in}.chopped.fastq\""
    
    def pychopper_command = pychopper_cmd_parts.join(' \\\n        ')

    """
    #!/bin/bash
    set -e
    set -o pipefail

    PYCHOPPER_PROCESS_LOG="${output_prefix_in}.pychopper_script.log"

    echo "=== Iniciando Pychopper para ${output_prefix_in} ===" | tee \$PYCHOPPER_PROCESS_LOG
    echo "Input FASTQ: ${fastq_input}" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Working directory: \$(pwd)" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Pipeline User Kit Param (params.cdna_kit): ${user_kit ?: "No especificado (o 'auto'/'none')"}" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Extra Args (input pychopper_args_in): ${pychopper_args_in ?: "Ninguno"}" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Threads para Pychopper: ${task.cpus}" | tee -a \$PYCHOPPER_PROCESS_LOG
    
    # Verificar que el archivo de entrada existe y es accesible
    echo "------------------------------------" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Verificando archivo de entrada..." | tee -a \$PYCHOPPER_PROCESS_LOG
    if [ ! -f "${fastq_input}" ]; then
        echo "ERROR: Archivo de entrada no encontrado: ${fastq_input}" | tee -a \$PYCHOPPER_PROCESS_LOG
        echo "Contenido del directorio actual:" | tee -a \$PYCHOPPER_PROCESS_LOG
        ls -lha | tee -a \$PYCHOPPER_PROCESS_LOG
        exit 1
    fi
    
    echo "Archivo encontrado:" | tee -a \$PYCHOPPER_PROCESS_LOG
    ls -lh "${fastq_input}" | tee -a \$PYCHOPPER_PROCESS_LOG
    
    # Verificar que no sea un enlace simbólico roto
    if [ -L "${fastq_input}" ]; then
        echo "Es un enlace simbólico a: \$(readlink -f "${fastq_input}" || echo "ROTO")" | tee -a \$PYCHOPPER_PROCESS_LOG
        REAL_FILE=\$(readlink -f "${fastq_input}" 2>/dev/null || echo "")
        if [ -z "\$REAL_FILE" ] || [ ! -f "\$REAL_FILE" ]; then
            echo "ERROR: Enlace simbólico roto" | tee -a \$PYCHOPPER_PROCESS_LOG
            exit 1
        fi
    fi
    
    echo "Tamaño del archivo: \$(stat -c%s "${fastq_input}" 2>/dev/null || stat -f%z "${fastq_input}" 2>/dev/null || echo "unknown") bytes" | tee -a \$PYCHOPPER_PROCESS_LOG
    
    echo "------------------------------------" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Comando Pychopper a ejecutar:" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "${pychopper_command}" | tee -a \$PYCHOPPER_PROCESS_LOG 
    echo "------------------------------------" | tee -a \$PYCHOPPER_PROCESS_LOG

    # Ejecutar pychopper con salida duplicada (log + stderr)
    ${pychopper_command} 2>&1 | tee -a \$PYCHOPPER_PROCESS_LOG
    
    PYCHOPPER_EXIT_CODE=\${PIPESTATUS[0]}

    if [ \$PYCHOPPER_EXIT_CODE -ne 0 ]; then
        echo "------------------------------------" | tee -a \$PYCHOPPER_PROCESS_LOG
        echo "ERROR: pychopper falló con código de salida \$PYCHOPPER_EXIT_CODE para ${output_prefix_in}" | tee -a \$PYCHOPPER_PROCESS_LOG
        echo "Archivos generados:" | tee -a \$PYCHOPPER_PROCESS_LOG
        ls -lh | tee -a \$PYCHOPPER_PROCESS_LOG
        exit \$PYCHOPPER_EXIT_CODE
    fi

    echo "------------------------------------" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Pychopper completado exitosamente para ${output_prefix_in}." | tee -a \$PYCHOPPER_PROCESS_LOG

    # Comprimir archivo principal
    if [ -f "${output_prefix_in}.chopped.fastq" ]; then
        echo "Comprimiendo ${output_prefix_in}.chopped.fastq..." | tee -a \$PYCHOPPER_PROCESS_LOG
        gzip -f "${output_prefix_in}.chopped.fastq"
        echo "Archivo comprimido exitosamente: \$(ls -lh ${output_prefix_in}.chopped.fastq.gz)" | tee -a \$PYCHOPPER_PROCESS_LOG
    else
        echo "ADVERTENCIA: Archivo principal ${output_prefix_in}.chopped.fastq no encontrado" | tee -a \$PYCHOPPER_PROCESS_LOG
        echo "Archivos presentes en el directorio:" | tee -a \$PYCHOPPER_PROCESS_LOG
        ls -lh | tee -a \$PYCHOPPER_PROCESS_LOG
        echo "Creando .gz vacío como fallback" | tee -a \$PYCHOPPER_PROCESS_LOG
        touch "${output_prefix_in}.chopped.fastq.gz"
    fi

    # Comprimir archivos opcionales
    if [ -f "${output_prefix_in}.unclassified.fastq" ]; then
        echo "Comprimiendo ${output_prefix_in}.unclassified.fastq..." | tee -a \$PYCHOPPER_PROCESS_LOG
        gzip -f "${output_prefix_in}.unclassified.fastq"
    fi
    
    if [ -f "${output_prefix_in}.rescued.fastq" ]; then
        echo "Comprimiendo ${output_prefix_in}.rescued.fastq..." | tee -a \$PYCHOPPER_PROCESS_LOG
        gzip -f "${output_prefix_in}.rescued.fastq"
    fi

    # Capturar versión de pychopper
    PYCHOPPER_TOOL_VERSION=\$(pychopper --version 2>&1 | grep -oP '(\\d+\\.)+\\d+' || echo "no_disponible")
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pychopper: \${PYCHOPPER_TOOL_VERSION}
    END_VERSIONS

    echo "------------------------------------" | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Proceso completado para ${output_prefix_in}." | tee -a \$PYCHOPPER_PROCESS_LOG
    echo "Archivos finales generados:" | tee -a \$PYCHOPPER_PROCESS_LOG
    ls -lh *.gz 2>/dev/null | tee -a \$PYCHOPPER_PROCESS_LOG || echo "No se encontraron archivos .gz" | tee -a \$PYCHOPPER_PROCESS_LOG
    """

    stub:
    """
    echo "STUB: Creando archivos dummy para PYCHOPPER ${output_prefix_in}"
    touch "${output_prefix_in}.chopped.fastq.gz"
    touch "${output_prefix_in}.chopper_stats.tsv"
    touch "${output_prefix_in}.pychopper_script.log"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pychopper: stub_version
    END_VERSIONS
    """
}