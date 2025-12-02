// =======================================
// MODULE: MULTIQC
// =======================================

nextflow.enable.dsl = 2

process MULTIQC {
    label 'process_low'
    tag "MultiQC Report"

    publishDir "${params.outdir}/multiqc", mode: 'copy', overwrite: true,
        saveAs: { filename -> filename == 'versions.yml' ? filename : null }

    input:
    val list_of_input_paths

    output:
    path "multiqc_report.html", emit: report
    path "multiqc_data"       , emit: data, optional: true
    path "versions.yml"       , emit: versions

    script:
    def args = task.ext.args ?: ''
    def config_file_path = params.multiqc_config ? file(params.multiqc_config) : null
    def multiqc_config_arg = (config_file_path && config_file_path.exists()) ? "--config ${config_file_path.toAbsolutePath()}" : ""
    def multiqc_title_arg = params.multiqc_title ? "--title \"${params.multiqc_title}\"" : ""

    def staging_cmds_list = []
    def link_name_counter = [:]

    list_of_input_paths.each { path_item ->
        def input_file_obj = file(path_item)

        if (input_file_obj.exists()) {
            def baseName = input_file_obj.getName()
            def original_path_for_prefix = input_file_obj.toAbsolutePath().toString()
            
            def finalLinkName = baseName

            if (baseName == "versions.yml" || baseName == "PRE_TRIM_Report" || baseName == "POST_TRIM_Report") {
                link_name_counter[baseName] = (link_name_counter[baseName] ?: 0) + 1
                if (link_name_counter[baseName] > 1) {
                    finalLinkName = "${baseName}_${link_name_counter[baseName]-1}"
                }
                if (input_file_obj.isDirectory()) {
            
                    def parentPathObj = input_file_obj.getParent()
                    def parentDirName = parentPathObj ? parentPathObj.getName() : "src_${System.currentTimeMillis()}"
                    def safeParentDirName = parentDirName.replaceAll("[^a-zA-Z0-9_.-]", "_")
                    finalLinkName = "${safeParentDirName}_${baseName}"
                } else if (baseName == "versions.yml") { 
                    def parentPathObj = input_file_obj.getParent()
                    def parentDirName = parentPathObj ? parentPathObj.getName() : "src_${System.currentTimeMillis()}"
                    def safeParentDirName = parentDirName.replaceAll("[^a-zA-Z0-9_.-]", "_")
                    finalLinkName = "${safeParentDirName}_${baseName}"
                }
            }

            finalLinkName = finalLinkName.replaceAll('[:/\\\\]','_').replaceAll(' ','_')
            def tempFinalLinkName = finalLinkName
            def count = 1
            
            staging_cmds_list << "ln -sfn '${input_file_obj.toAbsolutePath()}' ./${tempFinalLinkName}"
        } else {
            staging_cmds_list << "echo \"WARN: Input for MultiQC (path: ${path_item}) was passed but does not exist. Skipping link.\" >&2"
        }
    }

    // Unir los comandos
    def final_staging_cmds = staging_cmds_list.unique().join("\n    ")

    """
    echo "Creando enlaces simbólicos en el directorio de trabajo para MultiQC..."
    if [ -n "${final_staging_cmds}" ]; then
        ${final_staging_cmds}
    else
        echo "No hay archivos de entrada válidos para enlazar para MultiQC."
    fi

    echo "Contenido del directorio de trabajo de MultiQC después de enlazar (pwd: \$(pwd)):"
    ls -lah .

    multiqc \\
        ${multiqc_config_arg} \\
        ${multiqc_title_arg} \\
        ${args} \\
        .  // MultiQC escanea el directorio actual

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g" | cut -d ' ' -f 2)
    END_VERSIONS
    """
}