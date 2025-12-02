// =======================================
// SUBWORKFLOW: DO_REPORTING
// =======================================

nextflow.enable.dsl = 2

include { MULTIQC } from '../modules/multiqc.nf'

def collectValidPaths(input_ch, accumulator_ch, is_tuple_with_path) {
    def target_ch = accumulator_ch
    if (input_ch) {
        def processed_ch = input_ch.map { item ->
            // Extraer el path: si es una tupla, tomar el segundo elemento, sino tomar el ítem directamente.
            def path_to_check = is_tuple_with_path ? item[1] : item
            def file_obj = file(path_to_check) // Convertir a objeto Path de Nextflow
            return file_obj.exists() ? file_obj : null // Devolver el objeto Path si existe, sino null
        }.filter { it != null } // Filtrar los nulos (paths a archivos/directorios no existentes)

        target_ch = target_ch.mix(processed_ch) // Mezclar los paths válidos con el acumulador
    }
    return target_ch
}

workflow DO_REPORTING {
    take:
        // Canales de entrada:
        ch_nanoplot_reports_pre     // tuple(meta, path_report_dir)
        ch_nanoplot_reports_post    // tuple(meta, path_report_dir)
        ch_chopper_logs             // path_log (path directo)
        ch_flagstat_reports         // tuple(meta, path_flagstat_file)
        ch_idxstats_reports         // tuple(meta, path_idxstats_file)
        ch_stats_reports            // tuple(meta, path_samtools_stats_file)
        // ch_bambu_counts_transcript
        // ch_bambu_gtf_extended
        // ch_dexseq_dtu_results
        // ch_bigwig_files
        ch_software_versions        // Canal que emite UNA LISTA de paths a versions.yml

    main:
        // Inicializar el canal que recolectará todos los paths válidos para MultiQC
        ch_valid_paths_for_multiqc = Channel.empty()

        // Recolectar solo los paths que existen de los canales deseados
        // El segundo argumento de collectValidPaths es el canal acumulador
        // El tercer argumento (is_tuple_with_path) es true si el canal emite tuplas (meta, path)
        ch_valid_paths_for_multiqc = collectValidPaths(ch_nanoplot_reports_pre, ch_valid_paths_for_multiqc, true)
        ch_valid_paths_for_multiqc = collectValidPaths(ch_nanoplot_reports_post, ch_valid_paths_for_multiqc, true)
        ch_valid_paths_for_multiqc = collectValidPaths(ch_chopper_logs, ch_valid_paths_for_multiqc, false)
        ch_valid_paths_for_multiqc = collectValidPaths(ch_flagstat_reports, ch_valid_paths_for_multiqc, true)
        ch_valid_paths_for_multiqc = collectValidPaths(ch_idxstats_reports, ch_valid_paths_for_multiqc, true)
        ch_valid_paths_for_multiqc = collectValidPaths(ch_stats_reports, ch_valid_paths_for_multiqc, true)

        // Manejo de ch_software_versions (que emite una lista de paths)
        if (ch_software_versions) {
            // 1. flatMap(): Convierte el canal que emite una lista en un canal que emite los elementos de la lista individualmente.
            // 2. map(): Para cada path individual, verifica si existe y lo convierte a objeto File.
            // 3. filter(): Elimina los nulos (archivos de versión no encontrados).
            // 4. mix(): Añade los paths de archivos de versión válidos al canal acumulador.
            def valid_version_files = ch_software_versions
                                        .flatMap() // Emite cada path de la lista como un ítem de canal
                                        .map { path_val -> def f = file(path_val); return f.exists() ? f : null }
                                        .filter { it != null }
            ch_valid_paths_for_multiqc = ch_valid_paths_for_multiqc.mix(valid_version_files)
        }

        // Recolecta todos los paths únicos y válidos en una LISTA.
        def multiqc_input_list = ch_valid_paths_for_multiqc
                                     .unique()    // Asegura que cada path solo se incluya una vez.
                                     .collect()   // Agrupa todos los paths emitidos por el canal en una sola lista.
                                     .ifEmpty([]) // Si no se encontraron paths válidos, pasa una lista vacía a MultiQC.

        // Ejecutar MultiQC, pasándole la lista de paths recolectados.
        // El módulo MULTIQC tomará esta lista y creará enlaces simbólicos para cada path.
        MULTIQC(multiqc_input_list)

    emit:
        multiqc_report = MULTIQC.out.report // Path al informe HTML de MultiQC
        multiqc_data   = MULTIQC.out.data   // Path al directorio de datos de MultiQC (opcional)
}