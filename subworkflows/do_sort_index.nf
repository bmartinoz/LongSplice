// =============================
// SUBWORKFLOW: DO_SORT_INDEX
// =============================

nextflow.enable.dsl = 2

include { SAMTOOLS_SORT  } from '../modules/samtools_sort.nf'
include { SAMTOOLS_INDEX } from '../modules/samtools_index.nf'

workflow DO_SORT_INDEX {
    take:
        ch_bam    // tuples (meta, bam) desde el alineamiento

    main:
        // 1) Pasar un prefijo al sort (si tu módulo lo usa vía task.ext)
        ch_bam_ext = ch_bam.map { meta, bam_file ->
            def prefix = "${meta.id}.sorted"
            tuple(meta, bam_file, [ prefix: prefix ])
        }

        // 2) Ordenar
        SAMTOOLS_SORT(ch_bam_ext)
        ch_sorted_bam = SAMTOOLS_SORT.out.bam
        //ch_versions   = SAMTOOLS_SORT.out.versions.ifEmpty(Channel.empty())

        // 3) Indexar (emite triple OBLIGATORIA)
        SAMTOOLS_INDEX(ch_sorted_bam)
        ch_indexed = SAMTOOLS_INDEX.out.indexed   // (meta, bam_sorted, bai)
        //ch_versions = ch_versions.mix( SAMTOOLS_INDEX.out.versions.ifEmpty(Channel.empty()) )

    emit:
        sortbam  = ch_indexed
        //versions = ch_versions.unique()
}
