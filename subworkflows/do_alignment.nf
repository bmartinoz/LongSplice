// =============================
// SUBWORKFLOW: DO_ALIGNMENT
// =============================

nextflow.enable.dsl = 2

include { MINIMAP2_INDEX } from '../modules/minimap2_index.nf'
include { MINIMAP2_ALIGN } from '../modules/minimap2_align.nf'

workflow DO_ALIGNMENT {
    take:
        fasta_path  // Canal con el path único al FASTA
        reads       // Canal con tuple(meta_completo, fastq) para cada muestra

    main:
        //ch_versions = Channel.empty()

        // 1. Indexar FASTA (Se ejecuta una vez)
        MINIMAP2_INDEX(fasta_path)
        ch_mmi = MINIMAP2_INDEX.out.mmi.map { it[1] } // Extraer solo el path .mmi
        //ch_versions = ch_versions.mix(MINIMAP2_INDEX.out.versions.ifEmpty(null))

        // 2. Combinar lecturas con el índice
        ch_reads_to_align = reads.combine(ch_mmi) // tuple(meta, fastq, index)

        // 3. Alineación minimap2 directamente a BAM
        MINIMAP2_ALIGN(ch_reads_to_align)
        ch_bam = MINIMAP2_ALIGN.out.align_bam
        //ch_versions = ch_versions.mix(MINIMAP2_ALIGN.out.versions.ifEmpty(null))

    emit:
        align_bam = ch_bam
        //versions  = ch_versions.unique()
}
