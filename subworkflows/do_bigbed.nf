// =======================================
// SUBWORKFLOW: DO_BIGBED
// =======================================
nextflow.enable.dsl = 2

include { BAM_TO_BIGBED } from '../modules/bam_to_bigbed.nf'

workflow DO_BIGBED {
    take:
        ch_sortbam_indexed  // tuples: (meta, bam, idx) desde DO_SORT_INDEX
        ch_chrom_sizes      // path: *.sizes desde DO_GENOME

    main:
        // Quitar el índice del tuple
        ch_for_bb = ch_sortbam_indexed.map { meta, bam, idx -> tuple(meta, bam) }

        // Combinar cada (meta,bam) con el único sizes
        ch_in = ch_for_bb.combine(ch_chrom_sizes).map { meta, bam, sizes ->
            tuple(meta, bam, sizes)
        }

        BAM_TO_BIGBED(ch_in)

        ch_bigbed   = BAM_TO_BIGBED.out.bigbed
        //ch_versions = BAM_TO_BIGBED.out.versions

    emit:
        bigbed   = ch_bigbed
        //versions = ch_versions
}
