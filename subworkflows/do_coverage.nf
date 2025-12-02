// =======================================
// SUBWORKFLOW: DO_COVERAGE
// =======================================

nextflow.enable.dsl = 2

include { BAM_COVERAGE } from '../modules/bam_coverage.nf'

workflow DO_COVERAGE {
    take:
        ch_sortbam_indexed    // viene de DO_SORT_INDEX, con tres componentes: (meta, bam, index_path)

    main:
        // Mapear para quedarnos solo con (meta, bam), descartando el índice
        ch_for_cov = ch_sortbam_indexed.map { meta, bam, idx ->
            tuple(meta, bam)
        }

        // Lanzar BAM_COVERAGE recibiendo (meta, bam)
        BAM_COVERAGE(ch_for_cov)

        ch_bigwig  = BAM_COVERAGE.out.bigwig
        //ch_versions = BAM_COVERAGE.out.versions

    emit:
        bigwig  = ch_bigwig
        //versions = ch_versions
}
