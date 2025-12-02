// =======================================
// SUBWORKFLOW: DO_GENOME
// =======================================

nextflow.enable.dsl = 2

include { GET_CHROM_SIZES } from '../modules/get_chrom_sizes.nf'
include { GTF2BED }        from '../modules/gtf2bed.nf'
include { SAMTOOLS_FAIDX as SAMTOOLS_FAIDX_GENOME } from '../modules/samtools_faidx.nf'

workflow DO_GENOME {
    take:
    ch_samples // tuple(meta, reads)

    main:
    //ch_versions = Channel.empty()

    // 1. Obtener referencias únicas (FASTA, GTF) del samplesheet
    // Usamos una combinación de fasta.toString() y gtf.toString() como clave para la unicidad
    ch_unique_ref_files = ch_samples
        .map { meta, reads -> tuple(meta.fasta.toString(), meta.gtf.toString(), meta.fasta, meta.gtf) }
        .unique { it[0] + it[1] } // Clave única basada en paths de fasta y gtf
        .map { fasta_str, gtf_str, fasta_path, gtf_path ->
             // Usar baseName del FASTA como ID consistente
             def ref_id = fasta_path.baseName
             [
                 id: ref_id,
                 fasta: fasta_path,
                 gtf: gtf_path,
                 is_transcripts: false, // Asumiendo genoma
                 annot: ref_id
             ]
         }

    // Separar canales basados en el ID de referencia
    ch_ref_meta = ch_unique_ref_files

    // 2. Procesos: Indexación FASTA y conversión GTF->BED (por referencia única)

    // Input para FAIDX: tuple(id, fasta_path)
    ch_faidx_input = ch_ref_meta.map { ref -> tuple(ref.id, ref.fasta) }
    SAMTOOLS_FAIDX_GENOME(ch_faidx_input)
    ch_fai = SAMTOOLS_FAIDX_GENOME.out.fai // Canal: tuple(id, fai_path)
    //ch_versions = ch_versions.mix(SAMTOOLS_FAIDX_GENOME.out.versions.ifEmpty(null))

    // Input para GTF2BED: tuple(id, gtf_path)
    ch_gtf2bed_input = ch_ref_meta.map { ref -> tuple(ref.id, ref.gtf) }
    GTF2BED(ch_gtf2bed_input)
    ch_bed = GTF2BED.out.gtf_bed // Canal: tuple(id, bed_path)
    //ch_versions = ch_versions.mix(GTF2BED.out.versions.ifEmpty(null))


    // Input para GET_CHROM_SIZES: tuple(id, fai_path)
    GET_CHROM_SIZES(ch_fai)
    ch_sizes = GET_CHROM_SIZES.out.sizes // Canal: tuple(id, sizes_path)
    //ch_versions = ch_versions.mix(GET_CHROM_SIZES.out.versions.ifEmpty(null))


    // 3. Combinar toda la información de referencia usando el ID como clave
    ch_ref_data_combined = ch_ref_meta
        .map { ref -> tuple(ref.id, ref) } // (id, ref_map)
        .join(ch_fai, by: 0)     // Une por id -> (id, ref_map, fai_path)
        .join(ch_bed, by: 0)     // Une por id -> (id, ref_map, fai_path, bed_path)
        .join(ch_sizes, by: 0)   // Une por id -> (id, ref_map, fai_path, bed_path, sizes_path)
        .map { id, ref_meta, fai_path, bed_path, sizes_path ->
            // Crear el objeto final de información de referencia
            // Se mantienen los nombres originales del mapa y se añaden los nuevos paths
            ref_meta + [ fai: fai_path, bed: bed_path, sizes: sizes_path ]
        }

    // 4. Emitir canales individuales (asumiendo una sola referencia para simplificar downstream)
    emit:
    ref_info     = ch_ref_data_combined              // Canal con un mapa [id, fasta, gtf, sizes, bed, fai, ...]
    fasta        = ch_ref_data_combined.map { it.fasta }.first() // Asume 1 ref, emite el path
    gtf          = ch_ref_data_combined.map { it.gtf }.first()   // Asume 1 ref, emite el path
    bed          = ch_ref_data_combined.map { it.bed }.first()   // Asume 1 ref, emite el path
    fai          = ch_ref_data_combined.map { it.fai }.first()   // Asume 1 ref, emite el path
    chrom_sizes  = ch_ref_data_combined.map { it.sizes }.first() // Asume 1 ref, emite el path
    //versions     = ch_versions.unique()
}
