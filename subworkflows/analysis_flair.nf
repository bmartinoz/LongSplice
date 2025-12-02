#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ===== Módulos FLAIR core =====
include { FLAIR_ALIGN                   } from '../modules/analysis_flair/flair/align/main.nf'
include { FLAIR_CORRECT                 } from '../modules/analysis_flair/flair/correct/main.nf'
include { CAT_FASTQ                     } from '../modules/analysis_flair/cat_fastq/main.nf'
include { CONCAT_CORRECTED_BED          } from '../modules/analysis_flair/concat_corrected_bed/main.nf'
include { FLAIR_COLLAPSE                } from '../modules/analysis_flair/flair/collapse/main.nf'
include { FLAIR_QUANTIFY                } from '../modules/analysis_flair/flair/quantify/main.nf'
include { FLAIR_DIFFEXP                 } from '../modules/analysis_flair/flair/diffexp/main.nf'
include { FLAIR_DIFFSPLICE              } from '../modules/analysis_flair/flair/diffsplice/main.nf'


workflow ANALYSIS_FLAIR {

    take:
        reads             // tuple de meta, fastq
        genome                // genoma referencia
        genome_index           //indice del genoma
        gen_annot             // anotacion de referencia

    main:


    FLAIR_ALIGN(
        reads,
        genome,
        genome_index,
        gen_annot,
    )

    FLAIR_CORRECT(
        FLAIR_ALIGN.out.bed,
        genome,
        gen_annot,
    )

    // crear lo nuevos archivos concatenados
    //
    def new_meta = [id: "concat_fastq", single_end: true]

    def concat_fastq = reads
    .map { meta, fastq -> fastq }
    .collect()
    .map { fastq -> tuple(new_meta, fastq.flatten()) }

    def concat_corrected_bed = FLAIR_CORRECT.out.corrected_bed
    .map { meta, bed -> bed }
    .collect()
    .map { bed -> tuple(new_meta, bed.flatten()) }
    
    CAT_FASTQ(
        concat_fastq
    )

    CONCAT_CORRECTED_BED(
        concat_corrected_bed
    )

    FLAIR_COLLAPSE(
        genome,
        gen_annot,
        CAT_FASTQ.out.reads,
        CONCAT_CORRECTED_BED.out.combined_bed
    )

    def manifest_tsv = reads
        .collectFile(
            name: 'manifest.tsv',
            newLine: true
        ) { item ->
            def meta = item[0]
            def fastq = item[1]
            "${meta.id}\t${meta.condition}\t${meta.replicate}\t${fastq}"
        }

    FLAIR_QUANTIFY(
        FLAIR_COLLAPSE.out.isoforms_bed,
        FLAIR_COLLAPSE.out.isoforms_fa,
        manifest_tsv
    )

    FLAIR_DIFFEXP(
        FLAIR_QUANTIFY.out.quantify_counts
    )

    FLAIR_DIFFSPLICE(
        FLAIR_COLLAPSE.out.isoforms_bed,
        FLAIR_QUANTIFY.out.quantify_counts
    )







    emit:

    versions       = FLAIR_ALIGN.out.versions





}
