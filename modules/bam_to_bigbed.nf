// =======================================
// MODULE: BAM_TO_BIGBED
// =======================================
nextflow.enable.dsl = 2

process BAM_TO_BIGBED {
    tag   "${meta.id}"
    label 'process_medium'

    publishDir "${params.outdir ?: 'results'}/bigbed", mode: 'copy',
      saveAs: { fn -> fn }

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                 'docker://longsplice:latest' : 'longsplice:latest' }"

    input:
    tuple val(meta), path(bam), path(chrom_sizes)

    output:
    path "${meta.id}.reads.bb", emit: bigbed
    path "versions.yml",        emit: versions

    script:
    def threads = task.cpus ?: 2
    """
    set -euo pipefail
    export LC_ALL=C

    # 1) Asegurar índice BAI
    if [ ! -e "${bam}.bai" ]; then
      samtools index -@ ${threads} ${bam}
    fi

    # 2) BAM -> BED12 por bloques de alineamiento
    bedtools bamtobed -i ${bam} -bed12 -split > ${meta.id}.reads.bed

    # 3) Ordenar por (chrom, start)
    sort -k1,1 -k2,2n ${meta.id}.reads.bed > ${meta.id}.reads.sorted.bed

    # 4) Rellenar a BED12 si faltan columnas
    awk 'BEGIN{OFS="\\t"}{
      if (NF<12){
        \$7=\$2; \$8=\$3;
        \$9="0,0,0";
        \$10=1;
        \$11=\$3-\$2;
        \$12=0;
      }
      print
    }' ${meta.id}.reads.sorted.bed > ${meta.id}.reads.bed12

    # 5) BED12 -> bigBed
    bedToBigBed -type=bed12 ${meta.id}.reads.bed12 ${chrom_sizes} ${meta.id}.reads.bb

    # 6) Versiones
    {
      echo "\"${task.process}\":"
      echo "  bedtools: \$(bedtools --version 2>/dev/null | head -n1 || echo unknown)"
      echo "  bedToBigBed: \$(bedToBigBed 2>&1 | head -n1 | rev | cut -d' ' -f1 | rev || echo unknown)"
      echo "  samtools: \$(samtools --version 2>/dev/null | head -n1 || echo unknown)"
    } > versions.yml
    """
}
