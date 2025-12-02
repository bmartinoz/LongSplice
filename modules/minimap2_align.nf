// =======================================
// MODULE: MINIMAP2_ALIGN
// =======================================

nextflow.enable.dsl = 2

process MINIMAP2_ALIGN {
    tag "$meta.id"
    //label 'process_medium' 

    input:
    tuple val(meta), path(reads), path(index)

    output:
    tuple val(meta), path("${meta.id}.bam"), emit: align_bam
    path "versions.yml", emit: versions

    script:
    def preset = meta.technology.toUpperCase() == 'ONT' ? 'splice' :
                 (meta.technology.toUpperCase() == 'PACBIO' ? 'splice' : 'splice')
    def minimap2_opts = params.minimap2_opts ?: ''

    """
    minimap2 -ax ${preset} ${minimap2_opts} -t ${task.cpus} ${index} ${reads} \
        | samtools view -Sb - > ${meta.id}.bam

    if [ ! -s "${meta.id}.bam" ]; then
        echo "Error: Output BAM file ${meta.id}.bam is empty. Minimap2 or samtools view failed." >&2
        exit 1
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        minimap2: \$(minimap2 --version)
        samtools: \$(samtools --version | head -n1 | cut -d' ' -f2)
    END_VERSIONS
    """
    
    stub:
    """
    touch ${meta.id}.bam
    echo "${task.process}:" > versions.yml
    echo "  minimap2: stub" >> versions.yml
    echo "  samtools: stub" >> versions.yml
    """
}
