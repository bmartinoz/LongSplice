/*
 * Combine individual flair transcriptome outputs into a single set
 * Compatible with FLAIR quantify input requirements
 */
process COMBINE_FLAIR_TRANSCRIPTOME {
    publishDir "${params.outdir}/transcriptome", mode: 'copy'
    label 'process_low'
    
    input:
    path(bed_files)
    path(gtf_files)
    path(fa_files)
    
    output:
    path("combined_samples.flair.isoforms.bed"), emit: bed
    path("combined_samples.flair.isoforms.gtf"), emit: gtf
    path("combined_samples.flair.isoforms.fa"), emit: fa
    
    script:
    """
    # Combine BED files (remove duplicates by isoform ID - column 4)
    cat ${bed_files} | sort -k4,4 | awk '!seen[\$4]++' > combined_samples.flair.isoforms.bed
    
    # Combine GTF files (remove duplicates by attributes - column 9)
    cat ${gtf_files} | grep -v '^#' | sort -k1,1 -k4,4n | awk '!seen[\$9]++' > combined_samples.flair.isoforms.gtf
    
    # Combine FASTA files (remove duplicates manually)
    awk '/^>/ {if (seqid && !seen[seqid]++) {print header; print seq} header=\$0; seq=""; seqid=\$1; next} {seq=seq\$0} END {if (seqid && !seen[seqid]++) {print header; print seq}}' ${fa_files} > combined_samples.flair.isoforms.fa
    
    # Verificar archivos no vacíos
    for file in combined_samples.flair.isoforms.{bed,gtf,fa}; do
        if [ ! -s "\$file" ]; then
            echo "ERROR: \$file is empty!" >&2
            exit 1
        fi
    done
    
    echo "Combined isoforms successfully"
    """
}