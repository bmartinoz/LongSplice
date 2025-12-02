#!/usr/bin/env Rscript
# Author: Zifo Bioinformatics (nf-core/rnasplice)

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 4) stop("Usage: suppa_split_file.R <input_file> <samplesheet> <output_file_suffix> <calculate_ranges> <prefix>", call.=FALSE)

input_file         <- args[1]
samplesheet        <- args[2]
output_file_suffix <- args[3]
calculate_ranges   <- args[4]
prefix             <- ifelse(length(args) >= 5, args[5], "")

input_data <- read.csv(input_file, sep = "\t", header = TRUE, check.names = FALSE)
samplesheet <- read.csv(samplesheet, header = TRUE, check.names = FALSE)

if (!all(c("sample","condition") %in% colnames(samplesheet))) stop("Samplesheet must contain 'sample' and 'condition'.", call.=FALSE)
samplesheet <- samplesheet[,c("sample","condition")]
samplesheet <- samplesheet[!duplicated(samplesheet[,"sample"]),]
conditions  <- unique(samplesheet[,"condition"])

split_files <- function(condition, samplesheet, input_data, output_file_suffix, prefix, calculate_ranges){
    indices <- which(samplesheet$condition == condition)
    sample_names <- samplesheet[samplesheet$condition == condition,]$sample
    if (!all(samplesheet$sample %in% colnames(input_data))) stop("Input_file must contain samplesheet samples.", call.=FALSE)
    output_file <- ifelse(prefix=="", paste0(condition, output_file_suffix), paste0(prefix,"_",condition,output_file_suffix))
    write.table(input_data[,sample_names, drop=F], file = output_file, quote=FALSE, sep="\t")
    if (calculate_ranges == "TRUE") {
        range <- paste0(as.character(indices[1]), "-", as.character(indices[length(indices)]))
        if (indices[1] == 1) {
            cat(range, file = "ranges.txt", sep = "")
        } else {
            cat(c(",",range), file = "ranges.txt", sep = "", append=TRUE)
        }
    }
}

for (cond in conditions) split_files(cond, samplesheet, input_data, output_file_suffix, prefix, calculate_ranges)
sessionInfo()
