#!/usr/bin/env Rscript
#######################################################################
# make_validated_samplesheet.R
# Toma como entrada el samplesheet.csv original (con columnas group, replicate, …)
# y construye samplesheet_validated.csv con las columnas:
#   sample, condition, counts_file
#
# Uso:
#   Rscript make_validated_samplesheet.R \
#       --input <ruta/al/samplesheet_original.csv> \
#       --output <ruta/al/samplesheet_validated.csv>
#######################################################################

suppressPackageStartupMessages({
    library(optparse)
})

# ----------------------------------------------------
# (1) Definir y parsear argumentos
# ----------------------------------------------------
option_list <- list(
    make_option(c("-i", "--input"), type = "character",
                help = "Ruta al samplesheet original (debe contener columnas 'group' y 'replicate')."),
    make_option(c("-o", "--output"), type = "character",
                help = "Ruta de salida para samplesheet_validated.csv.")
)
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if (is.null(opt$input) || !file.exists(opt$input)) {
    stop("Error: no se encontró el samplesheet original (--input).")
}
if (is.null(opt$output)) {
    stop("Error: debe especificar nombre de salida (--output).")
}

# ----------------------------------------------------
# (2) Leer el samplesheet original
# ----------------------------------------------------
df <- read.csv(opt$input, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

# Asegurarnos de que existan las columnas 'group' y 'replicate'
req_cols <- c("group", "replicate")
if (!all(req_cols %in% colnames(df))) {
    stop("El samplesheet original debe tener columnas EXACTAS: 'group' y 'replicate'.")
}

# ----------------------------------------------------
# (3) Construir columnas 'sample' y 'condition'
# ----------------------------------------------------
#    sample    = "<group>_<replicate>",  e.g. "WT_1", "KO_2", …
#    condition = group,                   e.g. "WT" o "KO"
df$sample    <- paste0(df$group, "_", df$replicate)
df$condition <- df$group

# ----------------------------------------------------
# (4) Montar la ruta a counts_transcript.txt de Bambu
# ----------------------------------------------------
#    Se asume que Bambu deposita los conteos en: bambu/<sample>/counts_transcript.txt
df$counts_file <- file.path("bambu", df$sample, "counts_transcript.txt")

# ----------------------------------------------------
# (5) Seleccionar sólo las columnas requeridas
# ----------------------------------------------------
#    El formato que exige IsoformSwitchAnalyzeR es:
#      sample, condition, counts_file
out_df <- df[, c("sample", "condition", "counts_file")]

# ----------------------------------------------------
# (6) Escribir el fichero validado
# ----------------------------------------------------
write.table(
    out_df,
    file      = opt$output,
    sep       = ",",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
)

cat("-> Se ha generado:", opt$output, "\n")
