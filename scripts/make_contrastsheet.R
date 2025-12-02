#!/usr/bin/env Rscript
#######################################################################
# make_contrastsheet.R
# Genera "contrasts.csv" a partir de un samplesheet_validated.csv
# que ya contiene columnas "sample", "condition" y "counts_file".
#
# El resultado será un CSV con dos columnas:
#    treatment,control
# p. ej. “KO,WT” si los dos valores únicos en “condition” son WT y KO.
#
# Uso:
#   Rscript make_contrastsheet.R --input samplesheet_validated.csv --output contrasts.csv
#######################################################################

suppressPackageStartupMessages({
  library(optparse)
})

# ----------------------------------------------------
# 1) Definir y parsear argumentos
# ----------------------------------------------------
option_list <- list(
  make_option(c("-i","--input"),  type="character",
              help="Ruta a samplesheet_validated.csv (debe tener columnas: 'sample','condition','counts_file')"),
  make_option(c("-o","--output"), type="character",
              help="Ruta para el archivo de salida: contrasts.csv")
)
parser <- OptionParser(option_list = option_list)
opt <- parse_args(parser)

if ( is.null(opt$input) || !file.exists(opt$input) ) {
  stop("Error: no se encontró el archivo samplesheet_validated.csv (--input).")
}
if ( is.null(opt$output) ) {
  stop("Error: debe especificar nombre de salida (--output).")
}

# ----------------------------------------------------
# 2) Leer samplesheet_validated.csv
# ----------------------------------------------------
df <- read.csv(opt$input, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

req_cols <- c("sample", "condition", "counts_file")
if (!all(req_cols %in% colnames(df))) {
  stop("El samplesheet_validated.csv debe tener columnas EXACTAS: 'sample','condition','counts_file'.")
}

# ----------------------------------------------------
# 3) Obtener condiciones únicas
# ----------------------------------------------------
conds <- unique(df$condition)
if (length(conds) != 2) {
  stop("Se esperaban exactamente 2 condiciones en 'condition' (p.ej. WT y KO); se encontró: ",
       paste(conds, collapse = ", "))
}

# Vamos a asignar arbitrariamente:
#   control = el primer valor de `conds`,
#   treatment = el segundo
control   <- conds[1]
treatment <- conds[2]

# ----------------------------------------------------
# 4) Armar el data.frame final y escribirlo
# ----------------------------------------------------
contrasts_df <- data.frame(
  treatment = treatment,
  control   = control,
  stringsAsFactors = FALSE
)
write.table(
  contrasts_df,
  file      = opt$output,
  sep       = ",",
  row.names = FALSE,
  col.names = TRUE,
  quote     = FALSE
)

cat("-> Se ha generado:", opt$output, "\n")
