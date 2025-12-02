#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(openxlsx)
  library(ggplot2)
  library(scales)
})

option_list <- list(
  make_option(c("-i", "--input"), type="character", help="DEXSeq results file (.txt)"),
  make_option(c("-o", "--outprefix"), type="character", default="dtu_summary", help="Output file prefix [default %default]")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Verificar input
if (is.null(opt$input)) stop("Error: --input es obligatorio")

# Leer resultados
message("Leyendo archivo: ", opt$input)
dtu <- read_tsv(opt$input, show_col_types = FALSE)

# Chequear columnas requeridas
required_cols <- c("featureID", "groupID", "pvalue", "padj", "gene_adj_pvalue", "transcript_adj_pvalue")
missing <- setdiff(required_cols, colnames(dtu))
if (length(missing) > 0) stop("Faltan columnas necesarias: ", paste(missing, collapse=", "))

# Filtrar significativos
dtu_sig <- dtu %>% filter(transcript_adj_pvalue < 0.05)

# Top 20 por log2FC si está disponible
log2fc_col <- grep("^log2fold_", colnames(dtu), value = TRUE)
if (length(log2fc_col) > 0) {
  top20 <- dtu_sig %>% arrange(desc(abs(.data[[log2fc_col]]))) %>% head(20)
} else {
  top20 <- dtu_sig %>% head(20)
  log2fc_col <- NULL
}

# Guardar CSV y XLSX
write_csv(dtu, paste0(opt$outprefix, "_full.csv"))
write_csv(dtu_sig, paste0(opt$outprefix, "_significant.csv"))
write_csv(top20, paste0(opt$outprefix, "_top20.csv"))

wb <- createWorkbook()
addWorksheet(wb, "full"); writeData(wb, "full", dtu)
addWorksheet(wb, "significant"); writeData(wb, "significant", dtu_sig)
addWorksheet(wb, "top20"); writeData(wb, "top20", top20)
saveWorkbook(wb, paste0(opt$outprefix, ".xlsx"), overwrite = TRUE)

# Volcano plot si hay log2FC
if (!is.null(log2fc_col)) {
  dtu$log10padj <- -log10(dtu$transcript_adj_pvalue + 1e-10)
  dtu$is_significant <- dtu$transcript_adj_pvalue < 0.05

  p <- ggplot(dtu, aes_string(x=log2fc_col, y="log10padj", color="is_significant")) +
    geom_point(alpha=0.6) +
    scale_color_manual(values=c("grey", "red")) +
    labs(x="Log2 Fold Change", y="-log10 adjusted p-value", color="Significant") +
    theme_minimal()

  ggsave(paste0(opt$outprefix, "_volcano.png"), plot = p, width=7, height=5)
}


message("Análisis finalizado.")
