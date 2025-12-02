#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
})

# ---------- utils ----------
pkg_avail <- function(pkg) suppressWarnings(requireNamespace(pkg, quietly = TRUE))
pkg_ver   <- function(pkg) if (pkg_avail(pkg)) as.character(utils::packageVersion(pkg)) else NA_character_

safe_pdf <- function(fname, expr, width = 10, height = 8) {
  try({ grDevices::pdf(fname, width = width, height = height); force(expr); grDevices::dev.off() }, silent = TRUE)
}

clean_sample_names <- function(v) {
  v <- sub("\\.(fastq|fq|bam|sam|cram)(\\.gz)?$", "", v, ignore.case = TRUE)
  v <- sub("\\.(sorted|dedup|trimmed|filtered)$", "", v, ignore.case = TRUE)
  v <- sub("\\..*$", "", v)
  v
}

to_numeric_df <- function(df) {
  for (j in seq_len(ncol(df))) {
    x <- df[[j]]
    if (is.numeric(x) || is.integer(x)) next
    if (is.character(x)) {
      x2 <- suppressWarnings(as.numeric(x))
      if (!all(is.na(x2)) || all(is.na(x))) df[[j]] <- x2
    }
  }
  df
}

coerce_counts_matrix <- function(df) {
  num_cols <- vapply(df, is.numeric, logical(1))
  if (!any(num_cols)) stop("No hay columnas numéricas con conteos; revisa el archivo de conteos.")
  df <- df[, num_cols, drop = FALSE]
  df[is.na(df)] <- 0
  df[df < 0]   <- 0
  m <- as.matrix(round(df))
  storage.mode(m) <- "integer"
  m
}

read_counts_bambu_like <- function(path) {
  df <- tryCatch(
    read.delim(path, sep = "\t", header = TRUE, check.names = FALSE,
               quote = "", comment.char = "", na.strings = c("NA","NaN","","Inf","-Inf"),
               stringsAsFactors = FALSE),
    error = function(e) NULL
  )
  if (is.null(df) || ncol(df) == 1) {
    df <- read.csv(path, header = TRUE, check.names = FALSE,
                   quote = "", comment.char = "", na.strings = c("NA","NaN","","Inf","-Inf"),
                   stringsAsFactors = FALSE)
  }

  id_cols_known <- c("TXNAME","txname","feature_id","TX_ID","TXID",
                     "Geneid","GENEID","gene_id","GENE_ID","gene_name","GENENAME")
  rn_col <- intersect(colnames(df), id_cols_known)

  if (length(rn_col) >= 1) {
    pref   <- intersect(rn_col, c("TXNAME","txname","feature_id","TX_ID","TXID","Geneid","gene_id","GENE_ID"))
    rn_use <- if (length(pref) >= 1) pref[1] else rn_col[1]
    rownames(df) <- make.unique(as.character(df[[rn_use]]))
    df[, rn_col] <- list(NULL)
  } else {
    if (!is.numeric(df[[1]])) {
      rownames(df) <- make.unique(as.character(df[[1]]))
      df[[1]] <- NULL
    }
  }
  if (ncol(df) == 0) stop("Sin columnas de muestras tras quitar IDs; revisa el archivo.")

  df <- to_numeric_df(df)
  colnames(df) <- clean_sample_names(colnames(df))
  coerce_counts_matrix(df)
}

read_counts_stringtie2 <- function(path) {
  cm <- read.delim(path, sep = "\t", header = TRUE, check.names = FALSE, skip = 1,
                   quote = "", comment.char = "", stringsAsFactors = FALSE)
  dropc <- intersect(colnames(cm), c("Chr","Start","End","Length","Strand"))
  if (length(dropc)) cm[, dropc] <- list(NULL)

  gid_col <- colnames(cm)[1]
  cm[[gid_col]] <- as.character(cm[[gid_col]])
  cm <- aggregate(cm[, -1, drop=FALSE], by = list(Geneid = cm[[gid_col]]), FUN = sum, na.rm = TRUE)
  rownames(cm) <- make.unique(cm$Geneid)
  cm$Geneid <- NULL

  cm <- to_numeric_df(cm)
  colnames(cm) <- clean_sample_names(colnames(cm))
  coerce_counts_matrix(cm)
}

# siempre añade columna respetando 0 filas
safe_add_col <- function(df, name, value = NA_character_) {
  df[[name]] <- if (nrow(df) == 0) vector(mode = typeof(value), length = 0) else value
  df
}

# ---------- args ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Uso: run_deseq2.r <quant_method> <counts_file> [samplesheet_csv] [contrast_var] [ref] [tgt] [gtf] [alpha] [lfc_thr] [lfc_shrink_method]", call. = FALSE)
}

qmethod      <- args[1]
counts_file  <- args[2]
samplesheet  <- ifelse(length(args) >= 3 && nzchar(args[3]), args[3], NA_character_)
contrast_var <- ifelse(length(args) >= 4 && nzchar(args[4]), args[4], "group")
ref_level    <- ifelse(length(args) >= 5 && nzchar(args[5]), args[5], NA_character_)
tgt_level    <- ifelse(length(args) >= 6 && nzchar(args[6]), args[6], NA_character_)
gtf_path     <- ifelse(length(args) >= 7 && nzchar(args[7]), args[7], NA_character_)
alpha        <- ifelse(length(args) >= 8 && nzchar(args[8]), as.numeric(args[8]), 0.05)
lfc_thr      <- ifelse(length(args) >= 9 && nzchar(args[9]), as.numeric(args[9]), 1)
shrink_meth  <- ifelse(length(args) >= 10 && nzchar(args[10]), args[10], "")

outfile_txt  <- "deseq2_results.txt"
outfile_csv  <- "deseq2_results.csv"

# ---------- leer conteos ----------
countMat <- switch(qmethod,
  "stringtie2" = read_counts_stringtie2(counts_file),
  "bambu"      = read_counts_bambu_like(counts_file),
  stop("Unknown quant_method: ", qmethod, ". Use 'bambu' o 'stringtie2'.")
)

# ---------- colData ----------
build_coldata <- function(countMat, samplesheet, contrast_var) {
  if (!is.na(samplesheet) && file.exists(samplesheet)) {
    ss <- read.csv(samplesheet, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)

    sample_cols <- intersect(colnames(ss), c("sample","sample_id"))
    if (length(sample_cols) >= 1) {
      ss$sample_id <- clean_sample_names(as.character(ss[[sample_cols[1]]]))
    } else if (all(c("group","replicate") %in% colnames(ss))) {
      ss$sample_id <- clean_sample_names(paste(ss$group, ss$replicate, sep = "-"))
    } else {
      stop("Samplesheet debe contener 'sample'/'sample_id' o ('group','replicate').")
    }

    if (!all(colnames(countMat) %in% ss$sample_id)) {
      missing <- setdiff(colnames(countMat), ss$sample_id)
      stop("Muestras en conteos no encontradas en samplesheet: ", paste(missing, collapse = ", "))
    }
    ss <- ss[match(colnames(countMat), ss$sample_id), , drop = FALSE]

    if (contrast_var %in% colnames(ss)) {
      fac <- factor(ss[[contrast_var]])
      colname_use <- contrast_var
    } else {
      if (!"group" %in% colnames(ss))
        stop("No se encontró '", contrast_var, "' ni 'group' en samplesheet.")
      fac <- factor(ss$group)
      colname_use <- "group"
    }
    colData <- data.frame(row.names = ss$sample_id, check.names = FALSE)
    colData[[colname_use]] <- fac
    attr(colData, "design_factor") <- colname_use
    return(colData)
  } else {
    samples <- colnames(countMat)
    group   <- sub("^([^_]+)_.*$", "\\1", samples)
    colData <- data.frame(group = factor(group), row.names = samples, check.names = FALSE)
    attr(colData, "design_factor") <- "group"
    return(colData)
  }
}

colData <- build_coldata(countMat, samplesheet, contrast_var)
design_factor <- attr(colData, "design_factor")

# ---------- DESeq ----------
dds_full <- DESeqDataSetFromMatrix(countData = countMat, colData = colData, design = as.formula(paste0("~", design_factor)))
if (!is.na(ref_level)) dds_full[[design_factor]] <- stats::relevel(dds_full[[design_factor]], ref = ref_level)

min_count   <- 10L
min_samples <- max(2L, floor(ncol(dds_full) / 2))
keep <- rowSums(counts(dds_full) >= min_count) >= min_samples

if (sum(keep) == 0) {
  message("No hay entidades tras el filtrado (>=10 en suficientes muestras). Generando salidas vacías.")
  libsizes <- pmax(colSums(countMat), 1)
  norm_counts <- sweep(countMat, 2, libsizes, "/") * 1e6
  norm_df <- data.frame(feature_id = rownames(norm_counts), norm_counts, check.names = FALSE)
  rownames(norm_df) <- NULL

  res_df <- data.frame(
    baseMean=numeric(0), log2FoldChange=numeric(0), lfcSE=numeric(0), stat=numeric(0),
    pvalue=numeric(0), padj=numeric(0), feature_id=character(0), gene_name=character(0),
    significant=character(0), stringsAsFactors = FALSE
  )

  utils::write.csv(res_df, file = outfile_csv, row.names = FALSE)
  write.table(res_df, file = outfile_txt, sep = "\t", quote = FALSE, row.names = TRUE)
  utils::write.csv(norm_df, file = "normalized_counts.csv", row.names = FALSE)
  saveRDS(dds_full, "dds.rds")
  safe_pdf("ma_plot.pdf",      { plot.new(); title("MA Plot: sin entidades tras filtrado") })
  safe_pdf("volcano_plot.pdf", { plot.new(); title("Volcano: sin entidades tras filtrado") })
  safe_pdf("pca_plot.pdf",     { plot.new(); title("PCA: no disponible (VST no calculable)") })
  safe_pdf("heatmap.pdf",      { plot.new(); title("Heatmap: sin entidades") })

  writeLines(c(
    "Resumen del análisis DESeq2",
    "Total entidades analizadas (post-filtro): 0",
    paste0("alpha: ", alpha, " | LFC threshold: ", lfc_thr),
    "Significativos (padj < alpha): 0",
    "Up-regulados: 0",
    "Down-regulados: 0",
    paste0("Diseño: ~", design_factor),
    paste0("Contraste: ", ifelse(!is.na(tgt_level), tgt_level, "(auto)"), " vs ", ifelse(!is.na(ref_level), ref_level, "(auto)")),
    "LFC shrink: none (no aplicable)"
  ), "analysis_summary.txt")

  vers_lines <- c(
    "DESEQ2:",
    paste0("  r-base: ", R.version.string),
    paste0("  bioconductor-deseq2: ", pkg_ver("DESeq2")),
    paste0("  apeglm: ", pkg_ver("apeglm")),
    paste0("  ashr: ", pkg_ver("ashr")),
    paste0("  enhancedvolcano: ", pkg_ver("EnhancedVolcano")),
    paste0("  rtracklayer: ", pkg_ver("rtracklayer")),
    paste0("  ggplot2: ", pkg_ver("ggplot2")),
    paste0("  pheatmap: ", pkg_ver("pheatmap")),
    paste0("  rcolorbrewer: ", pkg_ver("RColorBrewer"))
  )
  writeLines(vers_lines, "versions.yml")
  quit(save = "no", status = 0)
}

dds <- dds_full[keep, ]
dds <- DESeq(dds)

get_results <- function(dds, design_factor, ref_level, tgt_level, alpha) {
  k <- nlevels(colData(dds)[[design_factor]])
  if (!is.na(ref_level) && !is.na(tgt_level)) {
    if (!all(c(ref_level, tgt_level) %in% levels(colData(dds)[[design_factor]])))
      stop("Niveles de contraste no presentes en '", design_factor, "'.")
    results(dds, contrast = c(design_factor, tgt_level, ref_level), alpha = alpha)
  } else if (k == 2) {
    lvls <- levels(colData(dds)[[design_factor]])
    results(dds, contrast = c(design_factor, lvls[2], lvls[1]), alpha = alpha)
  } else {
    rn <- resultsNames(dds)
    results(dds, name = tail(rn, 1), alpha = alpha)
  }
}

res <- get_results(dds, design_factor, ref_level, tgt_level, alpha)

# si por alguna razón res es de 0 filas, emitir salidas vacías limpias
res_df_tmp <- as.data.frame(res)
if (nrow(res_df_tmp) == 0) {
  message("DESeq2 devolvió 0 filas en 'res'; emitiendo salidas vacías seguras.")
  res_df <- safe_add_col(data.frame(
    baseMean=numeric(0), log2FoldChange=numeric(0), lfcSE=numeric(0), stat=numeric(0),
    pvalue=numeric(0), padj=numeric(0), feature_id=character(0),
    stringsAsFactors = FALSE
  ), "gene_name", NA_character_)
  utils::write.csv(res_df, file = outfile_csv, row.names = FALSE)
  write.table(res_df, file = outfile_txt, sep = "\t", quote = FALSE, row.names = TRUE)
  utils::write.csv(as.data.frame(counts(dds, normalized = TRUE)), "normalized_counts.csv")
  saveRDS(dds, "dds.rds")
  safe_pdf("ma_plot.pdf",      { plot.new(); title("MA Plot: sin entidades en 'res'") })
  safe_pdf("volcano_plot.pdf", { plot.new(); title("Volcano: sin entidades en 'res'") })
  safe_pdf("pca_plot.pdf",     { plot.new(); title("PCA: disponible si VST") })
  writeLines(c(
    "Resumen del análisis DESeq2",
    "Total entidades analizadas (res): 0",
    paste0("alpha: ", alpha, " | LFC threshold: ", lfc_thr),
    "Significativos (padj < alpha): 0",
    paste0("Diseño: ~", design_factor),
    paste0("Contraste: ",
           ifelse(!is.na(tgt_level), tgt_level, "(auto)"),
           " vs ",
           ifelse(!is.na(ref_level), ref_level, "(auto)")),
    "LFC shrink: none"
  ), "analysis_summary.txt")
  vers_lines <- c(
    "DESEQ2:",
    paste0("  r-base: ", R.version.string),
    paste0("  bioconductor-deseq2: ", pkg_ver("DESeq2")),
    paste0("  apeglm: ", pkg_ver("apeglm")),
    paste0("  ashr: ", pkg_ver("ashr")),
    paste0("  enhancedvolcano: ", pkg_ver("EnhancedVolcano")),
    paste0("  rtracklayer: ", pkg_ver("rtracklayer")),
    paste0("  ggplot2: ", pkg_ver("ggplot2")),
    paste0("  pheatmap: ", pkg_ver("pheatmap")),
    paste0("  rcolorbrewer: ", pkg_ver("RColorBrewer"))
  )
  writeLines(vers_lines, "versions.yml")
  quit(save="no", status=0)
}

shrink_used <- ""
try({
  if (nzchar(shrink_meth) && pkg_avail(shrink_meth)) {
    if (!is.na(ref_level) && !is.na(tgt_level))
      res <- lfcShrink(dds, contrast = c(design_factor, tgt_level, ref_level), type = shrink_meth)
    else
      res <- lfcShrink(dds, coef = tail(resultsNames(dds), 1), type = shrink_meth)
    shrink_used <- shrink_meth
  } else if (pkg_avail("apeglm")) {
    if (!is.na(ref_level) && !is.na(tgt_level))
      res <- lfcShrink(dds, contrast = c(design_factor, tgt_level, ref_level), type = "apeglm")
    else
      res <- lfcShrink(dds, coef = tail(resultsNames(dds), 1), type = "apeglm")
    shrink_used <- "apeglm"
  }
}, silent = TRUE)

# ---------- anotación con GTF (robusta 0 filas) ----------
annotate_results <- function(df, gtf_path) {
  df <- as.data.frame(df)
  df$feature_id <- rownames(df)

  if (nrow(df) == 0) return(safe_add_col(df, "gene_name", NA_character_))
  if (is.na(gtf_path) || !file.exists(gtf_path) || !pkg_avail("rtracklayer"))
    return(safe_add_col(df, "gene_name", NA_character_))

  g <- rtracklayer::import(gtf_path)
  gdf <- as.data.frame(g)
  has <- function(nm) nm %in% names(gdf)

  gene_tbl <- NULL
  if (has("gene_id")) {
    gene_df <- gdf[, intersect(c("type","gene_id","gene_name"), names(gdf)), drop = FALSE]
    if (has("type")) gene_df <- gene_df[is.na(gene_df$type) | gene_df$type == "gene", , drop = FALSE]
    if (!"gene_name" %in% names(gene_df)) gene_df <- safe_add_col(gene_df, "gene_name", NA_character_)
    gene_df <- gene_df[!is.na(gene_df$gene_id) & nzchar(gene_df$gene_id), , drop = FALSE]
    if (nrow(gene_df) > 0) gene_tbl <- unique(gene_df[, c("gene_id","gene_name")])
  }

  tx_tbl <- NULL
  if (has("transcript_id")) {
    tx_df <- gdf[, intersect(c("type","transcript_id","gene_id","gene_name"), names(gdf)), drop = FALSE]
    if (has("type")) tx_df <- tx_df[is.na(tx_df$type) | tx_df$type == "transcript", , drop = FALSE]
    if (!"gene_id"   %in% names(tx_df)) tx_df <- safe_add_col(tx_df, "gene_id", NA_character_)
    if (!"gene_name" %in% names(tx_df)) tx_df <- safe_add_col(tx_df, "gene_name", NA_character_)
    tx_df <- tx_df[!is.na(tx_df$transcript_id) & nzchar(tx_df$transcript_id), , drop = FALSE]
    if (nrow(tx_df) > 0) tx_tbl <- unique(tx_df[, c("transcript_id","gene_id","gene_name")])
  }

  out <- df
  if (!is.null(gene_tbl)) {
    out <- merge(out, gene_tbl, by.x = "feature_id", by.y = "gene_id", all.x = TRUE)
  }
  if ((!"gene_name" %in% names(out)) || all(is.na(out$gene_name))) {
    if (!is.null(tx_tbl)) {
      out <- merge(out, tx_tbl, by.x = "feature_id", by.y = "transcript_id", all.x = TRUE, suffixes = c("", ".tx"))
      if ("gene_name.tx" %in% names(out)) {
        if ("gene_name" %in% names(out)) {
          out$gene_name <- ifelse(is.na(out$gene_name), out$gene_name.tx, out$gene_name)
          out$gene_name.tx <- NULL
        } else {
          names(out)[names(out) == "gene_name.tx"] <- "gene_name"
        }
      }
    }
  }
  out <- safe_add_col(out, "gene_name", NA_character_)
  out
}

res_df <- as.data.frame(res)
res_df <- annotate_results(res_df, gtf_path)

# flags y orden
res_df$significant <- ifelse(!is.na(res_df$padj) &
                               res_df$padj < alpha &
                               abs(res_df$log2FoldChange) >= lfc_thr,
                             ifelse(res_df$log2FoldChange > 0, "Up", "Down"),
                             "NS")

ord <- order(res_df$pvalue, na.last = TRUE)
res_df <- res_df[ord, , drop = FALSE]

# ---------- export ----------
utils::write.csv(res_df, file = outfile_csv, row.names = FALSE)
write.table(res_df, file = outfile_txt, sep = "\t", quote = FALSE, row.names = TRUE)

norm_counts <- counts(dds, normalized = TRUE)
utils::write.csv(as.data.frame(norm_counts), "normalized_counts.csv")

saveRDS(dds, "dds.rds")
vsd <- tryCatch(vst(dds, blind = FALSE), error = function(e) NULL)
if (!is.null(vsd)) saveRDS(vsd, "vsd.rds")

# ---------- plots ----------
safe_pdf("ma_plot.pdf", { plotMA(res, main = "MA Plot") })

if (!is.null(vsd)) {
  safe_pdf("pca_plot.pdf", { plotPCA(vsd, intgroup = design_factor) })
} else {
  safe_pdf("pca_plot.pdf", { plot.new(); title("PCA no disponible (VST falló)") })
}

if (pkg_avail("EnhancedVolcano")) {
  suppressPackageStartupMessages(library(EnhancedVolcano))
  safe_pdf("volcano_plot.pdf", {
    print(EnhancedVolcano::EnhancedVolcano(
      res_df,
      lab       = if ("gene_name" %in% names(res_df)) res_df$gene_name else res_df$feature_id,
      x         = "log2FoldChange",
      y         = "padj",
      pCutoff   = alpha,
      FCcutoff  = lfc_thr,
      title     = paste0("DE: ", ifelse(is.na(tgt_level), "", tgt_level), " vs ", ifelse(is.na(ref_level), "", ref_level))
    ))
  }, width = 12, height = 10)
} else if (pkg_avail("ggplot2")) {
  suppressPackageStartupMessages(library(ggplot2))
  dfv <- res_df; dfv$neglog10padj <- -log10(dfv$padj)
  safe_pdf("volcano_plot.pdf", {
    print(ggplot2::ggplot(dfv, ggplot2::aes(x = log2FoldChange, y = neglog10padj, color = significant)) +
            ggplot2::geom_point(alpha = 0.6, size = 1.2) +
            ggplot2::geom_hline(yintercept = -log10(alpha), linetype = "dashed") +
            ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed") +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = "Volcano plot", x = "log2FC", y = "-log10(padj)"))
  }, width = 12, height = 10)
} else {
  dfv_ok <- res_df[is.finite(res_df$log2FoldChange) & is.finite(res_df$padj), , drop = FALSE]
  y <- -log10(dfv_ok$padj)
  safe_pdf("volcano_plot.pdf", {
    plot(dfv_ok$log2FoldChange, y, xlab = "log2FC", ylab = "-log10(padj)", main = "Volcano plot (base)", pch = 16, cex = 0.6)
    abline(h = -log10(alpha), v = c(-lfc_thr, lfc_thr), lty = 2)
  }, width = 12, height = 10)
}

# ---------- resumen ----------
sig_genes  <- sum(!is.na(res_df$padj) & res_df$padj < alpha)
up_genes   <- sum(res_df$significant == "Up",   na.rm = TRUE)
down_genes <- sum(res_df$significant == "Down", na.rm = TRUE)

writeLines(c(
  "Resumen del análisis DESeq2",
  paste0("Total entidades analizadas (post-filtro): ", nrow(res_df)),
  paste0("alpha: ", alpha, " | LFC threshold: ", lfc_thr),
  paste0("Significativos (padj < alpha): ", sig_genes),
  paste0("Up-regulados: ", up_genes),
  paste0("Down-regulados: ", down_genes),
  paste0("Diseño: ~", design_factor),
  paste0("Contraste: ", ifelse(!is.na(tgt_level), tgt_level, "(auto)"),
         " vs ", ifelse(!is.na(ref_level), ref_level, "(auto)")),
  paste0("LFC shrink: ", ifelse(nchar(shrink_used) > 0, shrink_used, "none"))
), "analysis_summary.txt")

vers_lines <- c(
  "DESEQ2:",
  paste0("  r-base: ", R.version.string),
  paste0("  bioconductor-deseq2: ", pkg_ver("DESeq2")),
  paste0("  apeglm: ", pkg_ver("apeglm")),
  paste0("  ashr: ", pkg_ver("ashr")),
  paste0("  enhancedvolcano: ", pkg_ver("EnhancedVolcano")),
  paste0("  rtracklayer: ", pkg_ver("rtracklayer")),
  paste0("  ggplot2: ", pkg_ver("ggplot2")),
  paste0("  pheatmap: ", pkg_ver("pheatmap")),
  paste0("  rcolorbrewer: ", pkg_ver("RColorBrewer"))
)
writeLines(vers_lines, "versions.yml")
