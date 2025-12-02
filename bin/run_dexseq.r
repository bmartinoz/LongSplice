#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(DRIMSeq))
suppressPackageStartupMessages(library(DEXSeq))
suppressPackageStartupMessages(library(stageR))
suppressPackageStartupMessages(library(BiocParallel))

## -------------------- CLI --------------------
option_list = list(
  make_option(c("-m","--method"), type="character", default="bambu"),
  make_option(c("-c","--counts"), type="character"),
  make_option(c("-s","--samplesheet"), type="character"),
  make_option(c("-a","--annotation"), type="character"),
  make_option(c("-o","--outfile"), type="character", default="dexseq_results.txt"),
  make_option(c("-n","--ncore"), type="integer", default=1),
  make_option(c("--min_gene_expr"), type="integer",  default=10),
  make_option(c("--min_feature_expr"), type="integer", default=10),
  make_option(c("--min_feature_prop"), type="double", default=0.1),
  make_option(c("--min_samps_gene_expr"), type="integer", default=NULL),
  make_option(c("--min_samps_feature_expr"), type="integer", default=NULL),
  make_option(c("--min_samps_feature_prop"), type="integer", default=NULL)
)
opt <- parse_args(OptionParser(option_list=option_list, description="Run DEXSeq/stageR DTU analysis."))

stopifnot(!is.null(opt$counts), !is.null(opt$samplesheet), !is.null(opt$annotation))
if(opt$method != "bambu") stop("Solo se admite --method bambu")
if(!file.exists(opt$counts)) stop("No existe --counts: ", opt$counts)
if(!file.exists(opt$samplesheet)) stop("No existe --samplesheet: ", opt$samplesheet)
if(!file.exists(opt$annotation)) stop("No existe --annotation: ", opt$annotation)

message("--- Parámetros ---")
message("method=", opt$method)
message("counts=", opt$counts)
message("samplesheet=", opt$samplesheet)
message("annotation=", opt$annotation)
message("outfile=", opt$outfile, "  ncore=", opt$ncore)

## -------------------- Utils: aplanar data.frame --------------------
flatten_df <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  for (j in seq_len(ncol(df))) {
    x <- df[[j]]
    # Columnas tipo lista -> colapsar a texto
    if (is.list(x)) {
      df[[j]] <- vapply(
        x,
        function(v) {
          if (length(v) == 0) return(NA_character_)
          if (length(v) == 1) return(as.character(v))
          paste(as.character(v), collapse = ";")
        },
        character(1)
      )
      next
    }
    # Clases comunes Bioconductor (*List)
    if (inherits(x, c("IntegerList","NumericList","CharacterList","CompressedList"))) {
      df[[j]] <- vapply(
        as.list(x),
        function(v) paste(as.character(v), collapse = ";"),
        FUN.VALUE = character(1)
      )
      next
    }
    # Otras clases no atómicas -> character
    if (!is.atomic(x)) {
      df[[j]] <- as.character(x)
    }
  }
  df
}

## -------------------- Leer counts --------------------
message("Leyendo matriz de conteos (Bambu)…")
count.matrix.raw <- data.frame(read.table(opt$counts, sep="\t", header=TRUE, check.names=FALSE, stringsAsFactors=FALSE))
if (ncol(count.matrix.raw) < 3)
  stop("Matriz de conteos con <3 columnas")

ok1 <- colnames(count.matrix.raw)[1] %in% c("TXNAME","feature_id")
ok2 <- colnames(count.matrix.raw)[2] %in% c("GENEID","gene_id")
if(!ok1 || !ok2) stop("Formato de Bambu inesperado (esperaba TXNAME/GENEID en 1-2)")

count.matrix <- count.matrix.raw
colnames(count.matrix)[1:2] <- c("feature_id","gene_id")
count_cols_start_idx <- 3
matrix_sample_ids_from_file <- colnames(count.matrix)[count_cols_start_idx:ncol(count.matrix)]
message("Muestras en counts (originales): ", paste(matrix_sample_ids_from_file, collapse=", "))

## limpiar sufijos típicos
clean_names <- function(x){
  x <- sub("^.*[\\\\/]", "", x)                                # basename
  x <- sub("\\.sorted(\\.bam)?$", "", x, ignore.case=TRUE)     # .sorted[.bam]
  x <- sub("\\.(bam|sam|cram)$", "", x, ignore.case=TRUE)      # .bam/.sam/.cram
  x <- sub("\\.dedup$", "", x, ignore.case=TRUE)               # .dedup
  x
}
matrix_sample_ids_cleaned <- clean_names(matrix_sample_ids_from_file)

if(!identical(matrix_sample_ids_from_file, matrix_sample_ids_cleaned)){
  message("Muestras en counts (limpias): ", paste(matrix_sample_ids_cleaned, collapse=", "))
  colnames(count.matrix)[count_cols_start_idx:ncol(count.matrix)] <- matrix_sample_ids_cleaned
}
matrix_sample_ids_for_validation <- colnames(count.matrix)[count_cols_start_idx:ncol(count.matrix)]

## asegurar numéricas
for(j in count_cols_start_idx:ncol(count.matrix)){
  count.matrix[,j] <- suppressWarnings(as.numeric(as.character(count.matrix[,j])))
}
count.matrix[is.na(count.matrix)] <- 0

## -------------------- Leer samplesheet (flexible) --------------------
message("Leyendo/procesando samplesheet…")
ss <- read.csv(opt$samplesheet, header=TRUE, stringsAsFactors=FALSE, check.names=FALSE)
if(nrow(ss) < 2) stop("Samplesheet con <2 filas")

cn <- tolower(colnames(ss))
pick_col <- function(opts){
  idx <- which(cn %in% tolower(opts))
  if(length(idx)) colnames(ss)[idx[1]] else NULL
}
col_sample_id <- pick_col(c("sample_id","sample","id","sampleid","sample_name","name","run"))
col_group     <- pick_col(c("group","condition","phenotype"))
col_repl      <- pick_col(c("replicate","rep","replica","replicate_id"))
if(is.null(col_group)) stop("El samplesheet necesita columna 'group' o 'condition'")

## construir sample_id_final
if(!is.null(col_sample_id)){
  sample_id_final <- trimws(ss[[col_sample_id]])
} else if(!is.null(col_repl)) {
  sample_id_final <- paste(trimws(ss[[col_group]]), trimws(ss[[col_repl]]), sep=".")
} else {
  g <- as.factor(trimws(ss[[col_group]]))
  rep_infer <- ave(seq_len(nrow(ss)), g, FUN=seq_along)
  sample_id_final <- paste(as.character(g), rep_infer, sep=".")
  message("No había sample ni replicate; se infirieron réplicas 1..n por grupo.")
}
condition_vec <- trimws(ss[[col_group]])
sampleinfo <- data.frame(
  sample_id = sample_id_final,
  condition = condition_vec,
  stringsAsFactors = FALSE
)

## Intersección con columnas de counts
counts_samples <- matrix_sample_ids_for_validation
missing <- setdiff(unique(sampleinfo$sample_id), counts_samples)
if(length(missing)){
  warning("Samplesheet contiene muestras no presentes en counts (se omiten): ",
          paste(missing, collapse=", "))
  sampleinfo <- sampleinfo[sampleinfo$sample_id %in% counts_samples, , drop=FALSE]
}
if(nrow(sampleinfo) < 2) stop("Menos de 2 muestras tras sheet∩counts")

## Orden y factores
sampleinfo <- sampleinfo[match(counts_samples[counts_samples %in% sampleinfo$sample_id],
                               sampleinfo$sample_id), , drop=FALSE]
sampleinfo <- sampleinfo[!is.na(sampleinfo$sample_id), , drop=FALSE]
sampleinfo$condition <- factor(sampleinfo$condition)
sampleinfo$sample    <- factor(sampleinfo$sample_id)  # requerido por la fórmula
rownames(sampleinfo) <- sampleinfo$sample_id

if(nlevels(sampleinfo$condition) < 2) stop("Se requieren >=2 niveles de condición")

message("sampleinfo:\n"); print(sampleinfo)
cond_levels <- levels(sampleinfo$condition)
lgcolName <- if (length(cond_levels) == 2)
  paste0("log2fold_", cond_levels[2], "_vs_", cond_levels[1]) else "log2FoldChange"
message("Columna LFC objetivo: ", lgcolName)

## -------------------- DRIMSeq (filtrado) --------------------
message("Filtrando con DRIMSeq…")
counts_for_drimseq <- count.matrix[, c("feature_id","gene_id", sampleinfo$sample_id)]
d <- dmDSdata(counts = counts_for_drimseq, samples = sampleinfo)

n_reps_min <- min(table(sampleinfo$condition))
n_samp_gene    <- if (!is.null(opt$min_samps_gene_expr))    opt$min_samps_gene_expr    else n_reps_min
n_samp_feature <- if (!is.null(opt$min_samps_feature_expr)) opt$min_samps_feature_expr else n_reps_min
n_samp_prop    <- if (!is.null(opt$min_samps_feature_prop)) opt$min_samps_feature_prop else n_reps_min

dFilter <- dmFilter(
  d,
  min_samps_feature_expr = n_samp_feature,
  min_samps_feature_prop = n_samp_prop,
  min_samps_gene_expr    = n_samp_gene,
  min_feature_expr       = opt$min_feature_expr,
  min_gene_expr          = opt$min_gene_expr,
  min_feature_prop       = opt$min_feature_prop
)

genes_before <- length(unique(counts(d)$gene_id))
genes_after  <- length(unique(counts(dFilter)$gene_id))
features_before <- nrow(counts(d))
features_after  <- nrow(counts(dFilter))
message(sprintf("Filtrado: Genes %d->%d, Features %d->%d", genes_before, genes_after, features_before, features_after))
if (features_after == 0) stop("Ningún feature quedó después del filtrado.")

counts_dexseq_input <- round(as.matrix(counts(dFilter)[, sampleinfo$sample_id]))
rownames(counts_dexseq_input) <- counts(dFilter)$feature_id
dex_feature_ids <- counts(dFilter)$feature_id
dex_gene_ids    <- counts(dFilter)$gene_id
dex_sample_data <- sampleinfo[colnames(counts_dexseq_input), , drop=FALSE]

## -------------------- DEXSeq --------------------
message("Ejecutando DEXSeq…")
BPPARAM <- if (opt$ncore > 1) MulticoreParam(workers = opt$ncore) else SerialParam()
formulaFullModel    <- ~ sample + exon + condition:exon
formulaReducedModel <- ~ sample + exon

if(!identical(colnames(counts_dexseq_input), rownames(dex_sample_data)))
  stop("colnames(counts) != rownames(sampleData)")

dxd <- DEXSeqDataSet(
  countData = counts_dexseq_input,
  sampleData = dex_sample_data,
  design = formulaFullModel,
  featureID = dex_feature_ids,
  groupID = dex_gene_ids
)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, formula=formulaFullModel, BPPARAM=BPPARAM)
dxd <- testForDEU(dxd, reducedModel=formulaReducedModel, fullModel=formulaFullModel, BPPARAM=BPPARAM)
dxd <- estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)
dxr_object <- DEXSeqResults(dxd, independentFiltering=FALSE)
if (is.null(dxr_object) || nrow(dxr_object) == 0) stop("DEXSeqResults vacío.")
dxr_df <- as.data.frame(dxr_object)

actual_lfc_col <- colnames(dxr_df)[grep("log2fold", tolower(colnames(dxr_df)))]
if (length(actual_lfc_col) >= 1) {
  colnames(dxr_df)[colnames(dxr_df) == actual_lfc_col[1]] <- lgcolName
}

if (!"featureID" %in% colnames(dxr_df)) {
  if (!is.null(rownames(dxr_df))) dxr_df$featureID <- rownames(dxr_df) else stop("No featureID en resultados.")
}
if (!"groupID" %in% colnames(dxr_df)) {
  gm <- data.frame(featureID = featureNames(dxd), groupID = groupIDs(dxd))
  dxr_df <- merge(dxr_df, unique(gm), by="featureID", all.x=TRUE)
}

## -------------------- stageR --------------------
message("Aplicando stageR…")
main_results_table <- dxr_df

if (all(c("pvalue","groupID","featureID") %in% colnames(dxr_df))) {
  dxr_ok <- dxr_df[!is.na(dxr_df$pvalue) & !is.na(dxr_df$groupID) & !is.na(dxr_df$featureID), ]
  if (nrow(dxr_ok) > 0) {
    pScreen <- tapply(dxr_ok$pvalue, dxr_ok$groupID, FUN = min, na.rm = TRUE)
    pScreen <- p.adjust(pScreen, method = "BH")

    pConfirmation <- matrix(dxr_ok$pvalue, ncol = 1)
    rownames(pConfirmation) <- as.character(dxr_ok$featureID)

    tx2gene <- unique(data.frame(
      TXNAME   = as.character(dxr_ok$featureID),
      GENENAME = as.character(dxr_ok$groupID),
      stringsAsFactors = FALSE
    ))

    commonF <- intersect(rownames(pConfirmation), tx2gene$TXNAME)
    pConfirmation <- pConfirmation[commonF, , drop = FALSE]
    tx2gene <- tx2gene[tx2gene$TXNAME %in% commonF, , drop = FALSE]

    commonG <- intersect(names(pScreen), tx2gene$GENENAME)
    pScreen <- pScreen[commonG]
    tx2gene <- tx2gene[tx2gene$GENENAME %in% commonG, , drop = FALSE]

    if (length(pScreen) && nrow(pConfirmation) && nrow(tx2gene)) {
      st <- stageR::stageRTx(
        pScreen = pScreen,
        pConfirmation = pConfirmation,
        pScreenAdjusted = TRUE,
        tx2gene = tx2gene
      )
      st <- stageR::stageWiseAdjustment(st, method = "dtu", alpha = 0.05)
      padj <- suppressWarnings(stageR::getAdjustedPValues(
        st, order = FALSE, onlySignificantGenes = FALSE
      ))

      # Normalizar tipos y nombres para evitar columnas-lista
      padj <- as.data.frame(padj, stringsAsFactors = FALSE)
      colnames(padj)[1:3] <- c("featureID_stageR","gene_adj_pvalue","transcript_adj_pvalue")
      padj$featureID_stageR      <- as.character(padj$featureID_stageR)
      padj$gene_adj_pvalue       <- suppressWarnings(as.numeric(padj$gene_adj_pvalue))
      padj$transcript_adj_pvalue <- suppressWarnings(as.numeric(padj$transcript_adj_pvalue))

      main_results_table$featureID <- as.character(main_results_table$featureID)
      main_results_table <- merge(
        as.data.frame(main_results_table, stringsAsFactors = FALSE),
        padj,
        by.x = "featureID",
        by.y = "featureID_stageR",
        all.x = TRUE,
        sort = FALSE
      )
    } else {
      main_results_table$gene_adj_pvalue       <- NA_real_
      main_results_table$transcript_adj_pvalue <- NA_real_
    }
  } else {
    main_results_table$gene_adj_pvalue       <- NA_real_
    main_results_table$transcript_adj_pvalue <- NA_real_
  }
} else {
  warning("Faltan columnas para stageR; se omite.")
  main_results_table$gene_adj_pvalue       <- NA_real_
  main_results_table$transcript_adj_pvalue <- NA_real_
}

## -------------------- Conteos normalizados (por SF) --------------------
message("Agregando conteos normalizados…")
raw_counts_matrix_for_norm <- counts_dexseq_input
sf <- tryCatch(sizeFactors(dxd), error=function(e) stop("No se pudieron obtener sizeFactors: ", e$message))
if(length(sf) != ncol(raw_counts_matrix_for_norm)){
  if(length(sf) > ncol(raw_counts_matrix_for_norm)){
    warning("Más sizeFactors que muestras; se truncarán.")
    sf <- sf[1:ncol(raw_counts_matrix_for_norm)]
  } else stop("Menos sizeFactors que muestras.")
}
names(sf) <- colnames(raw_counts_matrix_for_norm)
norm_counts <- sweep(raw_counts_matrix_for_norm, 2, sf, "/")
norm_df <- as.data.frame(norm_counts, check.names=FALSE)
norm_df$featureID <- rownames(norm_counts)

## merge final + APLANAR
if ("featureID" %in% colnames(main_results_table)) {
  final <- merge(main_results_table, norm_df, by="featureID", all.x=TRUE, sort=FALSE)
} else {
  final <- main_results_table
}
if (all(c("groupID","featureID") %in% colnames(final))) {
  final <- final[order(final$groupID, final$featureID), ]
} else if ("featureID" %in% colnames(final)) {
  final <- final[order(final$featureID), ]
}

# Aplanar para evitar 'unimplemented type list'
main_results_table <- flatten_df(main_results_table)
final <- flatten_df(final)

# dump rápido para depuración
capture.output(str(final), file = "dexseq_results.str.txt")

message("Escribiendo: ", opt$outfile)
write.table(final, file=opt$outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE, na="NA")
message("Listo.")
