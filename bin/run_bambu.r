#!/usr/bin/env Rscript
#######################################################################
# run_bambu.r
#######################################################################

suppressPackageStartupMessages({
    library(optparse)
    library(bambu)
})

# --------------------------------------------------------------------
# 1) Argumentos CLI
# --------------------------------------------------------------------
option_list <- list(
    make_option(c("-a", "--annotation"),  type="character", help="GTF de referencia (requerido)"),
    make_option(c("-f", "--fasta"),       type="character", help="Genoma FASTA (requerido)"),
    make_option(c("-o", "--tag"),         type="character",  default="bambu_output",
                help="Prefijo/carpeta de salida [default: %default]"),
    make_option(c("-n", "--ncore"),       type="integer",    default=1,
                help="Núcleos [default: %default]"),
    make_option(      "--ndr",            type="double",     default=NA,
                help="NDR manual 0-1 (AUTO si se omite)"),
    make_option(      "--discovery",   type="character", default=FALSE,
                help="Desactiva descubrimiento → cuantificación únicamente"),
    # ---------- UMBRALES ----------
    make_option(      "--min_readcount",  type="integer",    default=1,
                help="Lecturas mínimas por splice-site (descubrimiento)"),
    make_option(      "--min_txscore",    type="double",     default=0.2,
                help="Score mínimo para transcritos mono-exónicos"),
    make_option(      "--canonical_ss",   type="integer",    default=1,
                help="Soporte mínimo de splice-site canónico (si la versión lo soporta)")
)

opt       <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE)
args      <- opt$options
bam_files <- opt$args

if (args$discovery %in% c("TRUE", "True", "true", "T", "1")) {
    args$discovery <- TRUE
} else {
    args$discovery <- FALSE
}

# --------------------------------------------------------------------
# 2) Validaciones
# --------------------------------------------------------------------
if (!length(bam_files))             stop("Sin BAMs posicionales.")
if (!file.exists(args$annotation))  stop("GTF no encontrado.")
if (!file.exists(args$fasta))       stop("FASTA no encontrada.")

# --------------------------------------------------------------------
# 3) Log de inicio
# --------------------------------------------------------------------
cat("========== BAMBU ==========\n")
cat("Discovery :", ifelse(args$discovery, "TRUE", "FALSE"), "\n")
cat(sprintf("Umbrales   : min_readcount=%d | min_txscore=%.2f | canonical_ss=%d (si disponible)\n",
            args$min_readcount, args$min_txscore, args$canonical_ss))
cat("===========================\n")

# --------------------------------------------------------------------
# 4) Helper para armar opt.discovery y ejecutar bambu
# --------------------------------------------------------------------
run_bambu_once <- function(with_css = TRUE) {
    od <- list(
        min.readCount          = args$min_readcount,
        min.sampleNumber       = 1,
        min.txScore.singleExon = args$min_txscore,
        remove.subsetTx        = TRUE
    )
    if (with_css) {
        # Se intentará pasar canonicalSiteSupport solo si with_css==TRUE
        od$canonicalSiteSupport <- args$canonical_ss
    }
    bambu(
        reads        = bam_files,
        annotations  = args$annotation,
        genome       = args$fasta,
        discovery    = args$discovery,
        NDR          = if (!is.na(args$ndr)) args$ndr else NULL,
        ncore        = args$ncore,
        opt.discovery = od
    )
}

# --------------------------------------------------------------------
# 5) Ejecutar con fallback (CSS → sin CSS)
# --------------------------------------------------------------------
result <- tryCatch(
    run_bambu_once(with_css = TRUE),
    error = function(e) {
        msg <- conditionMessage(e)
        if (grepl("canonicalSiteSupport|Setting parameter that does not exist", msg, ignore.case = TRUE)) {
            message("Aviso: 'canonicalSiteSupport' no es reconocido por esta versión de bambu. Reintentando sin ese parámetro…")
            return(run_bambu_once(with_css = FALSE))
        }
        if (grepl("xgboost|Passed unrecognized parameters|Parameter 'data' has been renamed to 'x'|BiocParallel errors", msg, ignore.case = TRUE)) {
            stop(paste0(
                "Fallo compatible con xgboost≥2.0. ",
                "Asegura xgboost<2.0 (p.ej. 1.7.*) en el entorno Conda. ",
                "Error original: ", msg
            ))
        }
        stop(e)
    }
)

# --------------------------------------------------------------------
# 6) Filtro opcional (remover mono-exón novedosos)
# --------------------------------------------------------------------
row_data <- rowData(result)
monoexon_novel <- row_data$numExons == 1 & row_data$txClassDescription != "annotation"
result_filt <- result[!(row_data$TXNAME %in% row_data$TXNAME[monoexon_novel]), ]

# --------------------------------------------------------------------
# 7) Escritura de resultados
# --------------------------------------------------------------------
dir.create(args$tag, recursive = TRUE, showWarnings = FALSE)
writeBambuOutput(result,      path = args$tag)
writeBambuOutput(result_filt, path = paste0(args$tag, "_filtered"))
cat("Archivos exportados en: ", args$tag, "\n")
