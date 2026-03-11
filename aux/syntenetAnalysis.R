# ==============================================================================
# SYNTENY ANALYSIS
# ==============================================================================
# Description : Performs interspecies synteny analysis using syntenet.
#               Reads preprocessed proteome and annotation data alongside
#               DIAMOND similarity results, and identifies collinear regions
#               across species.
# Usage       : Rscript SyntenetAnalysis.R [options]
# Author      : Andres F. Lizcano Salas
# Date        : 11/Mar/2026
# ==============================================================================


# ==============================================================================
# DEPENDENCIES
# ==============================================================================

if (!suppressMessages(requireNamespace("optparse", quietly = TRUE))) {
  stop(
    "Package 'optparse' is not installed. Please install it before running this script.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(optparse))


# ==============================================================================
# ARGUMENT PARSING
# ==============================================================================

option_list <- list(
  make_option(
    c("-d", "--directory"),
    type = "character", action = "store", default = NULL,
    help = "Path to the working directory.",
    metavar = "PATH"
  ),
  make_option(
    c("-a", "--anchors"),
    type = "integer", action = "store", default = 8,
    help = "Minimum number of anchor points for syntenet to call a collinear region [default: %default].",
    metavar = "NUMBER"
  ),
  make_option(
    c("-g", "--gaps"),
    type = "integer", action = "store", default = 8,
    help = "Maximum number of allowed gaps between anchor points [default: %default].",
    metavar = "NUMBER"
  ),
  make_option(
    c("-e", "--evalue"),
    type = "numeric", action = "store", default = 1e-5,
    help = "E-value threshold for calling a collinear region [default: %default].",
    metavar = "NUMBER"
  ),
  make_option(
    c("-t", "--threads"),
    type = "integer", action = "store", default = 1,
    help = "Number of threads for parallel processing [default: %default].",
    metavar = "NUMBER"
  )
)

arguments <- parse_args(OptionParser(option_list = option_list))


# ==============================================================================
# ARGUMENT VALIDATION
# ==============================================================================

if (is.null(arguments$directory)) {
  stop("Argument '--directory' is required.", call. = FALSE)
}

if (!dir.exists(arguments$directory)) {
  stop(paste0("The specified directory does not exist: ", arguments$directory), call. = FALSE)
}

setwd(arguments$directory)

if (arguments$threads < 1) {
  stop("Argument '--threads' must be a positive integer.", call. = FALSE)
}


# ==============================================================================
# DEPENDENCIES (continued)
# ==============================================================================

# future is only required when parallel processing is requested
required_pkgs <- "syntenet"

if (arguments$threads > 1) {
  required_pkgs <- c(required_pkgs, "future")
}

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0("Package '", pkg, "' is not installed. Please install it before running this script."),
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(syntenet)
  if (arguments$threads > 1) library(future)
})

if (arguments$threads > 1) {
  plan("multicore", workers = arguments$threads)
}


# ==============================================================================
# HELPER FUNCTION: TIMESTAMP MESSAGE
# ==============================================================================

#' Print a timestamped message to the console
#'
#' @param msg Character string. Message to display.
#'
#' @return NULL
log_message <- function(msg) {
  cat(paste0("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", msg, "\n"))
}


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

log_message("Reading data...")

# Locate required subdirectories within the working directory
internal_dirs <- list.dirs(getwd(), recursive = FALSE)

gff_dir      <- grep("Gff",           internal_dirs, value = TRUE)
fasta_dir    <- grep("Protein",       internal_dirs, value = TRUE)
blast_dir    <- grep("DiamondResults", internal_dirs, value = TRUE)

if (length(gff_dir) == 0 || length(fasta_dir) == 0) {
  stop(
    "Required subdirectories 'Gff' and/or 'Protein' were not found in: ", arguments$directory,
    call. = FALSE
  )
}

if (length(blast_dir) == 0) {
  stop(
    "Required subdirectory 'DiamondResults' was not found in: ", arguments$directory,
    call. = FALSE
  )
}

proteomes  <- fasta2AAStringSetlist(fasta_dir)
annotation <- gff2GRangesList(gff_dir)
pdata      <- process_input(proteomes, annotation, gene_field = "ID")
blast_list <- read_diamond(blast_dir)

log_message("Running interspecies synteny analysis with syntenet...")

intersyn <- interspecies_synteny(
  blast_list,
  pdata$annotation,
  inter_dir = file.path(getwd(), "Collinearity"),
  anchors   = arguments$anchors,
  max_gaps  = arguments$gaps,
  e_value   = arguments$evalue
)

log_message("Done.")