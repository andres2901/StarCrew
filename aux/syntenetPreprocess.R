# ==============================================================================
# DATA PREPROCESSING
# ==============================================================================
# Description : Preprocesses proteome and annotation data for synteny analysis.
#               Reads FASTA protein sequences and GFF annotation files, runs
#               syntenet's input processing pipeline, and exports the results.
# Usage       : Rscript SyntenetPreprocess.R <working_directory>
# Author      : Andres F. Lizcano Salas
# Date        : 11/Mar/2026
# ==============================================================================


# ==============================================================================
# DEPENDENCIES
# ==============================================================================

required_pkgs <- c("syntenet")

for (pkg in required_pkgs) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop(
      paste0("Package '", pkg, "' is not installed. Please install it before running this script."),
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages(library(syntenet))


# ==============================================================================
# ARGUMENT PARSING AND VALIDATION
# ==============================================================================

args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("A working directory path must be provided as a positional argument.", call. = FALSE)
}

dir_path <- as.character(args[1])

if (!dir.exists(dir_path)) {
  stop(paste0("The specified directory does not exist: ", dir_path), call. = FALSE)
}

setwd(dir_path)


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

gff_dir   <- grep("Gff",     internal_dirs, value = TRUE)
fasta_dir <- grep("Protein", internal_dirs, value = TRUE)

# Both directories must be present to proceed
if (length(gff_dir) == 0 || length(fasta_dir) == 0) {
  stop(
    "Required subdirectories 'Gff' and/or 'Protein' were not found in: ", dir_path,
    call. = FALSE
  )
}

proteomes  <- fasta2AAStringSetlist(fasta_dir)
annotation <- gff2GRangesList(gff_dir)

log_message("Processing data...")

pdata <- process_input(proteomes, annotation, gene_field = "ID")

log_message("Saving protein data...")

empty <- export_sequences(
  pdata$seq,
  outdir = file.path(getwd(), "PreprocessData")
)

log_message("Done.")
