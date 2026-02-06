# ==============================================================================
# DATA PREPROCESSING SCRIPT
# ==============================================================================

# Argument parsing and check
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0) {
  stop("Error: The working directory path must be provided as a positional argument.", call. = FALSE)
}

dir_path <- as.character(args[1])

if (!dir.exists(dir_path)) {
  stop(paste0("Error: The specified directory does not exist: ", dir_path), call. = FALSE)
}

setwd(dir_path)

# Software Dependency Checks
suppressPackageStartupMessages(library(syntenet))

required_pkgs <- c("syntenet")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package \"", pkg, "\" not installed. Please install it to run this script."), call. = FALSE)
  }
}

# ==============================================================================
# MAIN EXECUTION PIPELINE
# ==============================================================================

internal_dir <- list.dirs(getwd(), recursive = FALSE)

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Reading data...", "\n", sep = ""))

gff_dir <- as.character(grep("Gff", internal_dir, value = TRUE))
fasta_dir <- as.character(grep("Protein", internal_dir, value = TRUE))

if (length(fasta_dir) == 0 | length(gff_dir) == 0) {
  stop("Error: Required 'Gff' or 'Protein' directories not found in the working directory.", call. = FALSE)
}

proteomes <- fasta2AAStringSetlist(fasta_dir)
annotation <- gff2GRangesList(gff_dir)

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Processing data...", "\n", sep = ""))
pdata <- process_input(proteomes, annotation, gene_field = "ID")

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Saving protein data...", "\n", sep = ""))
temp <- export_sequences(
  pdata$seq, 
  outdir = paste0(getwd(), "/PreprocessData")
)
