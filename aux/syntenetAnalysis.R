# ==============================================================================
# SYNTENY ANALYSIS SCRIPT
# ==============================================================================
# This script performs the synteny analysis of elements with syntenet 
# ==============================================================================

# Argument parsing and validation
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory [default %default]", metavar="PATH"),
  make_option(c("-a", "--anchors"), type="integer", action = "store", default=8,
              help="Number of minimum anchor points for syntenet to call a collinear regions [default %default]", metavar="number"),
  make_option(c("-g", "--gaps"), type="integer", action = "store", default=8,
              help="Number of maximum allowed gaps between anchor points for syntenet [default %default]", metavar="number"),
  make_option(c("-e", "--evalue"), type="numeric", action = "store", default=0.00001,
              help="e-value for syntenet to call a collinear region [default %default]", metavar="number"),
  make_option(c("-t", "--threads"), type="integer", action = "store", default=1,
              help="Number of threads for the analysis [default %default]", metavar="number")
)

arguments <- parse_args(OptionParser(option_list = option_list))

if (is.null(arguments$directory)) {
  stop("Error: directory must be provided.", call. = FALSE)
} else {
  setwd(arguments$directory)
}

# Software Dependency Checks
suppressPackageStartupMessages({
  library(syntenet)
  library(future)
})

required_pkgs <- c("syntenet")
if (arguments$threads > 1) {
  required_pkgs <- c(required_pkgs, "future")
}

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package \"", pkg, "\" not installed. Please install it to run this script."), call. = FALSE)
  }
}

if (arguments$threads > 1) {
  plan("multicore", workers = arguments$threads)
}

# ==============================================================================
# MAIN EXECUTION
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

pdata <- process_input(proteomes, annotation, gene_field = "ID")

blast_dir <- as.character(grep("DiamondResults", internal_dir, value = TRUE))
if (length(blast_dir) == 0) {
  stop("Error: 'DiamondResults' directory not found.", call. = FALSE)
}
blast_list <- read_diamond(blast_dir)

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] Processing data with syntenet...", "\n", sep = ""))

intersyn <- interspecies_synteny(
  blast_list, 
  pdata$annotation, 
  inter_dir = paste0(getwd(), "/Collinearity"), 
  anchors = arguments$anchors, 
  max_gaps = arguments$gaps, 
  e_value = arguments$evalue
)
