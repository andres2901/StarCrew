# Check and Call directory argument
args <- commandArgs(trailingOnly = TRUE)
dir <- as.character(args[1])

# Check software installation
suppressPackageStartupMessages(library(syntenet))

if (!requireNamespace("syntenet", quietly = TRUE)) {
  stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
}

# Call directories in the Working directory
internal_dir <- list.dirs(dir, recursive = F)
  
# Import data
gff_dir <- as.character(grep("Gff", internal_dir, value = TRUE))
fasta_dir <- as.character(grep("Protein", internal_dir, value = TRUE))

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Reading data..."," \n", sep=""))
proteomes <- fasta2AAStringSetlist(fasta_dir)
annotation <- gff2GRangesList(gff_dir)
  
# Preprocess files
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Processing data..."," \n", sep=""))
pdata <- process_input(proteomes, annotation, gene_field = "ID")
  
#export processed sequences
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Saving protein data..."," \n", sep=""))
temp <- export_sequences(pdata$seq,outdir = paste(dir,"/","PreprocessData", sep = ""))
