# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory [default %default]",metavar="PATH"),
  make_option(c("-a", "--anchors"), type="integer", action = "store", default=8,
              help=" Number of minimum anchor points for syntenet to call a collinear reagions [default %default]", metavar="number"),
  make_option(c("-g", "--gaps"), type="integer", action = "store", default=8,
              help="Number of maximum allowed gaps between anchor points for syntenet to call a collinear region [default %default]", metavar="number"),
  make_option(c("-t", "--threads"), type="integer", action = "store", default=1,
              help="Number of threads for the analysis [default %default]", metavar="number")
  make_option(c("-e", "--evalue"), type="numeric", action = "store", default=0.00001,
              help="e-value for syntenet to call a collinear region [default %default]", metavar="number")
)

# Parse the command-line arguments
arguments <- parse_args(OptionParser(option_list = option_list))

# Check for required arguments and flags
if (is.null(arguments$directory)) {
  stop("Error: directory must be provided.", call.=FALSE)
}

# Check for threads multicore
if (arguments$threads > 1) {
  suppressPackageStartupMessages(library(future))
  plan("multicore", workers = arguments$threads)
}

# Check software installation
suppressPackageStartupMessages(library(syntenet))

if (!requireNamespace("syntenet", quietly = TRUE)) {
   stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
}
  
# Call directories in the Working directory
internal_dir <- list.dirs(arguments$directory, recursive = F)

# Import data
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Reading data..."," \n", sep=""))
gff_dir <- as.character(grep("Gff", internal_dir, value = TRUE))
fasta_dir <- as.character(grep("Protein", internal_dir, value = TRUE))
  
proteomes <- fasta2AAStringSetlist(fasta_dir)
annotation <- gff2GRangesList(gff_dir)
  
#preprocess files
pdata <- process_input(proteomes, annotation, gene_field = "ID")
  
#import blast results
blast_dir <- as.character(grep("DiamondResults", internal_dir, value = TRUE))
blast_list <- read_diamond(blast_dir)
  
#perform synteny block analysis
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]"," Processing data with syntenet..."," \n", sep=""))
intersyn <- interspecies_synteny(blast_list, pdata$annotation, inter_dir = paste(arguments$directory,"/","Collinearity", sep = ""), anchors = arguments$anchors, max_gaps = arguments$gaps, e_value = arguments$evalue)
