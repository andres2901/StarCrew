# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory",metavar="PATH"),
  make_option(c("-m", "--mode"), type="character", action = "store", default=NULL,
              help="Mode of Overrepresentation analysis [accepted mode: Outliers, Enrichment]", metavar="string"),
  make_option(c("-a", "--approximation"), type="character", action = "store", default=NULL,
              help="Define the IQR approximation that is going to be used to defined outliers in 'Outliers' mode [accepted approximation: Standard, Skew].", metavar="string"),
  make_option(c("-n", "--coefficient"), type="float", action = "store", default=1.5,
              help="Define the coefficient for fence definition in 'Outliers' mode [default %default] [range = 1.5 - 3].", metavar="number"),
  make_option(c("-c", "--column"), type="character", action = "store", default=NULL,
              help="column name of the variable in the metadata file to be used.", metavar="string"),
  make_option(c("-v", "--value"), type="character", action = "store", default=NULL,
              help="value from the variable to be compare against the rest.", metavar="string")
)

# Parse the command-line arguments
arguments <- parse_args(OptionParser(option_list = option_list))

# Check for required arguments and flags
if (is.null(arguments$directory) | is.null(arguments$directory)) {
  stop("Error: a mandatory argument was not provided", call.=FALSE)
} else {
  setwd(arguments$directory)
}

# Check for argument mode and their sequential mandatory variables
if(arguments$mode != "Outliers" & arguments$mode != "Enrichment") {
  stop("Error: Provided mode is not accepted", call.=FALSE)
} else {
  if(arguments$mode == "Outliers") {
    if (is.null(arguments$approximation)) {
      stop("Error: the mandatory argument 'approximation' was not provided for 'Outliers' mode.", call.=FALSE)
    } else {
      if(arguments$approximation != "Standard" & arguments$approximation != "Skew") {
        stop("Error: Provided approximation is not accepted", call.=FALSE)
      }
    }
    if(arguments$coefficient < 1.5 | arguments$coefficient > 3){
      stop("Error: Provided coefficient for 'Otliers' mode is out of accepted range", call.=FALSE)
    }
  } else if(arguments$mode == "Enrichment"){
    if (is.null(arguments$column) | is.null(arguments$value)) {
      stop("Error: a mandatory argument was not provided for 'Enrichment' mode.", call.=FALSE)
    }
  }
}

# Check software installation
suppressPackageStartupMessages(library(syntenet))
suppressPackageStartupMessages(library(mrfDepth))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))


if (!requireNamespace("syntenet", quietly = TRUE)) {
   stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("mrfDepth", quietly = TRUE)) {
   stop("Package \"mrfDepth\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
   stop("Package \"ggplot2\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("svglite", quietly = TRUE)) {
   stop("Package \"svglite\" not installed. Please install it to run this script.", call. = FALSE)
}

load_and_preprocess_data <- function(
    orthofinder_count = "Orthogroups.GeneCount-captainless.tsv",
    orthofinder_mcl = "Orthogroups-captainlessID.txt",
    metadata_file = "metadata.csv"
) {

  # Read orthogrup counts and stay with the total number
  OrthoFinder <- read.table(
    orthofinder_count,
    header = TRUE,
    check.names = FALSE,
    row.names = 1
  )

  OrthoFinder2 <- as.vector(OrthoFinder$Total)
  names(OrthoFinder2) <- row.names(OrthoFinder)

  # Extract gene locations and format for gggenomes visualization
  annotation <- gff2GRangesList("Gff/")
  genes <- as.data.frame(unlist(annotation)) %>%
    filter(type == "gene") %>%
    mutate(
      seqnames = as.character(seqnames),
      type = "CDS",
      attribute = NA) %>%
    select(seq_id = seqnames, start, end, strand, type, attribute,ID) 

  # Read orthogroups MCL file
  Orthogroups <- read.table("Orthogroups.txt", sep = ":", col.names = c("Orthogroup","genes"))

  # Read metadata file if required
  if(arguments$mode == "Enrichment") {
    metadata <- read.csv(metadata_file, header = TRUE, sep = ";")
  } else {
    metadata <- vector() 
  }
  
  return(list(
    annotation = genes,
    metadata = metadata,
    orthogroups_mcl = Orthogroups,
    orthocounts = OrthoFinder2
  ))
}

process_outliers <- function(
  orthocounts,
  approximation,
  coefficient
) {

  # Calculate IQR
  Quartiles <- quantile(orthocounts, probs = c(0,0.25,0.5,0.75,1)) 
  IQR <- Quartiles[4] - Quartiles[2]

  # Calculate fence
  if(approximation == "Standard") {
    fence <- Quartiles[4] + (coefficient * IQR)
  } else if(approximation == "Skew") {
    MC <- medcouple(orthocounts, do.reflect = FALSE)
    fence <- Quartiles[4] + (coefficient * exp(3 * MC[1]) * IQR)
  }
  
  # Identify Outliers if any
  Outliers <- orthocounts[orthocounts >= fence]

  if(lenght(Outliers) >= ) {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","There has been outliers identified","\n", sep=""))
    Outliers_dataframe <- as.data.frame(Outliers)
    write.table(Outliers_dataframe, file = "Overrepresented_orthogroups.txt",
                        sep = '\t', row.names = T, col.names = F, quote = F)
  }

  

  # Create Figure
  orthocounts_dataframe <- as.data.frame(orthocounts)
  OrtSizeHist <- ggplot(orthocounts_dataframe, aes(x=orthocounts)) + geom_histogram(binwidth=1, fill="red") + geom_vline(aes(xintercept=fence), color="blue", linetype="dashed", linewidth=0.5) + xlab("Orthogroup total size")

  ggsave(OrtSizeHist, filename = "OrthogroupsSizeHistogram.svg", width = 14, height = 7)
}

# Start the process

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Reading data...","\n", sep=""))

data_list <- load_and_preprocess_data()

# Process data
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Processing data in '",arguments$mode, "' mode","\n", sep=""))

if(arguments$mode == "Outliers") {
  process_outliers(
    orthocounts = data_list$orthocounts,
    approximation= arguments$approximation,
    coefficient = arguments$coefficient)
} else if(arguments$mode == "Enrichment") {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Waiting to be develop","\n", sep=""))
}