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
  make_option(c("-c", "--coefficient"), type="numeric", action = "store", default=1.5,
              help="Define the coefficient for fence definition in 'Outliers' mode [default %default] [range = 1.5 - 3].", metavar="number"),
  make_option(c("-n", "--name"), type="character", action = "store", default=NULL,
              help="column name of the variable in the metadata file to be used.", metavar="string"),
  make_option(c("-v", "--value"), type="character", action = "store", default=NULL,
              help="value from the variable to be compare against the rest.", metavar="string"),
  make_option(c("-p", "--psignificant"), type="numeric", action = "store", default=0.05,
              help="value from the variable to be compare against the rest.", metavar="number")
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
    if (is.null(arguments$name) | is.null(arguments$value)) {
      stop("Error: a mandatory argument was not provided for 'Enrichment' mode.", call.=FALSE)
    }
  }
}

# Check software installation
suppressPackageStartupMessages(library(syntenet))
suppressPackageStartupMessages(library(mrfDepth))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(svglite))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(bc3net))

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
if (!requireNamespace("dplyr", quietly = TRUE)) {
   stop("Package \"dplyr\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("bc3net", quietly = TRUE)) {
   stop("Package \"bc3net\" not installed. Please install it to run this script.", call. = FALSE)
}

load_and_preprocess_data <- function(
    orthofinder_count = "Orthogroups.GeneCount-captainless.tsv",
    orthofinder_mcl = "Orthogroups-captainless.txt",
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
  Orthogroups <- read.table(orthofinder_mcl, sep = ":", col.names = c("Orthogroup","genes"))

  # Read metadata file if required
  if(arguments$mode == "Enrichment") {
    metadata <- read.csv(metadata_file, header = TRUE, sep = ";")
    col_classes <- c(rep("character", length(metadata)))
    metadata <- read.delim(
      metadata_file,
      header = TRUE,
      sep = ";",
      colClasses = col_classes
    )
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
  Outliers <- orthocounts[orthocounts > fence]

  if(length(Outliers) >= 1) {
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

process_enrichment <- function(
  annotation,
  orthogroups_mcl,
  metadata,
  variable_name,
  value
) {

  #Create gene set for Orthogroups enrichment 
  PreGene_list <- as.list(orthogroups_mcl$genes)
  names(PreGene_list) <- orthogroups_mcl$Orthogroup

  Gene_list <- lapply(PreGene_list, function(x) {
  general_vector <- unlist(strsplit(x, split = " "))
  general_vector <- general_vector[general_vector != ""]
  return(general_vector)
  })

  # Create reference gene vector
  Reference_gene <- unlist(Gene_list)

  # Create candidate gene vector
  ElementsIn <- metadata %>% filter(.data[[variable_name]] == value)
  Genes_elements <- annotation %>% filter( seq_id %in% ElementsIn$ElementID_updated) %>% select(ID) %>% unlist()
  Candidate_gene <- Reference_gene[Reference_gene %in% Genes_elements]

  # Run enrinchment analysis
  Enrichment_results <- enrichment(Candidate_gene, Reference_gene, Gene_list, adj = "fdr", verbose = FALSE)
  write.table(Enrichment_results, file = "Enrichment_results.txt",
    sep = '\t', row.names = F, col.names = T, quote = F)

  Significant_results <- Enrichment_results[Enrichment_results$padj <= arguments$psignificant,]

  if(length(Significant_results$padj) >= 1) {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","There has been enrich orthogroups identified","\n", sep=""))
    write.table(Significant_results, file = "Enrich_orthogroups.txt",
      sep = '\t', row.names = F, col.names = F, quote = F)
  }
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
  process_enrichment(
  annotation = data_list$annotation,
  orthogroups_mcl = data_list$orthogroups_mcl,
  metadata = data_list$metadata,
  variable_name = arguments$name,
  value = arguments$value)
}

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Overrepresentation analysis have finish.","\n", sep=""))
