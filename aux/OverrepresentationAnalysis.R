# ==============================================================================
# OVERREPRESENTATION ANALYSIS SCRIPT
# ==============================================================================
# This script performs the overrepresentation analysis of Orthogroups based on 
# Orthofinder results.
# ==============================================================================

# Argument parsing and validation
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory", metavar="PATH"),
  make_option(c("-m", "--mode"), type="character", action = "store", default=NULL,
              help="Mode of Overrepresentation analysis [accepted mode: Outliers, Enrichment]", metavar="string"),
  make_option(c("-a", "--approximation"), type="character", action = "store", default=NULL,
              help="Define the IQR approximation used for outliers [accepted: Standard, Skew].", metavar="string"),
  make_option(c("-c", "--coefficient"), type="numeric", action = "store", default=1.5,
              help="Coefficient for fence definition [default %default] [range 1.5 - 3].", metavar="number"),
  make_option(c("-n", "--name"), type="character", action = "store", default=NULL,
              help="column name of the variable in the metadata file to be used.", metavar="string"),
  make_option(c("-v", "--value"), type="character", action = "store", default=NULL,
              help="value from the variable to be compared against the rest.", metavar="string"),
  make_option(c("-p", "--psignificant"), type="numeric", action = "store", default=0.05,
              help="P-value significance threshold.", metavar="number"),
  make_option(c("-cm", "--countmode"), type="character", action = "store", default="",
              help="Counting mode [accepted: Gene, Ships].", metavar="string")
)

arguments <- parse_args(OptionParser(option_list = option_list))

if (is.null(arguments$directory)) {
  stop("Error: directory must be provided", call. = FALSE)
} else {
  setwd(arguments$directory)
}

if (arguments$mode != "Outliers" & arguments$mode != "Enrichment") {
  stop("Error: Provided mode is not accepted", call. = FALSE)
} else {
  if (arguments$mode == "Outliers") {
    if (is.null(arguments$approximation)) {
      stop("Error: 'approximation' is mandatory for 'Outliers' mode.", call. = FALSE)
    } else {
      if (arguments$approximation != "Standard" & arguments$approximation != "Skew") {
        stop("Error: Provided approximation is not accepted", call. = FALSE)
      }
    }
    if (arguments$coefficient < 1.5 | arguments$coefficient > 3) {
      stop("Error: Coefficient out of accepted range (1.5 - 3)", call. = FALSE)
    }
    if (arguments$countmode != "Gene" & arguments$countmode != "Ships") {
      stop("Error: Provided countmode is not accepted for Outliers", call. = FALSE)
    }
  } else if (arguments$mode == "Enrichment") {
    if (is.null(arguments$name) | is.null(arguments$value)) {
      stop("Error: 'name' and 'value' are mandatory for 'Enrichment' mode.", call. = FALSE)
    }
  }
}

# Software Dependency Checks
suppressPackageStartupMessages({
  library(syntenet)
  library(mrfDepth)
  library(ggplot2)
  library(svglite)
  library(dplyr)
  library(bc3net)
})

required_pkgs <- c("syntenet", "mrfDepth", "ggplot2", "svglite", "dplyr", "bc3net")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package \"", pkg, "\" not installed. Please install it to run this script."), call. = FALSE)
  }
}

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

load_and_preprocess_data <- function(
  orthofinder_count = "Orthogroups.GeneCount-captainless.tsv",
  orthofinder_mcl = "Orthogroups-captainless.txt",
  metadata_file = "metadata.csv",
  count_mode
) {

  OrthoFinder <- read.table(orthofinder_count, header = TRUE, check.names = FALSE, row.names = 1)

  if (count_mode == "Gene") {
    OrthoFinder2 <- as.vector(OrthoFinder$Total)
    names(OrthoFinder2) <- row.names(OrthoFinder)
  } else if (count_mode == "Ships") {
    OrthoFinder <- OrthoFinder[, !names(OrthoFinder) %in% c("Total")]
    OrthoFinder[OrthoFinder > 1] <- 1
    OrthoFinder2 <- rowSums(OrthoFinder)
  }

  annotation <- gff2GRangesList("Gff/")
  genes <- as.data.frame(unlist(annotation)) %>%
    filter(type == "gene") %>%
    mutate(seqnames = as.character(seqnames), type = "CDS", attribute = NA) %>%
    select(seq_id = seqnames, start, end, strand, type, attribute, ID) 

  Orthogroups <- read.table(orthofinder_mcl, sep = ":", col.names = c("Orthogroup", "genes"))

  if (arguments$mode == "Enrichment") {
    metadata <- read.csv(metadata_file, header = TRUE, sep = ";")
    col_classes <- c(rep("character", length(metadata)))
    metadata <- read.delim(metadata_file, header = TRUE, sep = ";", colClasses = col_classes)
  } else {
    metadata <- vector() 
  }
  
  return(list(annotation = genes, metadata = metadata, orthogroups_mcl = Orthogroups, orthocounts = OrthoFinder2))
}

process_outliers <- function(orthocounts, approximation, coefficient) {
  Quartiles <- quantile(orthocounts, probs = c(0, 0.25, 0.5, 0.75, 1)) 
  IQR_val <- Quartiles[4] - Quartiles[2]

  if (approximation == "Standard") {
    fence <- Quartiles[4] + (coefficient * IQR_val)
  } else if (approximation == "Skew") {
    MC <- medcouple(orthocounts, do.reflect = FALSE)
    fence <- Quartiles[4] + (coefficient * exp(3 * MC[1]) * IQR_val)
  }
  
  Outliers <- orthocounts[orthocounts > fence]

  if (length(Outliers) >= 1) {
    cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Outliers identified", "\n", sep = ""))
    Outliers_dataframe <- as.data.frame(Outliers)
    write.table(Outliers_dataframe, file = "Overrepresented_orthogroups.txt", sep = '\t', row.names = T, col.names = F, quote = F)
  }

  orthocounts_dataframe <- as.data.frame(orthocounts)
  OrtSizeHist <- ggplot(orthocounts_dataframe, aes(x = orthocounts)) + 
    geom_histogram(binwidth = 1, fill = "red") + 
    geom_vline(aes(xintercept = fence), color = "blue", linetype = "dashed", linewidth = 0.5) + 
    xlab("Orthogroup total size")

  ggsave(OrtSizeHist, filename = "OrthogroupsSizeHistogram.svg", width = 14, height = 7)
}

process_enrichment <- function(annotation, orthogroups_mcl, metadata, variable_name, value) {
  PreGene_list <- as.list(orthogroups_mcl$genes)
  names(PreGene_list) <- orthogroups_mcl$Orthogroup

  Gene_list <- lapply(PreGene_list, function(x) {
    vec <- unlist(strsplit(x, split = " "))
    return(vec[vec != ""])
  })

  Reference_gene <- unlist(Gene_list)
  ElementsIn <- metadata %>% filter(.data[[variable_name]] == value)
  Genes_elements <- annotation %>% filter(seq_id %in% ElementsIn$ElementID_updated) %>% select(ID) %>% unlist()
  Candidate_gene <- Reference_gene[Reference_gene %in% Genes_elements]

  Enrichment_results <- enrichment(Candidate_gene, Reference_gene, Gene_list, adj = "fdr", verbose = FALSE)
  write.table(Enrichment_results, file = "Enrichment_results.txt", sep = '\t', row.names = F, col.names = T, quote = F)

  Significant_results <- Enrichment_results[Enrichment_results$padj <= arguments$psignificant, ]

  if (nrow(Significant_results) >= 1) {
    cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Enriched orthogroups identified", "\n", sep = ""))
    write.table(Significant_results, file = "Enrich_orthogroups.txt", sep = '\t', row.names = F, col.names = F, quote = F)
  }
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Reading data...", "\n", sep = ""))

data_list <- load_and_preprocess_data(count_mode = arguments$countmode)

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Processing data in '", arguments$mode, "' mode", "\n", sep = ""))

if (arguments$mode == "Outliers") {
  process_outliers(
    orthocounts = data_list$orthocounts,
    approximation = arguments$approximation,
    coefficient = arguments$coefficient
  )
} else if (arguments$mode == "Enrichment") {
  process_enrichment(
    annotation = data_list$annotation,
    orthogroups_mcl = data_list$orthogroups_mcl,
    metadata = data_list$metadata,
    variable_name = arguments$name,
    value = arguments$value
  )
}
