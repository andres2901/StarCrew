# ==============================================================================
# OVERREPRESENTATION ANALYSIS
# ==============================================================================
# Description : Performs overrepresentation analysis of orthogroups based on
#               OrthoFinder results. Supports two modes:
#               - Outliers  : identifies statistically large orthogroups via IQR
#               - Enrichment: tests orthogroup enrichment in an element subset
# Usage       : Rscript OverrepresentationAnalysis.R [options]
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

required_pkgs <- c("syntenet", "mrfDepth", "ggplot2", "svglite", "dplyr", "bc3net")

for (pkg in required_pkgs) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop(
      paste0("Package '", pkg, "' is not installed. Please install it before running this script."),
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(syntenet)
  library(mrfDepth)
  library(ggplot2)
  library(svglite)
  library(dplyr)
  library(bc3net)
})

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
    c("-m", "--mode"),
    type = "character", action = "store", default = NULL,
    help = "Analysis mode. Accepted values: Outliers, Enrichment.",
    metavar = "STRING"
  ),
  make_option(
    c("-c", "--coefficient"),
    type = "numeric", action = "store", default = 1.5,
    help = "Fence coefficient for outlier detection [default: %default] [range: 1.5 - 3].",
    metavar = "NUMBER"
  ),
  make_option(
    c("-n", "--name"),
    type = "character", action = "store", default = NULL,
    help = "Column name of the grouping variable in the metadata file (required for Enrichment mode).",
    metavar = "STRING"
  ),
  make_option(
    c("-v", "--value"),
    type = "character", action = "store", default = NULL,
    help = "Value of the grouping variable to compare against all others (required for Enrichment mode).",
    metavar = "STRING"
  ),
  make_option(
    c("-p", "--psignificant"),
    type = "numeric", action = "store", default = 0.05,
    help = "P-value significance threshold [default: %default].",
    metavar = "NUMBER"
  ),
  make_option(
    c("--padjust"),
    type = "character", action = "store", default = "bonferroni",
    help = "Multiple testing correction method [default: %default]. Accepted values: BH, BY, bonferroni, fdr, hochberg, holm, hommel.",
    metavar = "STRING"
  ),
  make_option(
    c("--countmode"),
    type = "character", action = "store", default = NULL,
    help = "Gene counting mode (required for Outliers mode). Accepted values: Gene, Ship.",
    metavar = "STRING"
  )
)

arguments <- parse_args(OptionParser(option_list = option_list))

# ==============================================================================
# ARGUMENT VALIDATION
# ==============================================================================

valid_modes        <- c("Outliers", "Enrichment")
valid_countmodes   <- c("Gene", "Ship")
valid_padjust      <- c("BH", "BY", "bonferroni", "fdr", "hochberg", "holm", "hommel")

# --- Directory ----------------------------------------------------------------
if (is.null(arguments$directory)) {
  stop("Argument '--directory' is required.", call. = FALSE)
}
setwd(arguments$directory)

# --- Mode ---------------------------------------------------------------------
if (is.null(arguments$mode) || !arguments$mode %in% valid_modes) {
  stop(
    paste0("Argument '--mode' must be one of: ", paste(valid_modes, collapse = ", "), "."),
    call. = FALSE
  )
}

# --- Mode-specific arguments --------------------------------------------------
if (arguments$mode == "Outliers") {

  if (arguments$coefficient < 1.5 || arguments$coefficient > 3) {
    stop("Argument '--coefficient' must be between 1.5 and 3.", call. = FALSE)
  }

  if (is.null(arguments$countmode) || !arguments$countmode %in% valid_countmodes) {
    stop(
      paste0("Argument '--countmode' is required for Outliers mode and must be one of: ",
             paste(valid_countmodes, collapse = ", "), "."),
      call. = FALSE
    )
  }

} else if (arguments$mode == "Enrichment") {

  if (is.null(arguments$name) || is.null(arguments$value)) {
    stop("Arguments '--name' and '--value' are required for Enrichment mode.", call. = FALSE)
  }

}

# --- P-adjust method ----------------------------------------------------------
if (!arguments$padjust %in% valid_padjust) {
  stop(
    paste0("Argument '--padjust' must be one of: ", paste(valid_padjust, collapse = ", "), "."),
    call. = FALSE
  )
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
# FUNCTION DEFINITIONS
# ==============================================================================

#' Load and preprocess OrthoFinder data
#'
#' Reads orthogroup gene counts, GFF annotations, orthogroup gene lists,
#' and (for Enrichment mode) a metadata file.
#'
#' @param orthofinder_count Character. Path to the OrthoFinder gene count table.
#' @param orthofinder_mcl   Character. Path to the OrthoFinder MCL groups file.
#' @param metadata_file     Character. Path to the metadata CSV file (semicolon-delimited).
#' @param count_mode        Character. Counting mode: "Gene" (total genes) or
#'                          "Ship" (binary presence per element).
#'
#' @return A named list with elements:
#'   \item{annotation}{data.frame of gene coordinates from GFF files}
#'   \item{metadata}{data.frame of sample metadata, or empty vector if not needed}
#'   \item{orthogroups_mcl}{data.frame with columns Orthogroup and genes}
#'   \item{orthocounts}{Named numeric vector of orthogroup sizes}
load_and_preprocess_data <- function(
  orthofinder_count = "Orthogroups.GeneCount-captainless.tsv",
  orthofinder_mcl   = "Orthogroups-captainless.txt",
  metadata_file     = "metadata.csv",
  count_mode
) {

  # --- Orthogroup counts ------------------------------------------------------
  ortho_raw <- read.table(orthofinder_count, header = TRUE, check.names = FALSE, row.names = 1)

  if (count_mode == "Gene") {
    orthocounts <- setNames(as.vector(ortho_raw$Total), rownames(ortho_raw))

  } else if (count_mode == "Ship") {
    ortho_binary <- ortho_raw[, !names(ortho_raw) %in% "Total"]
    ortho_binary[ortho_binary > 1] <- 1
    orthocounts <- rowSums(ortho_binary)
  }

  # --- GFF annotation ---------------------------------------------------------
  annotation_raw <- gff2GRangesList("Gff/")
  annotation <- as.data.frame(unlist(annotation_raw)) %>%
    filter(type == "gene") %>%
    mutate(
      seqnames  = as.character(seqnames),
      type      = "CDS",
      attribute = NA
    ) %>%
    select(seq_id = seqnames, start, end, strand, type, attribute, ID)

  # --- Orthogroup gene lists --------------------------------------------------
  orthogroups_mcl <- read.table(
    orthofinder_mcl,
    sep       = ":",
    col.names = c("Orthogroup", "genes")
  )

  # --- Metadata (Enrichment mode only) ----------------------------------------
  if (arguments$mode == "Enrichment") {
    col_count  <- ncol(read.csv(metadata_file, header = TRUE, sep = ";"))
    metadata   <- read.delim(
      metadata_file,
      header    = TRUE,
      sep       = ";",
      colClasses = rep("character", col_count)
    )
  } else {
    metadata <- vector()
  }

  return(list(
    annotation     = annotation,
    metadata       = metadata,
    orthogroups_mcl = orthogroups_mcl,
    orthocounts    = orthocounts
  ))
}


#' Identify overrepresented orthogroups using IQR-based outlier detection
#'
#' Computes an upper fence based on the interquartile range with skewness
#' correction. Orthogroups above the fence are written to disk and a
#' histogram is saved as an SVG file.
#'
#' @param orthocounts   Named numeric vector. Orthogroup sizes.
#' @param coefficient   Numeric. Multiplier for the IQR fence (1.5 - 3).
#'
#' @return NULL
process_outliers <- function(orthocounts, coefficient) {

  quartiles <- quantile(orthocounts, probs = c(0, 0.25, 0.5, 0.75, 1))
  iqr_val   <- quartiles["75%"] - quartiles["25%"]

  # Compute upper fence
  mc    <- medcouple(orthocounts, do.reflect = FALSE)
  fence <- quartiles["75%"] + ((coefficient * exp(3 * mc[1])) * iqr_val)

  outliers <- orthocounts[orthocounts > fence]

  if (length(outliers) >= 1) {
    log_message("Outliers identified.")
    outliers_df <- as.data.frame(outliers)
    write.table(
      outliers_df,
      file      = "Overrepresented_orthogroups.txt",
      sep       = "\t",
      row.names = TRUE,
      col.names = FALSE,
      quote     = FALSE
    )
  } else {
    log_message("No outliers detected above the computed fence.")
  }

  orthocounts_df <- as.data.frame(orthocounts)

  hist_plot <- ggplot(orthocounts_df, aes(x = orthocounts)) +
    geom_histogram(binwidth = 1, fill = "red") +
    geom_vline(
      aes(xintercept = fence),
      color    = "blue",
      linetype = "dashed",
      linewidth = 0.5
    ) +
    xlab("Orthogroup total size")

  ggsave(hist_plot, filename = "OrthogroupsSizeHistogram.svg", width = 14, height = 7)
}

#' Test orthogroup enrichment in a subset of elements
#'
#' Filters genes belonging to the target element subset and runs Fisher's
#' exact test (via bc3net::enrichment) against all orthogroups as background.
#' Criteria to identify enrich OGs: adjusted p-value below threshold AND 
#' at least 50% of genes in orthogroup.
#' Results are written to disk; significant hits are saved separately.
#'
#' @param annotation     data.frame. Gene coordinate table (output of load_and_preprocess_data).
#' @param orthogroups_mcl data.frame. Orthogroup gene-list table with columns Orthogroup, genes.
#' @param metadata       data.frame. Sample metadata table.
#' @param variable_name  Character. Column name in metadata used to define the target subset.
#' @param value          Character. Value in variable_name that defines the target subset.
#'
#' @return NULL
process_enrichment <- function(annotation, orthogroups_mcl, metadata, variable_name, value) {

  # --- Build named gene list per orthogroup -----------------------------------
  gene_list_raw <- setNames(
    as.list(orthogroups_mcl$genes),
    orthogroups_mcl$Orthogroup
  )

  gene_list <- lapply(gene_list_raw, function(x) {
    tokens <- unlist(strsplit(x, split = " "))
    tokens[tokens != ""]
  })

  # --- Define reference and candidate gene sets -------------------------------
  reference_genes  <- unlist(gene_list)

  target_elements  <- metadata %>% filter(.data[[variable_name]] == value)
  candidate_genes  <- annotation %>%
    filter(seq_id %in% target_elements$ElementID_updated) %>%
    pull(ID)

  candidate_genes  <- reference_genes[reference_genes %in% candidate_genes]

  # --- Run enrichment test ----------------------------------------------------
  enrichment_results <- enrichment(
    candidate_genes,
    reference_genes,
    gene_list,
    adj     = arguments$padjust,
    verbose = FALSE
  )

  write.table(
    enrichment_results,
    file      = "Enrichment_results.txt",
    sep       = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote     = FALSE
  )

  # --- Filter significant results ---------------------------------------------
  significant_results <- enrichment_results %>%
    filter(padj <= arguments$psignificant) %>%
    filter((genes / all) >= 0.5)

  if (nrow(significant_results) >= 1) {
    log_message("Enriched orthogroups identified.")
    write.table(
      significant_results,
      file      = "Enrich_orthogroups.txt",
      sep       = "\t",
      row.names = FALSE,
      col.names = FALSE,
      quote     = FALSE
    )
  } else {
    log_message("No significantly enriched orthogroups found.")
  }
}


# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

log_message("Reading and preprocessing data...")

data_list <- load_and_preprocess_data(count_mode = arguments$countmode)

log_message(paste0("Processing data in '", arguments$mode, "' mode..."))

if (arguments$mode == "Outliers") {
  process_outliers(
    orthocounts   = data_list$orthocounts,
    coefficient   = arguments$coefficient
  )

} else if (arguments$mode == "Enrichment") {
  process_enrichment(
    annotation      = data_list$annotation,
    orthogroups_mcl = data_list$orthogroups_mcl,
    metadata        = data_list$metadata,
    variable_name   = arguments$name,
    value           = arguments$value
  )
}

log_message("Done.")
