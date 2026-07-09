# ==============================================================================
# CLUSTER ANALYSIS
# ==============================================================================
# Description : Performs comparative analysis of starships elements,
#               including hierarchical clustering of orthogroups, synteny
#               visualization, nesting detection, pangenome analysis,
#               and Captain vs Cargo discordance analysis.
# Usage       : Rscript ClusterAnalysis.R [options]
# Author      : Andres F. Lizcano Salas
# Date        : 07/Jul/2026
# ==============================================================================


# ==============================================================================
# DEPENDENCIES
# ==============================================================================

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop(
    "Package 'optparse' is not installed. Please install it before running this script.",
    call. = FALSE
  )
}
suppressPackageStartupMessages(library(optparse))

required_pkgs <- c(
  "ggplot2", "ape", "ggtree", "reshape2", "viridis", "dplyr",
  "gggenomes", "scales", "dendextend", "NbClust", "syntenet", "ggnewscale"
)

for (pkg in required_pkgs) {
  if (!suppressMessages(requireNamespace(pkg, quietly = TRUE))) {
    stop(
      paste0("Package '", pkg, "' is not installed. Please install it before running this script."),
      call. = FALSE
    )
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ape)
  library(reshape2)
  library(viridis)
  library(dplyr)
  library(gggenomes)
  library(scales)
  library(dendextend)
  library(NbClust)
  library(syntenet)
  library(ggnewscale)
  library(ggtree)
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
    c("-s", "--subclusters"),
    type = "integer", action = "store", default = 1,
    help = "Number of subclusters from the synteny analysis. Use 1 when no subclusters were identified [default: %default].",
    metavar = "NUMBER"
  ),
  make_option(
    c("-c", "--captainRemoval"),
    type = "integer", action = "store", default = 0,
    help = "Number of elements removed from the original cluster [default: %default].",
    metavar = "NUMBER"
  ),
  make_option(
    c("-p", "--pident"),
    type = "numeric", action = "store", default = 70,
    help = "Minimum percentage identity of BLAST results to include as links for nucleotide synteny visualization [default: %default] [range: 50 - 95].",
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

if (arguments$subclusters < 1) {
  stop("Argument '--subclusters' must be >= 1.", call. = FALSE)
}

if (arguments$pident < 50 || arguments$pident > 95) {
  stop("Argument '--pident' must be between 50 and 95", call. = FALSE)
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

#' Load and preprocess BLAST, annotation, and orthogroup data
#'
#' Reads BLAST results (filtered by identity), GFF annotations, orthogroup
#' tables, and optionally a subcluster definition file.
#'
#' @param blast_file Character. Path to the BLAST results file (tab-delimited, no header).
#' @param gff_file   Character. Path to the GFF annotation file (unused directly; GFFs
#'                   are loaded via gff2GRangesList from the "Gff/" subdirectory).
#' @param min_pident Numeric. Minimum percentage identity threshold for BLAST links.
#'
#' @return A named list with elements:
#'   \item{blast_results}{data.frame of filtered BLAST hits}
#'   \item{links}{data.frame of BLAST hits reformatted for gggenomes}
#'   \item{seqs}{data.frame of sequence lengths (one row per unique sequence)}
#'   \item{genes}{data.frame of gene coordinates from GFF files}
#'   \item{orthogroups_table}{data.frame with columns Orthogroup and genes}
#'   \item{subclusters}{data.frame of subcluster definitions, or empty vector}
load_and_preprocess_data <- function(
  blast_file = "Blast_CleanResults.txt",
  gff_file   = "Final_model.gff",
  min_pident
) {

  # --- BLAST results ----------------------------------------------------------
  blast_col_names   <- c("qseqid", "sseqid", "qstart", "qend", "sstart", "send",
                         "pident", "length", "qlen", "slen")
  blast_col_classes <- c("character", "character", rep("numeric", 8))

  blast_results <- read.delim(
    blast_file,
    header    = FALSE,
    sep       = "\t",
    col.names = blast_col_names,
    colClasses = blast_col_classes
  ) %>%
    filter(qseqid != sseqid, pident > min_pident)

  # Reformat BLAST hits for gggenomes link layer
  links <- blast_results %>%
    mutate(alignment_id = row_number()) %>%
    select(
      seq_id = qseqid, start = qstart, end = qend,
      seq_id2 = sseqid, start2 = sstart, end2 = send,
      pident, length
    )

  # Build sequence length table (one row per unique sequence, using max observed length)
  seqs <- bind_rows(
    blast_results %>% select(seq_id = qseqid, length = qlen),
    blast_results %>% select(seq_id = sseqid, length = slen)
  ) %>%
    distinct() %>%
    group_by(seq_id) %>%
    summarize(length = max(length)) %>%
    ungroup()

  # --- GFF annotation ---------------------------------------------------------
  annotation_raw <- gff2GRangesList("Gff/")
  genes <- as.data.frame(unlist(annotation_raw)) %>%
    filter(type == "gene") %>%
    mutate(
      seqnames  = as.character(seqnames),
      type      = "CDS",
      attribute = NA
    ) %>%
    select(seq_id = seqnames, start, end, strand, type, attribute, ID)

  # --- Orthogroup table -------------------------------------------------------
  orthogroups_table <- read.table(
    "Orthogroups.txt",
    sep       = ":",
    col.names = c("Orthogroup", "genes")
  )

  # --- Subcluster definitions (only when multiple subclusters requested) ------
  subcluster_definition <- vector()
  if (arguments$subclusters > 1) {
    subcluster_definition <- read.table(
      "Subclusters_MCL.txt",
      sep       = "\t",
      col.names = c("Subcluster", "Elements")
    )
  }

  return(list(
    blast_results   = blast_results,
    links           = links,
    seqs            = seqs,
    genes           = genes,
    orthogroups_table = orthogroups_table,
    subclusters     = subcluster_definition
  ))
}


#' Perform initial hierarchical clustering of orthogroup profiles
#'
#' Reads the OrthoFinder gene count table, computes a binary distance matrix
#' across elements (columns), and clusters them using average-linkage
#' hierarchical clustering. A fixed height cutoff (h = 0.99999999) is used to
#' identify top-level clusters.
#'
#' @param orthogroups_file Character. Path to the OrthoFinder gene count TSV file.
#'
#' @return A named list with elements:
#'   \item{transposed_counts}{Numeric matrix of orthogroup counts, elements as rows}
#'   \item{orthofinder_counts}{Numeric matrix without the "Total" column}
#'   \item{distance_mat}{dist object of binary distances between elements}
#'   \item{hierar_cl}{hclust object from average-linkage clustering}
#'   \item{cluster_fit}{Named integer vector of cluster assignments per element}
#'   \item{individual_clusters}{Integer. Number of top-level clusters identified}
perform_initial_clustering <- function(orthogroups_file = "Orthogroups.GeneCount.tsv") {

  ortho_raw <- read.table(orthogroups_file, header = TRUE, check.names = FALSE, row.names = 1)

  # Remove the "Total" summary column before transposing
  ortho_counts  <- ortho_raw[, !names(ortho_raw) %in% "Total"]
  transposed    <- t(ortho_counts)

  distance_mat  <- dist(transposed, method = "binary")
  hierar_cl     <- hclust(distance_mat, method = "average")

  # Ensure monotone height vector 
  hierar_cl$height <- sort(hierar_cl$height)

  cluster_fit        <- cutree(hierar_cl, h = 0.99999999)
  individual_clusters <- length(unique(cluster_fit))

  return(list(
    transposed_counts  = transposed,
    orthofinder_counts = ortho_counts,
    distance_mat       = distance_mat,
    hierar_cl          = hierar_cl,
    cluster_fit        = cluster_fit,
    individual_clusters = individual_clusters
  ))
}


#' Identify shared accessory genes
#'
#' Defines shared accesory genes (genes shared between subclusters) and subcluster accesory 
#' genes (only present in one subcluster). Results are written to disk as both .txt and .csv files.
#'
#' @param accessory        Numeric matrix. Orthogroup counts (orthogroups x elements).
#' @param accessory_matrix List of the Numeric matrix of accessory genes in each subcluster.
#' @param orthocounts      Numeric matrix. Orthogroup counts (orthogroups x elements).
#' @param cluster_number   Character. Suffix appended to output filenames. Use ""
#'                         for the top-level (non-clustered) analysis.
#'
#' @return A named list with elements:
#'   \item{shared_acc}{Character vector of share accessory orthogroup IDs}
#'   \item{subcluster_acc}{Character vector of subcluster accesory orthogroup IDs}
Shared_accessory_analysis <- function(accessory, accessory_matrix, orthocounts, cluster_number = "") {
  
  shared_acc <- character()
  sub_acc <- character()
  for (sub_id in names(accessory)) {
    subcluster_acc <- accessory[[paste(sub_id)]]
    subcluster_acc <- subcluster_acc[! subcluster_acc %in% shared_acc]
    matrix <- accessory_matrix[[paste(sub_id)]]
    elements_cluster <- colnames(matrix)
    for(orthogroupID in subcluster_acc) {
      Orthogroup <- as.data.frame(orthocounts[orthogroupID, ! colnames(orthocounts) %in% elements_cluster])
      Orthogroup <- Orthogroup[,colSums(Orthogroup) >= 1]
      if(is.null(ncol(Orthogroup))){ next }
      if(ncol(Orthogroup) >= 1) {
        shared_acc <- c(shared_acc, orthogroupID)
        subcluster_acc <- subcluster_acc[ ! subcluster_acc %in% orthogroupID]
      }
    }

    if(length(subcluster_acc) >= 1) {
      sub_acc <- c(sub_acc, subcluster_acc)
      subcluster_acc_matrix <- orthocounts[subcluster_acc, elements_cluster]
      write.table(
        subcluster_acc,
        file      = paste0("SubCluster", sub_id, cluster_number, "-Accessory", ".txt"),
        sep       = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
      )
      write.csv(
        subcluster_acc_matrix,
        file      = paste0("SubCluster", sub_id, cluster_number, "-Accessory",".csv"),
        row.names = TRUE, quote = FALSE
      )
    }
  }

  if( length(shared_acc) >= 1){
    share_acc_matrix <- orthocounts[shared_acc,]
    share_acc_matrix <- share_acc_matrix[rowSums(share_acc_matrix) >= 1,]
    write.table(
      shared_acc,
      file      = "Shared-Accessory.txt",
      sep       = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    write.csv(
      share_acc_matrix,
      file      = paste0("Shared-Accessory", ".csv"),
      row.names = TRUE, quote = FALSE
    )
  }

  return(list(shared_acc = shared_acc, subcluster_acc = sub_acc))
}


#' Identify core and accessory genes
#'
#' Defines core genes (present in >= 80% of all elements), shared accesory genes (genes shared 
#' between subclusters) and subcluster accesory genes (present in >= 80% of elements within
#' each subcluster). Results are written to disk as both .txt and .csv files.
#'
#' @param ortho_counts   Numeric matrix. Orthogroup counts (orthogroups x elements).
#' @param cluster        hclust object. Hierarchical clustering of the current set.
#' @param cluster_number Character. Suffix appended to output filenames. Use ""
#'                       for the top-level (non-clustered) analysis.
#' @param subcluster     data.frame or vector. Subcluster definitions from
#'                       load_and_preprocess_data()$subclusters.
#'
#' @return A named list with elements:
#'   \item{core}{Character vector of general core orthogroup IDs}
#'   \item{shared_acc}{Character vector of share accessory orthogroup IDs}
#'   \item{subcluster_acc}{Character vector of subcluster accesory orthogroup IDs}
pangenome_analysis <- function(ortho_counts, cluster, cluster_number = "", subcluster) {

  # --- General core genes (present in >= 80% of all elements) ----------------
  core_mat <- ortho_counts[rowSums(ortho_counts < 1) <= ncol(ortho_counts) * 0.2, ]
  core     <- rownames(core_mat)

  if (length(core) > 0) {
    write.table(
      core,
      file      = paste0("Core_genes", cluster_number, ".txt"),
      sep       = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )
    write.csv(
      core_mat,
      file      = paste0("Core_genes", cluster_number, ".csv"),
      row.names = TRUE, quote = FALSE
    )
  }

  subcluster_acc <- list()
  subcluster_acc_mat <- list()

  # --- Subcluster-specific core genes -----------------------------------------
  # Branch A: clustering is performed internally (cutree at h = 0.8)
  if (cluster_number != "" || arguments$subclusters == 1) {

    fit <- cutree(cluster, h = 0.8)

    if (length(unique(fit)) > 1) {
      for (sub_id in seq_along(unique(fit))) {
        sub_counts <- subset(ortho_counts, select = names(fit[fit == sub_id]))

        if (ncol(sub_counts) > 4) {
          sub_acc_mat <- sub_counts[rowSums(sub_counts < 1) <= ncol(sub_counts) * 0.2, ]
          sub_acc_mat <- sub_acc_mat[!(rownames(sub_acc_mat) %in% core), ]

          if (nrow(sub_acc_mat) > 0) {
            accesory         <- rownames(sub_acc_mat)
            subcluster_acc[[paste(sub_id)]] <- accesory
            subcluster_acc_mat[[paste(sub_id)]] <- sub_acc_mat
          }
        }
      }
    }

    if(length(subcluster_acc) >= 1) {
      accesory_definition <- Shared_accessory_analysis(subcluster_acc,subcluster_acc_mat, ortho_counts)
      shared_acc <- accesory_definition$shared_acc
      subcluster_acc <- accesory_definition$subcluster_acc
    } else {
      shared_acc <- character()
      subcluster_acc <- character()
    }

  # Branch B: subcluster definitions supplied externally via MCL file
  } else if (arguments$subclusters > 1) {

    for (sub_id in seq_along(subcluster$Subcluster)) {
      sub_elements <- unlist(strsplit(subcluster$Elements[sub_id], split = " "))
      sub_counts   <- subset(ortho_counts, select = sub_elements)

      if (ncol(sub_counts) > 4) {
        sub_acc_mat <- sub_counts[rowSums(sub_counts < 1) <= ncol(sub_counts) * 0.2, ]
        sub_acc_mat <- sub_acc_mat[!(rownames(sub_acc_mat) %in% core), ]

        if (nrow(sub_acc_mat) > 0) {
          accesory          <- rownames(sub_acc_mat)
          subcluster_acc[[paste(sub_id)]] <- accesory
          subcluster_acc_mat[[paste(sub_id)]] <- sub_acc_mat
        }
      }
    }

    if(length(subcluster_acc) >= 1) {
      accesory_definition <- Shared_accessory_analysis(subcluster_acc,subcluster_acc_mat, ortho_counts)
      shared_acc <- accesory_definition$shared_acc
      subcluster_acc <- accesory_definition$subcluster_acc
    } else {
      shared_acc <- character()
      subcluster_acc <- character()
    }
  }

  return(list(core = core, shared_acc = shared_acc, subcluster_acc = subcluster_acc))
}


#' Analyze and visualize individual top-level clusters
#'
#' For each top-level cluster with more than 2 orthogroups, builds a
#' phylogenetic tree and a heatmap of the orthogroup presence/absence profile.
#' Outputs are written as .nwk and .svg files.
#'
#' @param individual_clusters  Integer. Total number of top-level clusters.
#' @param transposed_counts    Numeric matrix. Orthogroup counts with elements as rows.
#' @param cluster_fit          Named integer vector. Cluster assignments per element.
#' @param seq_data             data.frame. Sequence lengths (unused here; reserved for
#'                             future extension).
#' @param gene_data            data.frame. Gene coordinates (unused here; reserved for
#'                             future extension).
#' @param blast_links          data.frame. BLAST results (unused here; reserved for
#'                             future extension).
#'
#' @return Numeric vector of cluster IDs that were successfully analyzed.
analyze_individual_clusters <- function(
  individual_clusters,
  transposed_counts,
  cluster_fit,
  seq_data,
  gene_data,
  blast_links
) {

  selected_clusters <- numeric()

  for (cluster_id in seq_len(individual_clusters)) {

    cluster_matrix <- transposed_counts[cluster_fit == cluster_id, ]

    if (is.matrix(cluster_matrix)) {
      # Remove orthogroups absent in all elements of this cluster
      cluster_matrix <- cluster_matrix[, colSums(abs(cluster_matrix)) > 0]
    } else {
      # Single-element cluster: create a placeholder 1x1 matrix
      cluster_matrix <- matrix(nrow = 1, ncol = 1)
      rownames(cluster_matrix) <- names(cluster_fit[cluster_fit == cluster_id])
    }

    if (ncol(cluster_matrix) > 2) {

      cluster_dist <- dist(cluster_matrix, method = "binary")
      cluster_hier <- hclust(cluster_dist, method = "average")
      my_tree      <- as.phylo(cluster_hier)

      # --- Write tree and element list ----------------------------------------
      if (individual_clusters >= 2) {
        selected_clusters <- c(selected_clusters, cluster_id)
        write.tree(my_tree, file = paste0("CargoHierarchicalTree_Cluster", cluster_id, ".nwk"))
        write.table(
          cluster_hier$labels,
          file      = paste0("Cluster", cluster_id, ".txt"),
          sep       = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
        )
      } else {
        write.tree(my_tree, file = "CargoHierarchicalTree.nwk")
      }

      # --- Heatmap of orthogroup presence/absence profile ---------------------
      cluster_melt <- reshape2::melt(cluster_matrix, as.is = TRUE)

      profile_plot <- ggplot(data = cluster_melt, aes(x = reorder(Var2, value), y = Var1, fill = value)) +
        geom_tile() +
        scale_fill_viridis(option = "rocket", direction = -1) +
        labs(y = "Element", x = "Orthogroup") +
        theme(axis.text.x = element_blank())

      # Cap dimensions at 49 inches (svglite limitation)
      plot_width  <- min(49, ncol(cluster_matrix) * 0.2 + 2)
      plot_height <- min(49, nrow(cluster_matrix))

      if (individual_clusters >= 2) {
        ggsave(profile_plot, filename = paste0("CargoHeatmap_Cluster", cluster_id, ".svg"),
               width = plot_width, height = plot_height)
      } else {
        ggsave(profile_plot, filename = "CargoHeatmap.svg",
               width = plot_width, height = plot_height)
      }

    } else {
      log_message(paste0("WARNING: Cluster ", cluster_id, " has too few orthogroups to analyze. ",
                         "Affected element(s): ", paste(rownames(cluster_matrix), collapse = ", ")))
    }
  }

  return(selected_clusters)
}

#' Annotate a gene data.frame with pangenome gene categories
#'
#' Helper used by plot_subcluster_synteny() to avoid repeating the same
#' three-step annotation logic in multiple code paths.
#'
#' @param genes_df         data.frame. Gene coordinates with an 'attribute' column.
#' @param orthogroup_table data.frame. Columns: Orthogroup, genes (space-delimited).
#' @param pangenome       Named list. Output of pangenome_analysis().
#'
#' @return data.frame. The input genes_df with the 'attribute' column updated.
annotate_gene_categories <- function(genes_df, orthogroup_table, pangenome) {

  extract_gene_ids <- function(ortho_ids) {
    tbl <- orthogroup_table[orthogroup_table$Orthogroup %in% ortho_ids, ]
    ids <- unlist(strsplit(tbl[, 2], split = " "))
    ids[ids != ""]
  }

  if (length(pangenome$shared_acc) > 0) {
    shared_genes <- extract_gene_ids(pangenome$shared_acc)
    genes_df   <- genes_df %>%
      mutate(attribute = ifelse(ID %in% shared_genes, "Shared Accessory", attribute))
  }

  if (length(pangenome$core) > 0) {
    core_genes <- extract_gene_ids(pangenome$core)
    genes_df  <- genes_df %>%
      mutate(attribute = ifelse(ID %in% core_genes, "Core", attribute))
  }

  if (length(pangenome$subcluster_acc) > 0) {
    sub_genes <- extract_gene_ids(pangenome$subcluster_acc)
    genes_df  <- genes_df %>%
      mutate(attribute = ifelse(ID %in% sub_genes, "Subcluster Accessory", attribute))
  }

  return(genes_df)
}

#' Annotate gene features and plot a synteny map for a set of elements
#'
#' Builds a ladderized phylogenetic tree from the cluster's distance matrix,
#' annotates genes as general core, specific core, or unannotated, and
#' combines the tree with a gggenomes synteny panel into a single composite plot.
#' The output SVG is named according to the provided cluster/subcluster identifiers.
#'
#' @param cluster_matrix    Numeric matrix. Orthogroup presence/absence with elements as rows.
#' @param orthogroup_table  data.frame. Columns: Orthogroup, genes (space-delimited gene IDs).
#' @param pangenome        Named list. Output of pangenome_analysis()
#' @param seq_data          data.frame. Sequence lengths (seq_id, length).
#' @param gene_data         data.frame. Gene coordinates (seq_id, start, end, strand, type, attribute, ID).
#' @param blast_links       data.frame. BLAST links reformatted for gggenomes (seq_id, start, end,
#'                          seq_id2, start2, end2, pident).
#' @param cluster_number    Character. Appended to output filename for multi-cluster runs.
#'                          Use "" for a single-cluster run.
#' @param subcluster_number Character. Appended to output filename for subcluster runs.
#'                          Use "" when not plotting a subcluster.
#'
#' @return NULL (called for side effects: writes SVG to disk)
plot_cluster_synteny <- function(
  cluster_matrix,
  orthogroup_table,
  pangenome,
  seq_data,
  gene_data,
  blast_links,
  cluster_number    = "",
  subcluster_number = ""
) {

  # --- Build ladderized tree from cluster distance matrix ---------------------
  cluster_dist <- dist(cluster_matrix, method = "binary")
  cluster_hier <- hclust(cluster_dist, method = "average")
  my_tree      <- ladderize(as.phylo(cluster_hier), right = FALSE)

  selected_seqs <- rownames(cluster_matrix)
  seqs_filtered <- seq_data  %>% filter(seq_id %in% selected_seqs)
  genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs)

  genes_filtered <- annotate_gene_categories(genes_filtered, orthogroup_table, pangenome)

  links_filtered <- blast_links %>%
    filter(seq_id %in% selected_seqs | seq_id2 %in% selected_seqs)


  # --- Order sequences to match tree tip order --------------------------------
  tree_base    <- ggtree::ggtree(my_tree, layout = "rectangular")
  tree_y_order <- tree_base$data %>%
    filter(isTip) %>%
    select(seq_id = label, y) %>%
    mutate(y = max(y) - y + 1)

  ordered_seqs <- seqs_filtered %>%
    inner_join(tree_y_order, by = "seq_id") %>%
    arrange(y)

  max_seq_len <- max(ordered_seqs$length, na.rm = TRUE)
  n_elements  <- nrow(cluster_matrix)

  # --- Tree panel -------------------------------------------------------------
  p_tree <- tree_base +
    ggtree::geom_tiplab(
      align  = TRUE,
      size   = ifelse(n_elements < 100, 2.5, 3),
      family = "mono"
    ) +
    ggtree::theme_tree2()

  # --- Genome / synteny panel -------------------------------------------------
  p_genome <- suppressMessages(
    gggenomes(seqs = ordered_seqs, links = links_filtered, genes = genes_filtered) +
      geom_seq(aes(y = y)) +
      geom_gene(aes(y = y, fill = attribute), show.legend = TRUE) +
      scale_fill_manual(
        name   = "Gene category",
        values = c("Core" = "red4", "Shared Accessory" = "green4", "Subcluster Accessory" = "darkorchid"),
        na.value = "cornsilk3",
        limits = c("Core", "Shared Accessory", "Subcluster Accessory")
      ) +
      new_scale_fill() +
      geom_link(aes(y = y, fill = pident), colour = NA) +
      scale_fill_continuous(name = "Alignment Identity (%)") +
      geom_bin_label(aes(y = y), x = -10, size = ifelse(n_elements < 50, 2.5, 3)) +
      scale_x_continuous(
        labels = label_number(accuracy = 1),
        limits = c(0, max_seq_len)
      ) +
      scale_y_continuous(
        expand = expansion(mult = 1.0001 * (n_elements ^ -1.1728))
      )
  )

  final_plot <- p_tree + p_genome

  # Cap dimensions at 49 inches (svglite limitation)
  plot_width  <- min(49, max(16, round(max_seq_len * 0.0001) +
                                  round(max(my_tree$edge.length, na.rm = TRUE) * 10)))
  plot_height <- min(49, n_elements)

  if (cluster_number != "") {
    if(subcluster_number != "") {
      filename <- paste0("CargoSynteny_Cluster", cluster_number, "-SubCluster", subcluster_number, ".svg")
    } else {
      filename <- paste0("CargoSynteny_Cluster", cluster_number, ".svg")
    }
  } else if (subcluster_number != "") {
    filename <- paste0("CargoSynteny_SubCluster", subcluster_number, ".svg")
  } else {
    filename <- "CargoSynteny.svg"
  }

  ggsave(final_plot, filename = filename, width = plot_width, height = plot_height, limitsize = FALSE)
}


#' Detect nesting events between starship elements
#'
#' For each element pair with shared orthogroups, tests whether one element
#' is fully contained within another based on BLAST alignment coverage.
#' A nesting event is flagged when the query element's alignments span
#' >= 80% of its own length, while the subject element's alignments
#' do not start near its own left boundary. Detected events are
#' visualized and saved as individual SVG files.
#'
#' @param ortho_counts Numeric matrix. Orthogroup counts (orthogroups x elements).
#' @param seq_data     data.frame. Sequence lengths (seq_id, length).
#' @param gene_data    data.frame. Gene coordinates.
#' @param blast_links  data.frame. Full BLAST results (qseqid, sseqid, qstart, qend,
#'                     sstart, send, pident, qlen, slen).
#'
#' @return NULL
check_nesting <- function(ortho_counts, seq_data, gene_data, blast_links) {

  for (elem_idx in seq_len(ncol(ortho_counts))) {

    elem_name  <- colnames(ortho_counts)[elem_idx]

    # Orthogroups present in this element
    elem_orthogroups <- rownames(
      ortho_counts[elem_name] %>% filter(!if_all(everything(), ~ .x == 0))
    )

    # Subset to orthogroups shared with any other element (>= 80% coverage threshold)
    reduced_ortho  <- ortho_counts[elem_orthogroups, ]
    reduced_shared <- reduced_ortho[
      , colSums(reduced_ortho < 1) < nrow(reduced_ortho) * 0.2,
      drop = FALSE
    ]

    if (!is.data.frame(reduced_shared)) next

    # Remove the focal element itself before comparing
    reduced_shared <- reduced_shared %>% select(-all_of(elem_name))
    if (!is.data.frame(reduced_shared)) next

    for (comp_idx in seq_len(ncol(reduced_shared))) {

      comp_name      <- colnames(reduced_shared)[comp_idx]
      comp_orthogroups <- rownames(
        ortho_counts[comp_name] %>% filter(!if_all(everything(), ~ .x == 0))
      )

      # Only test if the comparison element has substantially more orthogroups
      # (>= 1.5x), which is a prerequisite for the nested-inside relationship
      if (length(comp_orthogroups) <= length(elem_orthogroups) * 1.5) next

      high_id_links <- blast_links %>%
        filter(qseqid %in% elem_name, sseqid %in% comp_name) %>%
        select(
          seq_id = qseqid, start = qstart, end = qend,
          seq_id2 = sseqid, start2 = sstart, end2 = send, pident
        )

      if (nrow(high_id_links) == 0) next

      elem_len <- as.numeric(seq_data[seq_data$seq_id == unique(high_id_links$seq_id),  2])
      comp_len <- as.numeric(seq_data[seq_data$seq_id == unique(high_id_links$seq_id2), 2])

      # Nesting criteria:
      #   - Query (elem) alignments span at least 80% of its own length
      #   - Subject (comp) alignments do NOT start near its left boundary
      is_nested     <- any(high_id_links$start  < elem_len * 0.2) &&
                       any(high_id_links$end    > elem_len * 0.8)
      starts_inside <- any(high_id_links$start2 < comp_len * 0.2) | any(high_id_links$start2 > comp_len * 0.8)

      if (!is_nested || starts_inside) next

      log_message(paste0("Element '", elem_name, "' is nested inside element '", comp_name, "'."))

      selected_seqs  <- c(elem_name, comp_name)
      seqs_filtered  <- seq_data  %>% filter(seq_id %in% selected_seqs)
      genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs)

      nest_links <- blast_links %>%
        filter(qseqid %in% selected_seqs, sseqid %in% selected_seqs) %>%
        select(
          seq_id = qseqid, start = qstart, end = qend,
          seq_id2 = sseqid, start2 = sstart, end2 = send, pident
        )

      ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, selected_seqs))

      p_nest <- gggenomes(seqs = ordered_seqs, links = nest_links, genes = genes_filtered) +
        geom_seq(aes(y = y)) +
        geom_gene(aes(y = y)) +
        geom_link(aes(y = y, fill = pident), colour = NA) +
        geom_bin_label(aes(y = y), x = -10) +
        scale_x_continuous(
          labels = label_number(accuracy = 1),
          limits = c(0, max(ordered_seqs$length))
        ) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))

      ggsave(
        p_nest,
        filename  = paste0("IndividualNestingEvent-", elem_name, "in", comp_name, ".svg"),
        width     = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)),
        height    = min(49, length(selected_seqs)),
        limitsize = FALSE
      )
    }
  }
}


#' Visualize synteny per subcluster
#'
#' For each subcluster, plots the main synteny view.
#'
#' @param subcluster_number Integer. Total number of subclusters to iterate over.
#' @param orthogroup_table  data.frame. Orthogroup gene-list table.
#' @param dist_matrix       Numeric matrix. Full binary distance matrix (elements x orthogroups).
#' @param cluster_fit       Named integer vector. Subcluster assignments per element.
#' @param ortho_counts      Numeric matrix. Full orthogroup count matrix.
#' @param pangenome        Named list. Output of pangenome_analysis().
#' @param seq_data          data.frame. Sequence lengths.
#' @param gene_data         data.frame. Gene coordinates.
#' @param blast_links       data.frame. Full BLAST results.
#' @param cluster_matrix    Numeric matrix. Full transposed orthogroup count matrix.
#' @param cluster_number    Character. Prefix for output file names.
#' @param subcluster        data.frame or vector. External subcluster definitions.
#'
#' @return NULL
plot_subcluster_synteny <- function(
  subcluster_number,
  orthogroup_table,
  dist_matrix,
  cluster_fit,
  ortho_counts,
  pangenome,
  seq_data,
  gene_data,
  blast_links,
  cluster_matrix,
  cluster_number = "",
  subcluster
) {

  for (sub_id in seq_len(subcluster_number)) {

    if (cluster_number != "") {
      matrix           <- dist_matrix[cluster_fit == sub_id, , drop = FALSE]
      cluster_elements <- rownames(matrix)
      if (is.null(cluster_elements)) {
        cluster_elements <- names(cluster_fit[cluster_fit == sub_id])
      }
    } else {
      cluster_elements <- unlist(strsplit(subcluster$Elements[sub_id], split = " "))
      matrix           <- dist_matrix[cluster_elements, , drop = FALSE]
    }

    write.table(
      cluster_elements,
      file      = paste0(cluster_number, "SubCluster", sub_id, ".txt"),
      sep       = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
    )

    if (!is.matrix(matrix)) {
      log_message(paste0("Subcluster ", sub_id, " cannot be analyzed automatically (single element or malformed matrix)."))
      next
    }

    # --- Main synteny plot for this subcluster --------------------------------
    sub_cluster_matrix <- cluster_matrix[rownames(cluster_matrix) %in% cluster_elements, ]

    links_reformatted <- blast_links %>%
      select(
        seq_id = qseqid, start = qstart, end = qend,
        seq_id2 = sseqid, start2 = sstart, end2 = send, pident
      )

    plot_cluster_synteny(
      cluster_matrix    = sub_cluster_matrix,
      orthogroup_table  = orthogroup_table,
      pangenome        = pangenome,
      seq_data          = seq_data,
      gene_data         = gene_data,
      blast_links       = links_reformatted,
      subcluster_number = sub_id
    )
  }
}

#' Identify discordant elements between two cluster assignment vectors
#'
#' Aligns cluster labels between two partitions (to account for arbitrary label
#' permutations) and returns the names of elements whose assignments differ.
#'
#' @param tree1 Named integer vector. Cluster assignments from partition 1.
#' @param tree2 Named integer vector. Cluster assignments from partition 2.
#'
#' @return Character vector of element names with discordant assignments.
get_discordant <- function(tree1, tree2) {

  # Build a contingency table and map each label in tree1 to the most frequent
  # corresponding label in tree2 (greedy label alignment)
  ct  <- table(tree1, tree2)
  best_match <- apply(ct, 1, which.max)

  aligned_tree1 <- tree1
  for (i in seq_along(best_match)) {
    aligned_tree1[tree1 == names(best_match)[i]] <- as.numeric(best_match[i])
  }

  mismatched_indices <- which(tree2 != aligned_tree1)
  discordant_names <- names(tree2)[mismatched_indices]
  return(discordant_names)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

log_message("Reading and preprocessing data...")
data_list <- load_and_preprocess_data(min_pident = arguments$pident)

log_message("Performing initial hierarchical clustering...")
clustering_data <- perform_initial_clustering()

log_message("Checking for nesting events...")
check_nesting(
  ortho_counts = clustering_data$orthofinder_counts,
  seq_data     = data_list$seqs,
  gene_data    = data_list$genes,
  blast_links  = data_list$blast_results
)

log_message(paste0("Analyzing ", clustering_data$individual_clusters, " top-level cluster(s)..."))

# ==============================================================================
# BRANCH: SINGLE TOP-LEVEL CLUSTER
# ==============================================================================

if (clustering_data$individual_clusters == 1) {

  empty <- analyze_individual_clusters(
  individual_clusters = clustering_data$individual_clusters,
  transposed_counts   = clustering_data$transposed_counts,
  cluster_fit         = clustering_data$cluster_fit,
  seq_data            = data_list$seqs,
  gene_data           = data_list$genes,
  blast_links         = data_list$blast_results
  )

  if (arguments$subclusters >= 2) {

    log_message("Analyzing pangenome...")
    pangenome <- pangenome_analysis(
      ortho_counts   = clustering_data$orthofinder_counts,
      cluster        = clustering_data$hierar_cl,
      subcluster     = data_list$subclusters
    )

    log_message("Generating main synteny plot...")
    plot_cluster_synteny(
      cluster_matrix    = clustering_data$transposed_counts,
      orthogroup_table  = data_list$orthogroups_table,
      pangenome        = pangenome,
      seq_data          = data_list$seqs,
      gene_data         = data_list$genes,
      blast_links       = data_list$links
    )

    if (arguments$captainRemoval <= 1) {
      log_message(paste0("Analyzing ", arguments$subclusters, " subcluster(s) (user-defined)..."))
      fit <- cutree(clustering_data$hierar_cl, k = arguments$subclusters)

    } else {
      log_message("Optimizing subcluster number with NbClust...")
      nb_result <- NbClust(
        clustering_data$distance_mat,
        method = "average",
        min.nc = 1,
        max.nc = min(6, nrow(clustering_data$distance_mat) - 1),
        index  = "ball"
      )
      k_number <- round(mean(c(
        length(unique(nb_result$Best.partition)),
        arguments$subclusters
      )))
      log_message(paste0("Optimal K = ", k_number, ". Analyzing subclusters..."))
      fit <- cutree(clustering_data$hierar_cl, k = k_number)
      arguments$subclusters <- k_number
    }

    plot_subcluster_synteny(
      subcluster_number = arguments$subclusters,
      orthogroup_table  = data_list$orthogroups_table,
      dist_matrix       = as.matrix(clustering_data$distance_mat),
      cluster_fit       = fit,
      ortho_counts      = clustering_data$orthofinder_counts,
      pangenome        = pangenome,
      seq_data          = data_list$seqs,
      gene_data         = data_list$genes,
      blast_links       = data_list$blast_results,
      subcluster        = data_list$subclusters,
      cluster_matrix    = clustering_data$transposed_counts
    )

  } else {

    log_message("Generating single-cluster synteny plot...")
    pangenome <- pangenome_analysis(
      ortho_counts = clustering_data$orthofinder_counts,
      cluster      = clustering_data$hierar_cl,
      subcluster   = vector()
    )
    plot_cluster_synteny(
      cluster_matrix    = clustering_data$transposed_counts,
      orthogroup_table  = data_list$orthogroups_table,
      pangenome        = pangenome,
      seq_data          = data_list$seqs,
      gene_data         = data_list$genes,
      blast_links       = data_list$links
    )
  }

# ==============================================================================
# BRANCH: MULTIPLE TOP-LEVEL CLUSTERS
# ==============================================================================

} else if (clustering_data$individual_clusters > 1) {

  if (arguments$subclusters >= 2) {
    selected_cluster_ids <- analyze_individual_clusters(
      individual_clusters = clustering_data$individual_clusters,
      transposed_counts   = clustering_data$transposed_counts,
      cluster_fit         = clustering_data$cluster_fit,
      seq_data            = data_list$seqs,
      gene_data           = data_list$genes,
      blast_links         = data_list$blast_results
    )

    for (main_cluster_id in selected_cluster_ids) {

      cluster_matrix <- clustering_data$transposed_counts[
        clustering_data$cluster_fit == main_cluster_id, , drop = FALSE
      ]

      if (!is.matrix(cluster_matrix) || nrow(cluster_matrix) < 4) {
        log_message(paste0("Cluster ", main_cluster_id, " cannot be analyzed automatically (too few elements)."))
        next
      }

      # Remove orthogroups absent in all elements of this cluster
      cluster_matrix <- cluster_matrix[, colSums(abs(cluster_matrix)) > 0]

      log_message(paste0("Analyzing core genes for Cluster ", main_cluster_id, "..."))

      new_ortho_counts <- clustering_data$orthofinder_counts[
        , clustering_data$cluster_fit == main_cluster_id
      ]
      cluster_dist <- dist(cluster_matrix, method = "binary")
      cluster_hier <- hclust(cluster_dist, method = "average")

      cluster_suffix <- paste0("_Cluster", main_cluster_id)
      pangenome     <- pangenome_analysis(
        ortho_counts   = new_ortho_counts,
        cluster        = cluster_hier,
        cluster_number = cluster_suffix,
        subcluster     = vector()
      )

      plot_cluster_synteny(
        cluster_matrix    = cluster_matrix,
        orthogroup_table  = data_list$orthogroups_table,
        pangenome        = pangenome,
        seq_data          = data_list$seqs,
        gene_data         = data_list$genes,
        blast_links       = data_list$links,
        cluster_number    = main_cluster_id
      )

      nb_result <- NbClust(
        cluster_dist,
        method = "average",
        min.nc = 1,
        max.nc = min(6, nrow(cluster_matrix) - 1),
        index  = "ball"
      )
      k_number <- round(mean(c(
        length(unique(nb_result$Best.partition)),
        length(unique(cutree(cluster_hier, h = 0.6)))
      )))
      fit <- cutree(cluster_hier, k = k_number)

      plot_subcluster_synteny(
        subcluster_number = k_number,
        orthogroup_table  = data_list$orthogroups_table,
        dist_matrix       = as.matrix(cluster_dist),
        cluster_fit       = fit,
        ortho_counts      = clustering_data$orthofinder_counts,
        pangenome        = pangenome,
        seq_data          = data_list$seqs,
        gene_data         = data_list$genes,
        blast_links       = data_list$blast_results,
        cluster_matrix    = clustering_data$transposed_counts,
        cluster_number    = paste0("Cluster", main_cluster_id, "-"),
        subcluster        = vector()
      )
    }
  }

} else {
  log_message("No clusters to analyze.")
}

# ==============================================================================
# CAPTAIN VS CARGO DISCORDANCE ANALYSIS
# ==============================================================================

if (clustering_data$individual_clusters == 1) {

  log_message("Analyzing Captain vs Cargo discordance...")

  if (!file.exists("CaptainPhylogeny.nw")) {
    log_message("'CaptainPhylogeny.nw' not found. Skipping discordance analysis.")

  } else {

    phylo_tree   <- read.tree("CaptainPhylogeny.nw")
    dist_matrix  <- cophenetic.phylo(phylo_tree)

    # Normalize to [0, 1] before clustering to make the tree distances
    # comparable to the binary orthogroup distances used for the cargo tree
    norm_dist    <- (dist_matrix - min(dist_matrix)) / (max(dist_matrix) - min(dist_matrix))
    hclust_tree  <- hclust(as.dist(norm_dist), method = "average")
    hclust_tree$height <- sort(hclust_tree$height)

    # Determine optimal K: use NbClust for trees with > 4 tips, otherwise default to K = 2
    if (length(phylo_tree$tip.label) > 4) {
      nb_tree  <- NbClust(
        scale(norm_dist),
        method = "average",
        min.nc = 1,
        max.nc = min(6, nrow(norm_dist) - 1),
        index  = "ball"
      )
      k_num <- round(mean(c(
        length(unique(nb_tree$Best.partition)),
        length(unique(cutree(hclust_tree, h = 0.5)))
      )))
    } else {
      k_num <- 2
    }

    fit_tree  <- cutree(hclust_tree, k = k_num)
    fit_tree <- fit_tree[sort(names(fit_tree))]

    fit_cargo <- cutree(clustering_data$hierar_cl, k = k_num)
    fit_cargo <- fit_cargo[sort(names(fit_cargo))]

    discordant <- c(
      get_discordant(fit_tree,  fit_cargo),
      get_discordant(fit_cargo, fit_tree)
    )
    discordant <- discordant[ave(discordant, discordant, FUN = length) == 1]

    if (length(discordant) >= 1) {
      log_message(paste0(length(discordant), " discordant element(s) identified."))
      write.table(
        discordant,
        file      = "Discordant_elements.txt",
        sep       = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE
      )
    } else {
      log_message("No discordance detected between Captain and Cargo trees.")
    }

    # --- Tanglegram -----------------------------------------------------------
    svg(
      filename = "CaptainVsCargo_tanglegram.svg",
      width    = 16,
      height   = min(49, length(fit_cargo) / 2)
    )
    dendlist(clustering_data$hierar_cl, hclust_tree) %>%
      ladderize() %>%
      untangle(method = "step1side", k_seq = k_num:(length(fit_cargo) - 1)) %>%
      set("branches_k_color", k = k_num) %>%
      tanglegram(
        faster     = TRUE,
        main_left  = "Cargo Cluster",
        main_right = "Captain Tree"
      )
    invisible(dev.off())
  }
}

log_message("Done.")
