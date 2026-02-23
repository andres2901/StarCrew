# ==============================================================================
# CLUSTER ANALYSIS SCRIPT
# ==============================================================================
# This script performs comparative genomic analysis, including hierarchical 
# clustering of orthogroups, synteny visualization, and nesting detection.
# ==============================================================================

# Argument parsing and validation
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory", metavar="PATH"),
  make_option(c("-s", "--subclusters"), type="integer", action = "store", default=1,
              help=" Number of subclusters get in the synteny analysis. Minimum value is equal to 1 when no subcluster were identified [default %default]", metavar="number"),
  make_option(c("-c", "--captainRemoval"), type="integer", action = "store", default=0,
              help=" Number of elements removes from the original cluster [default %default]", metavar="number"),
  make_option(c("-p", "--pident"), type="numeric", action = "store", default=70,
              help=" Minimum percentage identity of blast results to be included as links for nucleotide synteny visualization [default %default]", metavar="number")
)

arguments <- parse_args(OptionParser(option_list = option_list))

if (is.null(arguments$directory)) {
  stop("Error: directory must be provided.", call. = FALSE)
} else {
  setwd(arguments$directory)
}

# Software Dependency Checks
suppressPackageStartupMessages({
  library(ggplot2)
  library(ape)
  library(ggtree)
  library(reshape2)
  library(viridis)
  library(dplyr)
  library(gggenomes)
  library(scales)
  library(dendextend)
  library(NbClust)
  library(syntenet)
  library(ggnewscale)
})

required_pkgs <- c("ggplot2", "ape", "ggtree", "reshape2", "viridis", "dplyr", 
                   "gggenomes", "scales", "dendextend", "NbClust", "syntenet", "ggnewscale")

for (pkg in required_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Package \"", pkg, "\" not installed. Please install it to run this script."), call. = FALSE)
  }
}

# ==============================================================================
# FUNCTION DEFINITIONS
# ==============================================================================

load_and_preprocess_data <- function(blast_file = "Blast_CleanResults.txt", gff_file = "Final_model.gff", min_pident) {
  blast_col_names <- c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "length", "qlen", "slen")
  blast_col_classes <- c("character", "character", rep("numeric", 8))
  
  blast_results <- read.delim(blast_file, header = FALSE, sep = "\t", col.names = blast_col_names, colClasses = blast_col_classes) %>%
    filter(qseqid != sseqid, pident > min_pident)
  
  links <- blast_results %>%
    mutate(alignment_id = row_number()) %>%
    select(seq_id = qseqid, start = qstart, end = qend, seq_id2 = sseqid, start2 = sstart, end2 = send, pident, length)
  
  seqs <- bind_rows(
    blast_results %>% select(seq_id = qseqid, length = qlen),
    blast_results %>% select(seq_id = sseqid, length = slen)
  ) %>%
    distinct() %>%
    group_by(seq_id) %>%
    summarize(length = max(length)) %>%
    ungroup()
  
  annotation <- gff2GRangesList("Gff/")
  genes <- as.data.frame(unlist(annotation)) %>%
    filter(type == "gene") %>%
    mutate(seqnames = as.character(seqnames), type = "CDS", attribute = NA) %>%
    select(seq_id = seqnames, start, end, strand, type, attribute, ID) 

  Orthogroups <- read.table("Orthogroups.txt", sep = ":", col.names = c("Orthogroup", "genes"))

  subcluster_definition <- vector()

  if(arguments$subclusters > 1) {
    subcluster_definition <- read.table("Subclusters_MCL.txt", sep = "\t", col.names = c("Subcluster", "Elements"))
  }
  return(list(blast_results = blast_results, links = links, seqs = seqs, genes = genes, orthogroups_table = Orthogroups, subclusters = subcluster_definition))
}

perform_initial_clustering <- function(orthogroups_file = "Orthogroups.GeneCount.tsv") {
  OrthoFinder <- read.table(orthogroups_file, header = TRUE, check.names = FALSE, row.names = 1)
  OrthoFinder2 <- OrthoFinder[, !names(OrthoFinder) %in% c("Total")]
  transposed_counts <- t(OrthoFinder2)

  distance_mat <- dist(transposed_counts, method = 'binary')
  Hierar_cl <- hclust(distance_mat, method = 'average')
  Hierar_cl$height <- sort(Hierar_cl$height)
  
  fit_Orthofinder <- cutree(Hierar_cl, h = 0.99999999)
  Individual_clusters <- length(unique(fit_Orthofinder))
  
  return(list(transposed_counts = transposed_counts, orthofinder_counts = OrthoFinder2, distance_mat = distance_mat, Hierar_cl = Hierar_cl, cluster_fit = fit_Orthofinder, Individual_clusters = Individual_clusters))
}

core_genes_analysis <- function(ortho_counts, Cluster, Cluster_number = "", subcluster) {
  Whole_core <- ortho_counts[rowSums(ortho_counts < 1) <= ncol(ortho_counts) * 0.2, ]
  core_genes <- rownames(Whole_core)

  if (length(core_genes) > 0) {
    write.table(core_genes, file = paste("Core_genes-General", Cluster_number, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(Whole_core, file = paste("Core_genes-General", Cluster_number, ".csv", sep = ""), row.names = T, quote = F)
  }

  core_genes2 <- character()

  if(Cluster_number != "" || arguments$subclusters == 1) {
    fit <- cutree(Cluster, h = 0.8)

    if (length(unique(fit)) > 1) {
      for (SubClusterId in 1:length(unique(fit))) {
        OrthoFinder_subcluster <- subset(ortho_counts, select = names(fit[grep(SubClusterId, fit)]))
        if (ncol(OrthoFinder_subcluster) > 4) {
          OrthoFinder_subcluster <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster < 1) <= ncol(OrthoFinder_subcluster) * 0.2, ]
          if (nrow(OrthoFinder_subcluster) > 0) {
            core <- rownames(OrthoFinder_subcluster)
            core_genes2 <- c(core_genes2, core)
            write.table(core, file = paste("Core_genes-Specific", SubClusterId, Cluster_number, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
            write.csv(OrthoFinder_subcluster, file = paste("Core_genes-Specific", SubClusterId, Cluster_number, ".csv", sep = ""), row.names = T, quote = F)
          }
        }
      }
    }
  } else if (arguments$subclusters > 1) {
    for (SubClusterId in 1:length(subcluster$Subcluster)) {
      Elements_in <- unlist(strsplit(subcluster$Elements[SubClusterId], split = " "))
      OrthoFinder_subcluster <- subset(ortho_counts, select = Elements_in)
      if (ncol(OrthoFinder_subcluster) > 4) {
        OrthoFinder_subcluster <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster < 1) <= ncol(OrthoFinder_subcluster) * 0.2, ]
        if (nrow(OrthoFinder_subcluster) > 0) {
          core <- rownames(OrthoFinder_subcluster)
          core_genes2 <- c(core_genes2, core)
          write.table(core, file = paste("Core_genes-SubCluster", SubClusterId, Cluster_number, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
          write.csv(OrthoFinder_subcluster, file = paste("Core_genes-SubCluster", SubClusterId, Cluster_number, ".csv", sep = ""), row.names = T, quote = F)
        }
      }
    }
  }

  return(list(general_core = core_genes, specific_core = core_genes2))
}

analyze_individual_clusters <- function(Individual_clusters, transposase_OrthoFinder, fit_Orthofinder, seq_data, gene_data, blast_links) {
  Selected_clusters <- numeric()
  for (ClusterId in 1:Individual_clusters) {
    Cluster_matrix <- transposase_OrthoFinder[grep(ClusterId, fit_Orthofinder), ]
    if (is.array(Cluster_matrix)) {
      Cluster_matrix <- Cluster_matrix[, colSums(abs(Cluster_matrix)) > 0]
    } else {
      Cluster_matrix <- matrix(, nrow = 1, ncol = 1)
      rownames(Cluster_matrix) <- names(fit_Orthofinder[grep(ClusterId, fit_Orthofinder)])
    }
    
    if (ncol(Cluster_matrix) > 2) {
      Cluster_distance <- dist(Cluster_matrix, method = 'binary')
      Cluster_Hier <- hclust(Cluster_distance, method = 'average')
      my_tree <- as.phylo(Cluster_Hier)
      
      if (Individual_clusters >= 2) {
        write.tree(phy = my_tree, file = paste("CargoHierarchicalTree_Cluster", ClusterId, ".nwk", sep = ""))
        write.table(Cluster_Hier$labels, file = paste("Cluster", ClusterId, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
      } else {
        write.tree(phy = my_tree, file = "CargoHierarchicalTree.nwk")
      }
      
      Cluster_melt <- reshape2::melt(Cluster_matrix, as.is = T)
      Profile_plot <- ggplot(data = Cluster_melt, aes(x = reorder(Var2, value), y = Var1, fill = value)) +
        geom_tile() +
        scale_fill_viridis(option = "rocket", direction = -1) +
        labs(y = "element", x = "Orthogroup") + theme(axis.text.x = element_blank())

      plot_width <- min(49, ncol(Cluster_matrix) * 0.2 + 2)
      plot_height <- min(49, nrow(Cluster_matrix))
      if (Individual_clusters >= 2) {
        ggsave(Profile_plot, filename = paste("CargoHeatmap_Cluster", ClusterId, ".svg", sep = ""), width = plot_width, height = plot_height)
      } else {
        ggsave(Profile_plot, filename = "CargoHeatmap.svg", width = plot_width, height = plot_height)
      }
    } else {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "WARNING: There are issues with the following element: ", as.character(rownames(Cluster_matrix)), "\n", sep = ""))
    }
  }
  return(Selected_clusters)
}

plot_cluster_synteny <- function(Cluster_matrix, orthogroup_table, core_genes, seq_data, gene_data, blast_links, Cluster_number = "", SubCluster_number = "") {
  Cluster_distance <- dist(Cluster_matrix, method = 'binary')
  Cluster_Hier <- hclust(Cluster_distance, method = 'average')
  my_tree <- as.phylo(Cluster_Hier)
  tree_sorted <- ladderize(my_tree, right = FALSE)
  selected_seqs <- row.names(Cluster_matrix)
  seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs)
  genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs)

  if (length(core_genes$specific_core) > 0) {
    specific_table <- orthogroup_table[orthogroup_table$Orthogroup %in% core_genes$specific_core, ]
    specific_vector <- unlist(strsplit(specific_table[, 2], split = " "))
    specific_vector <- specific_vector[specific_vector != ""]
    genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% specific_vector, "specific", attribute))
  }

  if (length(core_genes$general_core) > 0) {
    general_table <- orthogroup_table[orthogroup_table$Orthogroup %in% core_genes$general_core, ]
    general_vector <- unlist(strsplit(general_table[, 2], split = " "))
    general_vector <- general_vector[general_vector != ""]
    genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% general_vector, "general", attribute))
  }

  links_filtered <- blast_links %>% filter(seq_id %in% selected_seqs | seq_id2 %in% selected_seqs)
  p_tree_base <- ggtree::ggtree(tree_sorted, layout = "rectangular")
  tree_y_coords <- p_tree_base$data %>% filter(isTip) %>% select(seq_id = label, y) %>% mutate(y = max(y) - y + 1)
  ordered_seqs <- seqs_filtered %>% inner_join(tree_y_coords, by = "seq_id") %>% arrange(y)
  max_seq_len <- max(ordered_seqs$length, na.rm = TRUE)

  p_tree <- p_tree_base + ggtree::geom_tiplab(align = TRUE, size = ifelse(nrow(Cluster_matrix) < 100, 2.5, 3), family = "mono") + ggtree::theme_tree2()

  p_genome <- suppressMessages(gggenomes(seqs = ordered_seqs, links = links_filtered, genes = genes_filtered) +
    geom_seq(aes(y = y)) +
    geom_gene(aes(y = y, fill = attribute), show.legend = T) +
    scale_fill_manual(name = "Core genes", values = c("general" = "red4", "specific" = "green4"), na.value = "cornsilk3", limits = c("general", "specific")) +
    new_scale_fill() +
    geom_link(aes(y = y, fill = pident), colour = NA) +
    scale_fill_continuous(name = "Alignment Identity (%)") +
    geom_bin_label(aes(y = y), x = -10, size = ifelse(nrow(Cluster_matrix) < 50, 2.5, 3)) +
    scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max_seq_len)) +
    scale_y_continuous(expand = expansion(mult = 1.0001 * (nrow(Cluster_matrix) ^ -1.1728))))
      
  final_plot <- p_tree + p_genome
  plot_width_final <- min(49, max(16, round(max_seq_len * 0.0001) + round(max(tree_sorted$edge.length, na.rm = TRUE) * 10)))
  plot_height_final <- min(49, nrow(Cluster_matrix))
      
  if (Cluster_number != "") {
    ggsave(final_plot, filename = paste("CargoSynteny_Cluster", Cluster_number, ".svg", sep = ""), width = plot_width_final, height = plot_height_final, limitsize = FALSE)
  } else if (SubCluster_number != "") {
    ggsave(final_plot, filename = paste("CargoSynteny_SubCluster", SubCluster_number, ".svg", sep = ""), width = plot_width_final, height = plot_height_final, limitsize = FALSE)
  } else {
    ggsave(final_plot, filename = "CargoSynteny.svg", width = plot_width_final, height = plot_height_final, limitsize = FALSE)
  }
}

check_nesting <- function(ortho_counts, seq_data, gene_data, blast_links) {
  for (Element in 1:ncol(ortho_counts)) { 
    Element_name <- colnames(ortho_counts[Element])
    Genes_element <- rownames(ortho_counts[Element_name] %>% filter(!if_all(everything(), ~ .x == 0)))
    reduced_orthofinder <- ortho_counts[Genes_element, ]
    reduced_orthofinder2 <- reduced_orthofinder[, colSums(reduced_orthofinder < 1) < nrow(reduced_orthofinder) * 0.2]
    
    if (is.data.frame(reduced_orthofinder2)) {
      reduced_orthofinder2 <- reduced_orthofinder2 %>% select(-all_of(Element_name))
      
      if (is.data.frame(reduced_orthofinder2)) {
        for (comparison in 1:ncol(reduced_orthofinder2)) {
          Element_compare <- colnames(reduced_orthofinder2[comparison])
          Genes_compare <- rownames(ortho_counts[Element_compare] %>% filter(!if_all(everything(), ~ .x == 0)))
          
          if (length(Genes_compare) > (length(Genes_element) * 1.5)) {
            links_filtered2 <- blast_links %>% filter(qseqid %in% Element_name & sseqid %in% Element_compare) %>%
              select(seq_id = qseqid, start = qstart, end = qend, seq_id2 = sseqid, start2 = sstart, end2 = send, pident)
            
            if (nrow(links_filtered2) == 0) {
              next
            }
            
            seq_idLength <- as.numeric(seq_data[seq_data$seq_id %in% unique(links_filtered2$seq_id), 2])
            seq_id2Length <- as.numeric(seq_data[seq_data$seq_id %in% unique(links_filtered2$seq_id2), 2])
            is_nested_q <- any(links_filtered2$start < seq_idLength * 0.2) && any(links_filtered2$end > seq_idLength * 0.8)
            is_not_nested_s <- any(links_filtered2$start2 < seq_id2Length * 0.2) 

            if (is_nested_q && ! is_not_nested_s) {
              cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Element ", Element_name, " is nested in element ", Element_compare, "\n", sep = ""))
              selected_seq <- c(Element_name, Element_compare)
              seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seq)
              genes_filtered <- gene_data %>% filter(seq_id %in% selected_seq)
              links_filtered <- blast_links %>% filter(qseqid %in% selected_seq & sseqid %in% selected_seq) %>%
                select(seq_id = qseqid, start = qstart, end = qend, seq_id2 = sseqid, start2 = sstart, end2 = send, pident)
              ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, selected_seq))
              p_genome <- gggenomes(seqs = ordered_seqs, links = links_filtered, genes = genes_filtered) +
                geom_seq(aes(y = y)) + geom_gene(aes(y = y)) + geom_link(aes(y = y, fill = pident), colour = NA) +
                geom_bin_label(aes(y = y), x = -10) + scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
                theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
              ggsave(p_genome, filename = paste("IndividualNestingEvent-", Element_name, "in", Element_compare, ".svg", sep = ""), 
                     width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)), height = min(49, length(selected_seq)), limitsize = FALSE)
            }
          }
        }
      }
    }
  }
}

gene_movement_analysis <- function(subcluster_number, matrix, ortho_counts, core_genes, big = FALSE, Cluster_number = "") {
  test_NbClust <- vector()
  if (big) {
    elements_cluster <- rownames(matrix)
    elements_other <- colnames(matrix)
  } else {
    elements_cluster <- rownames(matrix)
    elements_other <- tail(colnames(matrix), n = 1)
  }
  selected_seqs2 <- c(elements_cluster, elements_other)
  OrthoFinder_subcluster_cluster <- subset(ortho_counts, select = elements_cluster) %>% filter(!if_all(everything(), ~ .x == 0))
  OrthoFinder_subcluster_other <- subset(ortho_counts, select = elements_other) %>% filter(!if_all(everything(), ~ .x == 0))
  Gene_movement <- intersect(rownames(OrthoFinder_subcluster_cluster), rownames(OrthoFinder_subcluster_other))
  Nesting_test <- Gene_movement
  Gene_movement <- Gene_movement[!Gene_movement %in% core_genes$general_core]
  
  if (length(Nesting_test) >= nrow(OrthoFinder_subcluster_cluster) * 0.8 || length(Nesting_test) >= nrow(OrthoFinder_subcluster_other) * 0.8) {
    single_movements <- TRUE
    multiple_movements <- FALSE
    write.table(selected_seqs2, file = paste(Cluster_number, "SubCluster", subcluster_number, "_nestingEvent.txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
  } else if (length(Gene_movement) > 0 && big) {
    single_movements <- TRUE
    OrthoFinder_subcluster <- ortho_counts[Gene_movement, selected_seqs2]
    OrthoFinder_subcluster <- OrthoFinder_subcluster[, colSums(OrthoFinder_subcluster) > 0]
    elements_cluster <- elements_cluster[elements_cluster %in% colnames(OrthoFinder_subcluster)]
    elements_other <- elements_other[elements_other %in% colnames(OrthoFinder_subcluster)]
    trans_OF <- t(OrthoFinder_subcluster)
    
    if (nrow(trans_OF) > 3 && ncol(trans_OF) > 3) {
      distance_mat <- dist(trans_OF, method = 'binary')
      test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(trans_OF) - 1), index = "ball")
      multiple_movements <- TRUE
    } else {
      elements_cluster <- elements_cluster[elements_cluster %in% colnames(OrthoFinder_subcluster)]
      elements_other <- elements_other[elements_other %in% colnames(OrthoFinder_subcluster)]
      multiple_movements <- FALSE
    }
    write.table(Gene_movement, file = paste(Cluster_number, "SubCluster", subcluster_number, "_moveOrthologs.txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(OrthoFinder_subcluster, file = paste(Cluster_number, "SubCluster", subcluster_number, "_moveOrthologsTable.csv", sep = ""), row.names = T, quote = F)
  } else if (length(Gene_movement) > 0) {
    OrthoFinder_subcluster <- ortho_counts[Gene_movement, selected_seqs2]
    elements_cluster <- elements_cluster[elements_cluster %in% colnames(OrthoFinder_subcluster)]
    elements_other <- elements_other[elements_other %in% colnames(OrthoFinder_subcluster)]
    single_movements <- TRUE
    multiple_movements <- FALSE
    write.table(Gene_movement, file = paste(Cluster_number, "SubCluster", subcluster_number, "_moveOrthologs.txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(OrthoFinder_subcluster, file = paste(Cluster_number, "SubCluster", subcluster_number, "_moveOrthologsTable.csv", sep = ""), row.names = T, quote = F)
  } else {
    single_movements <- FALSE
    multiple_movements <- FALSE
  }
  return(list(movement = single_movements, multiple = multiple_movements, orthogroups = Gene_movement, elements_in = elements_cluster, elements_out = elements_other, K_analysis = test_NbClust))
}

plot_subcluster_synteny <- function(subcluster_number, orthogroup_table, dist_matrix, cluster_fit, ortho_counts, core_genes, seq_data, gene_data, blast_links, cluster_matrix, Cluster_number = "", subcluster) {
  for (ClusterId in 1:subcluster_number) {
    if(Cluster_number != "") {
      matrix <- dist_matrix[grep(ClusterId, cluster_fit), ]
      Cluster_elements <- rownames(matrix)
      if (is.null(Cluster_elements)) {
        Cluster_elements <- names(grep(ClusterId, cluster_fit, value = TRUE))
      }
    } else {
      Cluster_elements <- unlist(strsplit(subcluster$Elements[ClusterId], split = " "))
      matrix <- dist_matrix[Cluster_elements,]
    }
    
    write.table(Cluster_elements, file = paste(Cluster_number, "SubCluster", ClusterId, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
    selected_seqs2_defined <- FALSE
    multiple_movements <- FALSE
    
    if (is.matrix(matrix)) {
      subcluster_matrix <- cluster_matrix[rownames(cluster_matrix) %in% Cluster_elements,]

      links <- blast_links %>%
          select(seq_id = qseqid, start = qstart, end = qend, seq_id2 = sseqid, start2 = sstart, end2 = send, pident)

      plot_cluster_synteny(
        Cluster_matrix = subcluster_matrix, 
        orthogroup_table = orthogroup_table, 
        core_genes = core_genes, 
        seq_data = seq_data, 
        gene_data = gene_data, 
        blast_links = links,
        SubCluster_number = ClusterId
      )

      reduce_matrix <- matrix[, colSums(matrix) < nrow(matrix) * 0.98]
      
      if (ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
        reduce_matrix2 <- reduce_matrix[, (nrow(reduce_matrix) + 1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2) * 0.98, ]
        
        if (is.matrix(reduce_matrix3)) {
          reduce_matrix3 <- reduce_matrix3[, colSums(reduce_matrix3 < 1) > 1]
          reduce_matrix3 <- reduce_matrix3[, names(sort(colSums(reduce_matrix3), decreasing = F))]
          Movement <- gene_movement_analysis(subcluster_number = ClusterId, matrix = reduce_matrix3, ortho_counts = ortho_counts, core_genes = core_genes, big = TRUE, Cluster_number = Cluster_number)
          selected_seqs2 <- c(Movement$elements_in, Movement$elements_out)
          selected_seqs2_defined <- Movement$movement
          multiple_movements <- Movement$multiple
          elements_cluster <- Movement$elements_in
          elements_other <- Movement$elements_out
          test_NbClust <- Movement$K_analysis
        }
      } else if (ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
        reduce_matrix <- reduce_matrix[order(reduce_matrix[, ncol(reduce_matrix)], decreasing = TRUE), ]
        Movement <- gene_movement_analysis(subcluster_number = ClusterId, matrix = reduce_matrix, ortho_counts = ortho_counts, core_genes = core_genes, Cluster_number = Cluster_number)
        selected_seqs2 <- c(Movement$elements_in, Movement$elements_out)
        selected_seqs2_defined <- Movement$movement
      } else {
        reduce_matrix <- matrix[, colSums(matrix) < nrow(matrix) * 0.99]
        
        if (ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
          reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
          reduce_matrix2 <- reduce_matrix[, (nrow(reduce_matrix) + 1):ncol(reduce_matrix)]
          reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2) * 0.99, ]
          
          if (is.matrix(reduce_matrix3)) {
            reduce_matrix3 <- reduce_matrix3[, colSums(reduce_matrix3 < 1) > 1]
            reduce_matrix3 <- reduce_matrix3[, names(sort(colSums(reduce_matrix3), decreasing = F))]
            Movement <- gene_movement_analysis(subcluster_number = ClusterId, matrix = reduce_matrix3, ortho_counts = ortho_counts, core_genes = core_genes, big = TRUE, Cluster_number = Cluster_number)
            selected_seqs2 <- c(Movement$elements_in, Movement$elements_out)
            selected_seqs2_defined <- Movement$movement
            multiple_movements <- Movement$multiple
            elements_cluster <- Movement$elements_in
            elements_other <- Movement$elements_out
            test_NbClust <- Movement$K_analysis
          }  
        } else if (ncol(reduce_matrix) > nrow(reduce_matrix)) {
          reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
          reduce_matrix <- reduce_matrix[order(reduce_matrix[, ncol(reduce_matrix)], decreasing = TRUE), ]
          Movement <- gene_movement_analysis(subcluster_number = ClusterId, matrix = reduce_matrix, ortho_counts = ortho_counts, core_genes = core_genes, Cluster_number = Cluster_number)
          selected_seqs2 <- c(Movement$elements_in, Movement$elements_out)
          selected_seqs2_defined <- Movement$movement
        }
      }
      
      if (selected_seqs2_defined & !multiple_movements) {
        seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs2)
        genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs2)
        
        if (length(core_genes$specific_core) > 0) {
          spec_table <- orthogroup_table[orthogroup_table$Orthogroup %in% core_genes$specific_core, ]
          spec_vector <- unlist(strsplit(spec_table[, 2], split = " "))
          spec_vector <- spec_vector[spec_vector != ""]
          genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% spec_vector, "specific", attribute))
        }
        
        if (length(core_genes$general_core) > 0) {
          gen_table <- orthogroup_table[orthogroup_table$Orthogroup %in% core_genes$general_core, ]
          gen_vector <- unlist(strsplit(gen_table[, 2], split = " "))
          gen_vector <- gen_vector[gen_vector != ""]
          genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% gen_vector, "general", attribute))
        }
        
        mov_table <- orthogroup_table[orthogroup_table$Orthogroup %in% Movement$orthogroups, ]
        mov_vector <- unlist(strsplit(mov_table[, 2], split = " "))
        mov_vector <- mov_vector[mov_vector != ""]
        genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% mov_vector, "Movement associated", attribute))
        
        links_filtered <- blast_links %>% filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(seq_id = qseqid, start = qstart, end = qend, seq_id2 = sseqid, start2 = sstart, end2 = send, pident)
        
        ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, selected_seqs2))
        p_genome <- gggenomes(seqs = ordered_seqs, links = links_filtered, genes = genes_filtered) +
          geom_seq(aes(y = y)) + geom_gene(aes(y = y, fill = attribute), show.legend = T) +
          scale_fill_manual(name = "Core genes", values = c("general" = "red4", "specific" = "green4", "Movement associated" = "darkorchid"), na.value = "cornsilk3", limits = c("general", "specific", "Movement associated")) +
          new_scale_fill() + geom_link(aes(y = y, fill = pident), colour = NA) + scale_fill_continuous(name = "Alignment Identity (%)") +
          geom_bin_label(aes(y = y), x = -10) + scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
        
        ggsave(p_genome, filename = paste(Cluster_number, "MovementSynteny_SubCluster", ClusterId, ".svg", sep = ""), 
               width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)), height = min(49, length(selected_seqs2)), limitsize = FALSE)
      } else if (selected_seqs2_defined & multiple_movements) {
        k_out_cluster <- unique(test_NbClust$Best.partition[elements_other])
        k_cluster <- unique(test_NbClust$Best.partition[elements_cluster])
        
        for (i in k_cluster) {
          for (j in k_out_cluster) {
            cl_seqs <- names(test_NbClust$Best.partition[grep(i, test_NbClust$Best.partition)][elements_cluster])
            cl_seqs <- cl_seqs[!is.na(cl_seqs)]
            
            oth_seqs <- names(test_NbClust$Best.partition[grep(j, test_NbClust$Best.partition)][elements_other])
            oth_seqs <- oth_seqs[!is.na(oth_seqs)]
            
            OF_sub_cl <- subset(ortho_counts, select = cl_seqs) %>% filter(!if_all(everything(), ~ .x == 0))
            OF_sub_oth <- subset(ortho_counts, select = oth_seqs) %>% filter(!if_all(everything(), ~ .x == 0))
            
            G_mov <- intersect(rownames(OF_sub_cl), rownames(OF_sub_oth))
            OF_sub_cl <- OF_sub_cl[G_mov,]
            OF_sub_cl <- OF_sub_cl[,  colSums(abs(OF_sub_cl)) > 0]
            OF_sub_oth <- OF_sub_oth[G_mov,]
            OF_sub_oth <- OF_sub_oth[,  colSums(abs(OF_sub_oth)) > 0]

            sel_seq2 <- c(colnames(OF_sub_cl), colnames(OF_sub_oth))
            #sel_seq2 <- c(cl_seqs, oth_seqs)
            sel_seq2 <- sel_seq2[!is.na(sel_seq2)]
            
            if (length(sel_seq2) >= 2) {
              seqs_filtered <- seq_data %>% filter(seq_id %in% sel_seq2)
              genes_filtered <- gene_data %>% filter(seq_id %in% sel_seq2)
              links_filtered <- blast_links %>% filter(qseqid %in% sel_seq2 & sseqid %in% sel_seq2) %>%
                select(seq_id = qseqid, start = qstart, end = qend, seq_id2 = sseqid, start2 = sstart, end2 = send, pident)
              
              if (length(core_genes$specific_core) > 0) {
                spec_table <- orthogroup_table[orthogroup_table$Orthogroup %in% core_genes$specific_core, ] 
                spec_vector <- unlist(strsplit(spec_table[, 2], split = " "))
                spec_vector <- spec_vector[spec_vector != ""]
                genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% spec_vector, "specific", attribute))
              }
        
              if (length(core_genes$general_core) > 0) {
                gen_table <- orthogroup_table[orthogroup_table$Orthogroup %in% core_genes$general_core, ]
                gen_vector <- unlist(strsplit(gen_table[, 2], split = " "))
                gen_vector <- gen_vector[gen_vector != ""]
                genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% gen_vector, "general", attribute))
              }
        
              mov_table <- orthogroup_table[orthogroup_table$Orthogroup %in% G_mov, ]
              mov_vector <- unlist(strsplit(mov_table[, 2], split = " "))
              mov_vector <- mov_vector[mov_vector != ""]
              genes_filtered <- genes_filtered %>% mutate(attribute = ifelse(ID %in% mov_vector, "Movement associated", attribute))
              
              if (nrow(links_filtered) >= 1) {
                ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, sel_seq2))
                p_genome <- gggenomes(seqs = ordered_seqs, links = links_filtered, genes = genes_filtered) +
                  geom_seq(aes(y = y)) + geom_gene(aes(y = y, fill = attribute), show.legend = T) +
                  scale_fill_manual(name = "Core genes", values = c("general" = "red4", "specific" = "green4", "Movement associated" = "darkorchid"), na.value = "cornsilk3", limits = c("general", "specific", "Movement associated")) +
                  new_scale_fill() + geom_link(aes(y = y, fill = pident), colour = NA) + scale_fill_continuous(name = "Alignment Identity (%)") +
                  geom_bin_label(aes(y = y), x = -10) + scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
                  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
                
                ggsave(p_genome, filename = paste(Cluster_number, "MovementSynteny_SubCluster", ClusterId, "-", i, "vs", j, ".svg", sep = ""), 
                       width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)), height = min(49, length(sel_seq2)), limitsize = FALSE)

                New_ortho_count <- ortho_counts[G_mov,sel_seq2]
                write.csv(New_ortho_count, file = paste(Cluster_number, "SubCluster", ClusterId, "-", i, "vs", j, "_moveOrthologsMatrix.csv", sep = ""), row.names = T, quote = F)
              }
            }
          }
        }
      } else {
        cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Subcluster ", ClusterId, " do not have identifiable putative movement event.", "\n", sep = ""))
      }
    } else {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Subcluster ", ClusterId, " can not be analyzed automatically.", "\n", sep = ""))
    }
  }
}

get_discordant <- function(tree1, tree2) {
  ct <- table(tree1, tree2)
  bm <- apply(ct, 1, which.max)
  t_aln <- tree1
  
  for (i in 1:length(bm)) {
    original_label <- names(bm)[i]
    new_label <- as.numeric(bm[i])
    t_aln[tree1 == original_label] <- new_label
  }
  
  mismatched_indices <- which(tree2 != t_aln)
  discordant_names <- names(tree2)[mismatched_indices]
  return(discordant_names)
}

# ==============================================================================
# MAIN EXECUTION
# ==============================================================================

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Reading data...", "\n", sep = ""))
data_list <- load_and_preprocess_data(min_pident = arguments$pident)

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Identifying if there are separate clusters", "\n", sep = ""))
clustering_data <- perform_initial_clustering()

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Checking for nesting events ", "\n", sep = ""))
check_nesting(ortho_counts = clustering_data$orthofinder_counts, seq_data = data_list$seqs, gene_data = data_list$genes, blast_links = data_list$blast_results)

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Analyzing ", clustering_data$Individual_clusters, " individual cluster(s) identified", "\n", sep = ""))
Selected_Cluster_ID <- analyze_individual_clusters(
  Individual_clusters = clustering_data$Individual_clusters, 
  transposase_OrthoFinder = clustering_data$transposed_counts, 
  fit_Orthofinder = clustering_data$cluster_fit, 
  seq_data = data_list$seqs, 
  gene_data = data_list$genes, 
  blast_links = data_list$blast_results
)

if (clustering_data$Individual_clusters == 1) {
  if (arguments$subclusters >= 2) {
    cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Analyzing if there are possible core genes", "\n", sep = ""))
    core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl, subcluster= data_list$subclusters)
    
    cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Generating main plot", "\n", sep = ""))
    plot_cluster_synteny(
      Cluster_matrix = clustering_data$transposed_counts, 
      orthogroup_table = data_list$orthogroups_table, 
      core_genes = core_genes, 
      seq_data = data_list$seqs, 
      gene_data = data_list$genes, 
      blast_links = data_list$links
    )

    if (arguments$captainRemoval <= 1) {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Analyzing ", arguments$subclusters, " subclusters", "\n", sep = ""))
      fit <- cutree(clustering_data$Hierar_cl, k = arguments$subclusters)
      plot_subcluster_synteny(
        subcluster_number = arguments$subclusters, 
        orthogroup_table = data_list$orthogroups_table, 
        dist_matrix = as.matrix(clustering_data$distance_mat), 
        cluster_fit = fit, 
        ortho_counts = clustering_data$orthofinder_counts, 
        core_genes = core_genes, 
        seq_data = data_list$seqs, 
        gene_data = data_list$genes, 
        blast_links = data_list$blast_results,
        subcluster = data_list$subclusters,
        cluster_matrix = clustering_data$transposed_counts
      )
    } else {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Analyzing subclusters with NbClust optimization", "\n", sep = ""))
      test_NbClust <- NbClust(clustering_data$distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(clustering_data$distance_mat) - 1), index = "ball")
      K_number <- round(mean(c(length(unique(test_NbClust$Best.partition)), arguments$subclusters)))
      fit <- cutree(clustering_data$Hierar_cl, k = K_number)
      plot_subcluster_synteny(
        subcluster_number = K_number, 
        orthogroup_table = data_list$orthogroups_table, 
        dist_matrix = as.matrix(clustering_data$distance_mat), 
        cluster_fit = fit, 
        ortho_counts = clustering_data$orthofinder_counts, 
        core_genes = core_genes, 
        seq_data = data_list$seqs, 
        gene_data = data_list$genes, 
        blast_links = data_list$blast_results,
        subcluster = data_list$subclusters,
        cluster_matrix = clustering_data$transposed_counts
      )
    }
  } else {
    cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Generating single cluster plot", "\n", sep = ""))
    core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl)
    plot_cluster_synteny(
      Cluster_matrix = clustering_data$transposed_counts, 
      orthogroup_table = data_list$orthogroups_table, 
      core_genes = core_genes, 
      seq_data = data_list$seqs, 
      gene_data = data_list$genes, 
      blast_links = data_list$links
    )
  }
} else if (clustering_data$Individual_clusters > 1) {
  if (arguments$subclusters >= 2) {
    for (MainClusterID in Selected_Cluster_ID) {
      Cluster_matrix <- clustering_data$transposed_counts[grep(MainClusterID, clustering_data$cluster_fit), ]
      
      if (!is.matrix(Cluster_matrix) || nrow(Cluster_matrix) < 4) {
        cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Cluster ", MainClusterID, " can not be analyzed automatically", "\n", sep = ""))
        next
      }
      
      Cluster_matrix <- Cluster_matrix[, colSums(abs(Cluster_matrix)) > 0]
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Analyzing core genes in Cluster ", MainClusterID, "\n", sep = ""))
      
      New_ortho_counts <- clustering_data$orthofinder_counts[, grep(MainClusterID, clustering_data$cluster_fit)]
      distance_mat <- dist(Cluster_matrix, method = 'binary')
      Hierar_cl <- hclust(distance_mat, method = 'average')
      
      core_genes <- core_genes_analysis(ortho_counts = New_ortho_counts, Cluster = Hierar_cl, Cluster_number = paste("_Cluster", MainClusterID, sep = ""))
      plot_cluster_synteny(Cluster_matrix = Cluster_matrix, orthogroup_table = data_list$orthogroups_table, core_genes = core_genes, seq_data = data_list$seqs, gene_data = data_list$genes, blast_links = data_list$links, Cluster_number = MainClusterID)
      
      test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(Cluster_matrix) - 1), index = "ball")
      K_number <- round(mean(c(length(unique(test_NbClust$Best.partition)), length(unique(cutree(Hierar_cl, h = 0.6))))))
      fit <- cutree(Hierar_cl, k = K_number)
      
      plot_subcluster_synteny(
        subcluster_number = K_number, 
        orthogroup_table = data_list$orthogroups_table, 
        dist_matrix = as.matrix(distance_mat), 
        cluster_fit = fit, 
        ortho_counts = clustering_data$orthofinder_counts, 
        core_genes = core_genes, 
        seq_data = data_list$seqs, 
        gene_data = data_list$genes, 
        blast_links = data_list$blast_results, 
        Cluster_number = paste("Cluster", MainClusterID, "-", sep = "")
      )
    }
  }
} else {
  cat(paste("    [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "No subcluster to analyzed", "\n", sep = ""))
}

if (clustering_data$Individual_clusters == 1) {
  cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Analyzing Captain vs Cargo discordance", "\n", sep = ""))
  
  if (file.exists("CaptainPhylogeny.nw")) {
    phylo_tree <- read.tree("CaptainPhylogeny.nw")
    dist_matrix <- cophenetic.phylo(phylo_tree)
    normalized_dist_matrix <- (dist_matrix - min(dist_matrix)) / (max(dist_matrix) - min(dist_matrix))
    hclust_tree <- hclust(as.dist(normalized_dist_matrix), method = "average")
    hclust_tree$height <- sort(hclust_tree$height)
    
    if (length(phylo_tree$tip.label) > 4) {
      test_nb <- NbClust(scale(normalized_dist_matrix), method = "average", min.nc = 1, max.nc = min(6, nrow(normalized_dist_matrix) - 1), index = "ball")
      K_num <- round(mean(c(length(unique(test_nb$Best.partition)), length(unique(cutree(hclust_tree, h = 0.5))))))
    } else { 
      K_num <- 2 
    }
    
    fit_tree <- cutree(hclust_tree, k = K_num)
    fit_tree <- fit_tree[order(names(fit_tree))]
    
    fit_cargo <- cutree(clustering_data$Hierar_cl, k = K_num)
    fit_cargo <- fit_cargo[order(names(fit_cargo))]

    discordant_elements <- unique(c(get_discordant(fit_tree, fit_cargo), get_discordant(fit_cargo, fit_tree)))

    if (length(discordant_elements) >= 1) {
      cat(paste("    [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] There ", length(discordant_elements), " discordant element(s) identified\n", sep = ""))
      write.table(discordant_elements, file = "Discordant_elements.txt", sep = '\t', row.names = F, col.names = F, quote = F)
    } else {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] No discordance detected\n", sep = ""))
    }

    svg(filename = "CaptainVsCargo_tanglegram.svg", width = 16, height = min(49, length(fit_cargo) / 2))
    dendlist(clustering_data$Hierar_cl, hclust_tree) %>% 
      ladderize %>% 
      untangle(method = "step1side", k_seq = K_num:(length(fit_cargo) - 1)) %>% 
      set("branches_k_color", k = K_num) %>% 
      tanglegram(faster = TRUE, main_left = "Cargo Cluster", main_right = "Captain tree")
    invisible(dev.off())
  }
}
