# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory [default %default]",metavar="PATH"),
  make_option(c("-s", "--subclusters"), type="integer", action = "store", default=1,
              help=" Number of subclusters get in the synteny analysis. Minimum value is equal to 1 when no subcluster were identified [default %default]", metavar="number"),
  make_option(c("-c", "--captainRemoval"), type="integer", action = "store", default=0,
              help=" Number of elements removes from the original cluster [default %default]", metavar="number")
  make_option(c("-p", "--pident"), type="numeric", action = "store", default=70,
              help=" Minimum percentage identity of blast results to be included as links for nucleotide synteny visualization [default %default]", metavar="number")
)

# Parse the command-line arguments
arguments <- parse_args(OptionParser(option_list = option_list))

# Check for required arguments and flags
if (is.null(arguments$directory)) {
  stop("Error: directory must be provided.", call.=FALSE)
} else {
  setwd(arguments$directory)
}

# Check software installation
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(viridis))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(gggenomes))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(NbClust))
suppressPackageStartupMessages(library(syntenet))
suppressPackageStartupMessages(library(ggnewscale))

if (!requireNamespace("ggplot2", quietly = TRUE)) {
   stop("Package \"ggplot2\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("ape", quietly = TRUE)) {
   stop("Package \"ape\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("ggtree", quietly = TRUE)) {
   stop("Package \"ggtree\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("reshape2", quietly = TRUE)) {
   stop("Package \"reshape2\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("viridis", quietly = TRUE)) {
   stop("Package \"viridis\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
   stop("Package \"dplyr\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("gggenomes", quietly = TRUE)) {
   stop("Package \"gggenomes\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("scales", quietly = TRUE)) {
   stop("Package \"scales\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("dendextend", quietly = TRUE)) {
   stop("Package \"dendextend\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("NbClust", quietly = TRUE)) {
   stop("Package \"NbClust\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("syntenet", quietly = TRUE)) {
   stop("Package \"syntenet\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("ggnewscale", quietly = TRUE)) {
   stop("Package \"ggnewscale\" not installed. Please install it to run this script.", call. = FALSE)
}

# Function block

load_and_preprocess_data <- function(
    blast_file = "Blast_CleanResults.txt",
    gff_file = "Final_model.gff",
    min_pident
) {
  
  #Define column names and specifications based on standard BLAST outfmt 6 + qlen and slen
  blast_col_names <- c("qseqid", "sseqid", "qstart", "qend", "sstart", "send",
                       "pident", "length", "qlen", "slen")
  blast_col_classes <- c("character", "character", rep("numeric", 8))
  
  # Read and Preprocess BLAST Results (blast_results)
  blast_results <- read.delim(
    blast_file,
    header = FALSE,
    sep = "\t",
    col.names = blast_col_names,
    colClasses = blast_col_classes
  ) %>%
    # Filter out self-hits and low-identity hits
    filter(qseqid != sseqid, pident > min_pident)
  
  # Create 'links' for gggenomes
  links <- blast_results %>%
    mutate(alignment_id = row_number()) %>%
    select(
      seq_id = qseqid, start = qstart, end = qend,
      seq_id2 = sseqid, start2 = sstart, end2 = send,
      pident, length
    )
  
  # Create 'seqs' for gggenomes
  seqs <- bind_rows(
    blast_results %>% select(seq_id = qseqid, length = qlen),
    blast_results %>% select(seq_id = sseqid, length = slen)
  ) %>%
    distinct() %>%
    group_by(seq_id) %>%
    summarize(length = max(length)) %>%
    ungroup()
  
  # Extract gene locations and format for gggenomes visualization
  annotation <- gff2GRangesList("Gff/")
  genes <- as.data.frame(unlist(annotation)) %>%
    filter(type == "gene") %>%
    mutate(
      seqnames = as.character(seqnames),
      type = "CDS",
      attribute = NA) %>%
    select(seq_id = seqnames, start, end, strand, type, attribute,ID) 

  Orthogroups <- read.table("Orthogroups.txt", sep = ":", col.names = c("Orthogroup","genes"))
  
  return(list(
    blast_results = blast_results,
    links = links,
    seqs = seqs,
    genes = genes,
    orthogroups_table = Orthogroups
  ))
}

perform_initial_clustering <- function(
    orthogroups_file = "Orthogroups.GeneCount.tsv"
) {
  
  # Read OrthoFinder Gene Counts
  OrthoFinder <- read.table(
    orthogroups_file,
    header = TRUE,
    check.names = FALSE,
    row.names = 1
  )

  # Remove the 'Total' column
  OrthoFinder2 <- OrthoFinder[, !names(OrthoFinder) %in% c("Total")]
  
  # Transpose: We want genomes as ROWS and orthogroups as COLUMNS
  transposed_counts <- t(OrthoFinder2)

  # Calculate distance using 'binary' (Jaccard-like): looks at presence vs absence
  distance_mat <- dist(transposed_counts, method = 'binary')
  
  # Perform UPGMA clustering
  Hierar_cl <- hclust(distance_mat, method = 'average')
  Hierar_cl$height <- sort(Hierar_cl$height)
  
  # Cut tree at nearly 1.0 height to identify if there has been any break in the cluster connection
  fit_Orthofinder <- cutree(Hierar_cl, h = 0.99999999)
  Individual_clusters <- length(unique(fit_Orthofinder))
  
  return(list(
    transposed_counts = transposed_counts,
    orthofinder_counts = OrthoFinder2,
    distance_mat = distance_mat,
    Hierar_cl = Hierar_cl,
    cluster_fit = fit_Orthofinder,
    Individual_clusters = Individual_clusters
  ))
}

core_genes_analysis <- function(    
    ortho_counts,
    Cluster,
    Cluster_number = ""
) {
  
  

  # Define 'core' as genes absent in 20% or less of the total element in this cluster.
  Whole_core <- ortho_counts[rowSums(ortho_counts < 1 ) <= ncol(ortho_counts)*0.2,]
  core_genes <- rownames(Whole_core)

  # Save results if 'core' genes were identified

  if(length(core_genes) > 0) {

    write.table(core_genes, file = paste("Core_genes-General",Cluster_number, ".txt", sep = ""),
                        sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(Whole_core, file = paste("Core_genes-General",Cluster_number, ".csv", sep = ""),
                        row.names = T, quote = F)
  }

  # Use a high height (h=0.8) to split the cluster into subclusters
  fit <- cutree(Cluster, h = 0.8)
  # Initialize vector of 'specific core' genes
  core_genes2 <- character()

  # Only look for specific core genes if the cluster actually splits into sub-clusters.
  if(length(unique(fit)) > 1) {
    for(SubClusterId in 1:length(unique(fit))) {
      # Extract the elements belonging to the current sub-cluster
      OrthoFinder_subcluster <- subset(ortho_counts, select = names(fit[grep(SubClusterId,fit)]))

      # Only analyze groups with more than 4 elements to ensure the "core" definition is meaningful.
      if(ncol(OrthoFinder_subcluster) > 4) {
        # Apply the same 20% absence threshold within this sub-cluster
        OrthoFinder_subcluster <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster < 1 ) <= ncol(OrthoFinder_subcluster)*0.2,]

        # Save results for this specific sub-cluster if 'core' genes was identified
        if(nrow(OrthoFinder_subcluster) > 0) {
          core <- rownames(OrthoFinder_subcluster)
          core_genes2 <- c(core_genes2,core)
          write.table(core, file = paste("Core_genes-Specific",SubClusterId,Cluster_number, ".txt", sep = ""),
                        sep = '\t', row.names = F, col.names = F, quote = F)
          write.csv(OrthoFinder_subcluster, file = paste("Core_genes-Specific",SubClusterId,Cluster_number, ".csv", sep = ""),
                        row.names = T, quote = F)
        }
      }
    }
  }

  return(list(
    general_core = core_genes,
    specific_core = core_genes2
    ))
}

analyze_individual_clusters <- function(
    Individual_clusters,
    transposase_OrthoFinder,
    fit_Orthofinder,
    seq_data,
    gene_data,
    blast_links
) {
  
  # Initialize the vector to track clusters that were successfully analyze
  Selected_clusters <- numeric()
  
  for (ClusterId in 1:Individual_clusters) {

    # Extract genomes belonging to the current cluster
    Cluster_matrix <- transposase_OrthoFinder[grep(ClusterId, fit_Orthofinder), ]
    
    # Logic to handle R's behavior of turning single-row matrices into vectors
    if (is.array(Cluster_matrix)) {
      # Remove orthogroups with zero-count in this cluster
      Cluster_matrix <- Cluster_matrix[, colSums(abs(Cluster_matrix)) > 0]
    } else {
      # Fallback for single-element clusters to prevent script failure
      Cluster_matrix <- matrix(, nrow = 1, ncol = 1)
      rownames(Cluster_matrix) <- names(fit_Orthofinder[grep(ClusterId,fit_Orthofinder)])
    }
    
    # Only proceed if we have enough data (orthogroups) to build a hierarchical tree
    if (ncol(Cluster_matrix) > 2) {
      
      # Build Hierarchical tree
      Cluster_distance <- dist(Cluster_matrix, method = 'binary')
      Cluster_Hier <- hclust(Cluster_distance, method = 'average')
      my_tree <- as.phylo(Cluster_Hier)
      
      # Save Tree and Cluster list, if multiple clusters wher identified
      if (Individual_clusters >= 2) {
        write.tree(phy = my_tree, file = paste("CargoHierarchicalTree_Cluster", ClusterId, ".nwk", sep = ""))
        write.table(Cluster_Hier$labels, file = paste("Cluster", ClusterId, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
      } else {
        write.tree(phy = my_tree, file = "CargoHierarchicalTree.nwk")
      }
      
      # Generated Heatmap of Orthogroups Presence/Absence
      Cluster_melt <- reshape2::melt(Cluster_matrix, as.is = T)
      Profile_plot <- ggplot(data = Cluster_melt, aes(x = reorder(Var2, value), y = Var1, fill = value)) +
        geom_tile() +
        scale_fill_viridis(option = "rocket", direction = -1) +
        labs(y = "element", x = "Orthogroup") +
        theme(axis.text.x = element_blank())


      # Save heatmap: modified dimensions based on the number of orthogroups and elements.
      plot_width <- min(49, ncol(Cluster_matrix) * 0.2 + 2)
      plot_height <- min(49, nrow(Cluster_matrix))
      if (Individual_clusters >= 2) {
        ggsave(Profile_plot, filename = paste("CargoHeatmap_Cluster", ClusterId, ".svg", sep = ""),
               width = plot_width, height = plot_height)
      } else {
        ggsave(Profile_plot, filename = "CargoHeatmap.svg",
               width = plot_width, height = plot_height)
      }
    } else {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "WARNING: There are issues with the following element: ", as.character(rownames(Cluster_matrix)), "\n", sep = ""))
    }
  }

  return(Selected_clusters)
}

plot_cluster_synteny <- function(
    Cluster_matrix, 
    orthogroup_table,
    core_genes,   
    seq_data,        
    gene_data,      
    blast_links,
    Cluster_number = ""
) {
  
  # Generate hierarchical cluster tree
  Cluster_distance <- dist(Cluster_matrix, method = 'binary')
  Cluster_Hier <- hclust(Cluster_distance, method = 'average')
  my_tree <- as.phylo(Cluster_Hier)
      
  # Data Preparation for Synteny Plot
  tree_sorted <- ladderize(my_tree, right = FALSE)
  selected_seqs <- row.names(Cluster_matrix)
      
  seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs)
  genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs)

  if( length(core_genes$specific_core) > 0 ) {
    # Extract Gene IDs belonging to the specific core orthogroups
    specific_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% core_genes$specific_core,]
    specific_vector <- unlist(strsplit(specific_table[,2],split = " "))
    specific_vector <- specific_vector[ specific_vector != ""]

    # Update the gene table
    genes_filtered <- genes_filtered %>%
      mutate(attribute = ifelse(ID %in% specific_vector, "specific", attribute))
  }

  if( length(core_genes$general_core) > 0 ) {
    # Extract Gene IDs belonging to the general core orthogroups
    general_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% core_genes$general_core,]
    general_vector <- unlist(strsplit(general_table[,2],split = " "))
    general_vector <- general_vector[ general_vector != ""]

    # Update the gene table: general core takes priority over 'specific' for coloring
    genes_filtered <- genes_filtered %>%
      mutate(attribute = ifelse(ID %in% general_vector, "general", attribute))
  }

  links_filtered <- blast_links %>%
    filter(seq_id %in% selected_seqs | seq_id2 %in% selected_seqs)
      
  # Align Sequence Order with Phylogenetic Tree
  p_tree_base <- ggtree::ggtree(tree_sorted, layout = "rectangular")
      
  # Determine sequence order based on tree tip position
  tree_y_coords <- p_tree_base$data %>%
    filter(isTip) %>%
    select(seq_id = label, y) %>%
    mutate(y = max(y) - y + 1) # Flip Y axis to match top-down order
      
  ordered_seqs <- seqs_filtered %>%
    inner_join(tree_y_coords, by = "seq_id") %>%
    arrange(y)
      
  max_seq_len <- max(ordered_seqs$length, na.rm = TRUE)

  p_tree <- p_tree_base +
  ggtree::geom_tiplab(align = TRUE, size = ifelse(nrow(Cluster_matrix) < 100, 2.5, 3), family = "mono") +
  ggtree::theme_tree2()

  p_genome <- gggenomes(
    seqs = ordered_seqs,
    links = links_filtered,
    genes = genes_filtered
    ) +
    geom_seq(aes(y = y)) +
    geom_gene(aes(y = y, fill = attribute), show.legend = T) +
    scale_fill_manual(name = "Core genes", values = c("general" = "red4", "specific" = "green4"), na.value = "cornsilk3", limits = c("general","specific")) + # Color genes based on the 'attribute'
    new_scale_fill() + # Use a second scale for the link (BLAST) identity
    geom_link(aes(y = y, fill = pident), colour = NA) +
    scale_fill_continuous(name = "Alignment Identity (%)") +
    geom_bin_label(aes(y = y), x = -10, size = ifelse(nrow(Cluster_matrix) < 50, 2.5, 3)) +
    scale_x_continuous(
     labels = label_number(accuracy = 1),
     limits = c(0, max_seq_len)) +
    scale_y_continuous(expand = expansion(mult = 1.0001 * (nrow(Cluster_matrix) ^ -1.1728))) # Y scale is dependant on number of sequences to match the three vizualitation
      
  # Combine plots
  final_plot <- p_tree + p_genome
      
  # Save Final Plot
  plot_width_final <- min(49, max(16, round(max_seq_len * 0.0001) + round(max(tree_sorted$edge.length, na.rm = TRUE) * 10)))
  plot_height_final <- min(49, nrow(Cluster_matrix))
      
  if ( Cluster_number != "") {
    ggsave(final_plot, filename = paste("CargoSynteny_Cluster", ClusterId, ".svg", sep = ""),
      width = plot_width_final, height = plot_height_final, limitsize = FALSE)
    # Record successful cluster
    Selected_clusters <- c(Selected_clusters, ClusterId)
  } else {
    ggsave(final_plot, filename = "CargoSynteny.svg",
      width = plot_width_final, height = plot_height_final, limitsize = FALSE)
  }
}

check_nesting <- function(
  ortho_counts,   
  seq_data,        
  gene_data,      
  blast_links) {

  # Iteration Through Every Element
  for(Element in 1:ncol(ortho_counts)) { 
    Element_name <- colnames(ortho_counts[Element])
    # Identify orthogroups present in the current element (non-zero rows)
    Genes_element <- rownames(ortho_counts[Element_name] %>%  filter(!if_all(everything(), ~ .x == 0)))

    # Create a reduced matrix containing only the orthogroups found in this Element
    reduced_orthofinder <- ortho_counts[Genes_element,]

    # Filter for other elements that share at least 80% of these orthogroups (Candidates for being 'containers' of this element)
    reduced_orthofinder2 <- reduced_orthofinder[,colSums(reduced_orthofinder < 1) < nrow(reduced_orthofinder) * 0.2]
    if(is.data.frame(reduced_orthofinder2)){
      # Remove the element itself from the comparison set
      reduced_orthofinder2 <- reduced_orthofinder2 %>% select(-all_of(Element_name))
      if(is.data.frame(reduced_orthofinder2)){
        for(comparison in 1:ncol(reduced_orthofinder2)) {
          Element_compare <- colnames(reduced_orthofinder2[comparison])
          Genes_compare <- rownames(ortho_counts[Element_compare] %>%  filter(!if_all(everything(), ~ .x == 0)))

          # A 'container' element should be significantly larger (1.5x more orthogroups)
          if(length(Genes_compare) > (length(Genes_element) * 1.5)) {

            # Extract links between the small element (query) and large element (subject)
            links_filtered2 <- blast_links %>%
            filter(qseqid %in% Element_name & sseqid %in% Element_compare) %>%
            select(seq_id = qseqid, start = qstart, end = qend,
                 seq_id2 = sseqid, start2 = sstart, end2 = send, pident)

            # Retrieve total lengths of both sequences
            seq_idLength <- as.numeric(seq_data[seq_data$seq_id %in% unique(links_filtered2$seq_id),2])
            seq_id2Length <- as.numeric(seq_data[seq_data$seq_id %in% unique(links_filtered2$seq_id2),2])

            # Check if the hits cover the small element from start to end (0.2 to 0.8 [allow some issues with the boundaries of the small element])
            is_nested_q <- any(links_filtered2$start < seq_idLength * 0.2) && any(links_filtered2$end > seq_idLength * 0.8)
            is_nested_s <- any(links_filtered2$start2 < seq_id2Length * 0.2) && any(links_filtered2$end2 > seq_id2Length * 0.8)

            if (xor(is_nested_q, is_nested_s)) {
              cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Element ", Element_name, " is nested in element ", Element_compare, "\n", sep = ""))

              selected_seq <- c(Element_name, Element_compare)
 
              seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seq)
              genes_filtered <- gene_data %>% filter(seq_id %in% selected_seq)
            
              links_filtered <- blast_links %>%
              filter(qseqid %in% selected_seq & sseqid %in% selected_seq) %>%
              select(seq_id = qseqid, start = qstart, end = qend,
                 seq_id2 = sseqid, start2 = sstart, end2 = send, pident)

              # Ordering sequences
              ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, selected_seq))
         
              # Plotting (Synteny map)
              p_genome <- gggenomes(
              seqs = ordered_seqs,
              links = links_filtered,
              genes = genes_filtered
              ) +
              geom_seq(aes(y = y)) +
              geom_gene(aes(y = y)) +
              geom_link(aes(y = y, fill = pident), colour = NA) +
              geom_bin_label(aes(y = y), x = -10) +
              scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
              theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
        
              # Saving the plot
              ggsave(p_genome, filename = paste("IndividualNestingEvent-", Element_name, "in",Element_compare, ".svg", sep = ""),
               width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)),
               height = min(49, length(selected_seq)), limitsize = FALSE)
            }
          }
        }

      }
    }
  }
}

gene_movement_analysis <- function(
    subcluster_number,  
    matrix,     
    ortho_counts, 
    core_genes,
    big=FALSE,   
    Cluster_number = ""
  ) {

  test_NbClust <- vector()

  # Define Comparison Groups: Elements in the subcluster that are in the row; 'other' elements are potential donors/recipients from different subcluster
  if(big){
    elements_cluster <- rownames(matrix)
    elements_other <- colnames(matrix)
  } else {
    elements_cluster <- rownames(matrix)
    elements_other <- tail(colnames(matrix),n=1)
  }

  selected_seqs2 <- c(elements_cluster,elements_other)

  # Filter Orthogroup Profiles: Isolate genes present in the subcluster vs genes present in the 'other' group
  OrthoFinder_subcluster_cluster <- subset(ortho_counts, select = elements_cluster) %>%
    filter(!if_all(everything(), ~ .x == 0))
  OrthoFinder_subcluster_other <- subset(ortho_counts, select = elements_other) %>%
    filter(!if_all(everything(), ~ .x == 0))
  #OrthoFinder_subcluster_cluster <- subset(ortho_counts, select = elements_cluster)
  #OrthoFinder_subcluster_cluster <- OrthoFinder_subcluster_cluster %>%  filter(!if_all(everything(), ~ .x == 0))
  #OrthoFinder_subcluster_other <- subset(ortho_counts, select = elements_other)
  #OrthoFinder_subcluster_other <- OrthoFinder_subcluster_other %>%  filter(!if_all(everything(), ~ .x == 0))

  # Movement = Genes in both groups MINUS the core genes of the cluster
  Gene_movement <- intersect(rownames(OrthoFinder_subcluster_cluster),rownames(OrthoFinder_subcluster_other))
  Nesting_test <- Gene_movement # Save raw intersection for nesting check
  Gene_movement <- Gene_movement[ ! Gene_movement %in% core_genes$general_core]
  
  # Logic Branching: Nesting vs. Simple Movement vs. Multiple Events
  # CASE A: Nesting Event (High overlap of total gene content)
  if(length(Nesting_test) >= nrow(OrthoFinder_subcluster_cluster)*0.8 || length(Nesting_test) >= nrow(OrthoFinder_subcluster_other)*0.8) {
    single_movements <- TRUE
    multiple_movements <- FALSE
    write.table(selected_seqs2, file = paste(Cluster_number,"SubCluster", subcluster_number, "_nestingEvent.txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
  # CASE B: Movement with Multiple Potential Origins. Meaning that different genes belong to different putative movement donod/acceptors
  } else if(length(Gene_movement) > 0 && big) {
    single_movements <- TRUE
    OrthoFinder_subcluster <- ortho_counts[Gene_movement,selected_seqs2]
    OrthoFinder_subcluster <- OrthoFinder_subcluster[,colSums(OrthoFinder_subcluster) > 0]

    # Refine element lists based on those actually carrying the moving genes
    elements_cluster <- elements_cluster[elements_cluster %in% colnames(OrthoFinder_subcluster)]
    elements_other <- elements_other[elements_other %in% colnames(OrthoFinder_subcluster)]

    transposase_OrthoFinder <- t(OrthoFinder_subcluster)

    # Use NbClust to see if 'moving' genes cluster into distinct groups
    if(nrow(transposase_OrthoFinder) > 3 && ncol(transposase_OrthoFinder) > 3 ){
      distance_mat <- dist(transposase_OrthoFinder, method='binary')
      Hierar_cl <- hclust(distance_mat, method = 'average')
      # Determine optimal number of clusters via the "Ball" index
      test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(transposase_OrthoFinder)-1), index = "ball")
 
      if(length(unique(test_NbClust$Best.partition)) > 1) {
        multiple_movements <- TRUE
      } else {
        multiple_movements <- FALSE
      }
    } else {
      multiple_movements <- FALSE
    }
    # Save results
    write.table(Gene_movement, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologs.txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(OrthoFinder_subcluster, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologsTable.csv", sep = ""),
                row.names = T, quote = F)

  # CASE C: Simple Single Movement
  } else if(length(Gene_movement) > 0) {
    OrthoFinder_subcluster <- ortho_counts[Gene_movement,selected_seqs2]
    single_movements <- TRUE
    multiple_movements <- FALSE

    # Save results
    write.table(Gene_movement, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologs.txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(OrthoFinder_subcluster, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologsTable.csv", sep = ""),
                row.names = T, quote = F)
  } else {
    single_movements <- FALSE
    multiple_movements <- FALSE
  }

  # Return Results
  return(list(
    movement = single_movements,
    multiple = multiple_movements,
    orthogroups= Gene_movement,
    elements_in = elements_cluster,
    elements_out = elements_other,
    K_analysis = test_NbClust
    ))
}

plot_subcluster_synteny <- function(
    subcluster_number,
    orthogroup_table,  
    dist_matrix,    
    cluster_fit,    
    ortho_counts, 
    core_genes,   
    seq_data,        
    gene_data,      
    blast_links,
    Cluster_number = ""
) {
  
  # Loop through the number of requested subclusters
  for (ClusterId in 1:subcluster_number) {
    
    # Extract the distance matrix rows for the current cluster
    matrix <- dist_matrix[grep(ClusterId, cluster_fit), ]
    Cluster_elements <- rownames(matrix)
    
    # Handle the case where matrix might be a vector (not a matrix)
    if (is.null(Cluster_elements)) {
      Cluster_elements <- names(grep(ClusterId, cluster_fit, value = TRUE))
    }
    
    # Write the cluster elements to a file
    write.table(Cluster_elements, file = paste(Cluster_number,"SubCluster", ClusterId, ".txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
    
    # Initialize flag to track if we successfully selected sequences for plotting
    selected_seqs2_defined <- FALSE
    multiple_movements <- FALSE
    
    # Main Cluster Processing Logic
    if (is.matrix(matrix)) {
      # First reduction step
      reduce_matrix <- matrix[, colSums(matrix) < nrow(matrix) * 0.98]
      
      # Block 1: There at least two sequence related to this subcluster
      if (ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
        reduce_matrix2 <- reduce_matrix[, (nrow(reduce_matrix) + 1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2) * 0.98, ]
        
        if (is.matrix(reduce_matrix3)) {

          reduce_matrix3 <- reduce_matrix3[, colSums(reduce_matrix3 < 1) > 1]
          reduce_matrix3 <- reduce_matrix3[, names(sort(colSums(reduce_matrix3), decreasing = F))]
          
          # OrthoFinder analysis
          Movement <- gene_movement_analysis(
          subcluster_number = ClusterId,  
          matrix = reduce_matrix3,     
          ortho_counts = ortho_counts, 
          core_genes = core_genes,
          big = TRUE,   
          Cluster_number = Cluster_number)

          selected_seqs2 <- c(Movement$elements_in,Movement$elements_out)
          selected_seqs2_defined <- Movement$movement
          multiple_movements <- Movement$multiple
          elements_cluster <- Movement$elements_in
          elements_other <- Movement$elements_out
          test_NbClust <- Movement$K_analysis
        }
        
      # Block 2: There's only one sequence related
      } else if (ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
        reduce_matrix <- reduce_matrix[order(reduce_matrix[,ncol(reduce_matrix)], decreasing = TRUE),]

        Movement <- gene_movement_analysis(
          subcluster_number = ClusterId,  
          matrix = reduce_matrix,     
          ortho_counts = ortho_counts, 
          core_genes = core_genes,   
          Cluster_number = Cluster_number)

          selected_seqs2 <- c(Movement$elements_in,Movement$elements_out)
          selected_seqs2_defined <- Movement$movement
        
      # Block 3: Fallback using a 0.99 threshold
      } else {
        reduce_matrix <- matrix[, colSums(matrix) < nrow(matrix) * 0.99]
        
        if (ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
          reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
          reduce_matrix2 <- reduce_matrix[, (nrow(reduce_matrix) + 1):ncol(reduce_matrix)]
          reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2) * 0.99, ]
          
          if (is.matrix(reduce_matrix3)) {

            reduce_matrix3 <- reduce_matrix3[, colSums(reduce_matrix3 < 1) > 1]
            reduce_matrix3 <- reduce_matrix3[, names(sort(colSums(reduce_matrix3), decreasing = F))]
            
            # OrthoFinder analysis
            gene_movement_analysis(
            subcluster_number = ClusterId,  
            matrix = reduce_matrix3,     
            ortho_counts = ortho_counts, 
            core_genes = core_genes, 
            big = TRUE,    
            Cluster_number = Cluster_number)

            selected_seqs2 <- c(Movement$elements_in,Movement$elements_out)
            selected_seqs2_defined <- Movement$movement
            multiple_movements <- Movement$multiple
            elements_cluster <- Movement$elements_in
            elements_other <- Movement$elements_out
            test_NbClust <- Movement$K_analysis
          }  
        } else if (ncol(reduce_matrix) > nrow(reduce_matrix)) {
          reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
          reduce_matrix <- reduce_matrix[order(reduce_matrix[,ncol(reduce_matrix)], decreasing = TRUE),]
          
          Movement <- gene_movement_analysis(
          subcluster_number = ClusterId,  
          matrix = reduce_matrix,     
          ortho_counts = ortho_counts, 
          core_genes = core_genes,   
          Cluster_number = Cluster_number)

          selected_seqs2 <- c(Movement$elements_in,Movement$elements_out)
          selected_seqs2_defined <- Movement$movement
        }
      }
      
      # Data preparation and Plotting (Runs only if sequences were selected)
      if (selected_seqs2_defined & ! multiple_movements) {
        
        # Filtering data
        seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs2)
        genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs2)

        if( length(core_genes$specific_core) > 0 ) {
          # Extract Gene IDs belonging to the specific core orthogroups
          specific_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% core_genes$specific_core,]
          specific_vector <- unlist(strsplit(specific_table[,2],split = " "))
          specific_vector <- specific_vector[ specific_vector != ""]

          genes_filtered <- genes_filtered %>%
            mutate(attribute = ifelse(ID %in% specific_vector, "specific", attribute))
        }

        # Extract Gene IDs belonging to the movement orthogroups
        movement_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% Movement$orthogroups,]
        movement_vector <- unlist(strsplit(movement_table[,2],split = " "))
        movement_vector <- movement_vector[ movement_vector != ""]

        genes_filtered <- genes_filtered %>%
          mutate(attribute = ifelse(ID %in% movement_vector, "Movement associated", attribute))

        if( length(core_genes$general_core) > 0 ) {
          # Extract Gene IDs belonging to the general core orthogroups
          general_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% core_genes$general_core,]
          general_vector <- unlist(strsplit(general_table[,2],split = " "))
          general_vector <- general_vector[ general_vector != ""]

          genes_filtered <- genes_filtered %>%
            mutate(attribute = ifelse(ID %in% general_vector, "general", attribute))
        }

        links_filtered <- blast_links %>%
          filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(seq_id = qseqid, start = qstart, end = qend,
                 seq_id2 = sseqid, start2 = sstart, end2 = send, pident)
        
        # Ordering sequences
        ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, selected_seqs2))
        
        # Plotting (Synteny map)
        p_genome <- gggenomes(
          seqs = ordered_seqs,
          links = links_filtered,
          genes = genes_filtered
        ) +
          geom_seq(aes(y = y)) +
          geom_gene(aes(y = y, fill = attribute), show.legend = T) +
          scale_fill_manual(name = "Core genes", values = c("general" = "red4", "specific" = "green4", "Movement associated" = "darkorchid"), na.value = "cornsilk3", limits = c("general","specific", "Movement associated")) +
          new_scale_fill() +
          geom_link(aes(y = y, fill = pident), colour = NA) +
          scale_fill_continuous(name = "Alignment Identity (%)") +
          geom_bin_label(aes(y = y), x = -10) +
          scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
        
        # Saving the plot
        ggsave(p_genome, filename = paste(Cluster_number,"CargoSynteny_SubCluster", ClusterId, ".svg", sep = ""),
               width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)),
               height = min(49, length(selected_seqs2)), limitsize = FALSE)
      # Data preparation and Plotting when multiple origins of movement where identify
      } else if(selected_seqs2_defined & multiple_movements) {

        k_out_cluster <- unique(test_NbClust$Best.partition[elements_other])
        k_cluster <- unique(test_NbClust$Best.partition[elements_cluster])

        # Plot each comparison
        for(i in k_cluster) {
          for(j in k_out_cluster) {
              cluster_seqs <- names(test_NbClust$Best.partition[grep(i,test_NbClust$Best.partition)][elements_cluster])
              cluster_seqs <- cluster_seqs[!is.na(cluster_seqs)]
              other_seqs <- names(test_NbClust$Best.partition[grep(j,test_NbClust$Best.partition)][elements_other])
              other_seqs <- other_seqs[!is.na(other_seqs)]

              OrthoFinder_subcluster_cluster <- subset(ortho_counts, select = cluster_seqs)
              OrthoFinder_subcluster_cluster <- OrthoFinder_subcluster_cluster %>%  filter(!if_all(everything(), ~ .x == 0))
              OrthoFinder_subcluster_other <- subset(ortho_counts, select = other_seqs)
              OrthoFinder_subcluster_other <- OrthoFinder_subcluster_other %>%  filter(!if_all(everything(), ~ .x == 0))

              Gene_movement <- intersect(rownames(OrthoFinder_subcluster_cluster),rownames(OrthoFinder_subcluster_other))

              selected_seqs2 <- c(names(test_NbClust$Best.partition[grep(i,test_NbClust$Best.partition)][elements_cluster]),names(test_NbClust$Best.partition[grep(j,test_NbClust$Best.partition)][elements_other]))
              selected_seqs2 <- selected_seqs2[!is.na(selected_seqs2)]

              OrthoFinder_subcluster <- ortho_counts[Gene_movement,selected_seqs2]
              OrthoFinder_subcluster <- OrthoFinder_subcluster[,colSums(OrthoFinder_subcluster) > 0]

              Gene_movement <- rownames(OrthoFinder_subcluster)

              selected_seq2 <- colnames(OrthoFinder_subcluster)

              if(length(selected_seq2) >= 2) {
                # Filtering data
                seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs2)
                genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs2)

                if( length(core_genes$specific_core) > 0 ) {
                  # Extract Gene IDs belonging to the specific core orthogroups
                  core_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% core_genes$specific_core,]
                  core_vector <- unlist(strsplit(core_table[,2],split = " "))
                  core_vector <- core_vector[ core_vector != ""]

                  genes_filtered <- genes_filtered %>%
                  mutate(attribute = ifelse(ID %in% core_vector, "specific", attribute))
                }

                if( length(core_genes$general_core) > 0 ) {
                  # Extract Gene IDs belonging to the general core orthogroups
                  general_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% core_genes$general_core,]
                  general_vector <- unlist(strsplit(general_table[,2],split = " "))
                  general_vector <- general_vector[ general_vector != ""]
 
                  genes_filtered <- genes_filtered %>%
                  mutate(attribute = ifelse(ID %in% general_vector, "general", attribute))
                }

                # Extract Gene IDs belonging to the movement orthogroups
                movement_table <- orthogroup_table[ orthogroup_table$Orthogroup %in% Gene_movement,]
                movement_vector <- unlist(strsplit(movement_table[,2],split = " "))
                movement_vector <- movement_vector[ movement_vector != ""]
 
                genes_filtered <- genes_filtered %>%
                  mutate(attribute = ifelse(ID %in% movement_vector, "Movement associated", attribute))

                links_filtered <- blast_links %>%
                filter(qseqid %in% selected_seqs2 & sseqid %in% selected_seqs2) %>%
                select(seq_id = qseqid, start = qstart, end = qend,
                 seq_id2 = sseqid, start2 = sstart, end2 = send, pident)
        
                # Ordering sequences
                ordered_seqs <- seqs_filtered %>% arrange(match(seq_id, selected_seqs2))

                if(nrow(links_filtered) >= 1){
                  # Plotting (Synteny map)
                  p_genome <- gggenomes(
                  seqs = ordered_seqs,
                  links = links_filtered,
                  genes = genes_filtered
                  ) +
                  geom_seq(aes(y = y)) +
                  geom_gene(aes(y = y, fill = attribute), show.legend = T) +
                  scale_fill_manual(name = "Core genes", values = c("general" = "red4", "specific" = "green4", "Movement associated" = "darkorchid"), na.value = "cornsilk3", limits = c("general","specific", "Movement associated")) +
                  new_scale_fill() +
                  geom_link(aes(y = y, fill = pident), colour = NA) +
                  scale_fill_continuous(name = "Alignment Identity (%)") +
                  geom_bin_label(aes(y = y), x = -10) +
                  scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
                  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
        
                  # Saving the plot
                  ggsave(p_genome, filename = paste(Cluster_number,"CargoSynteny_SubCluster", ClusterId, "-", i,"vs",j,".svg", sep = ""),
                  width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)),
                  height = min(49, length(selected_seqs2)), limitsize = FALSE)
                }
              }
          }
        }
      } else {
        cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Subcluster ", ClusterId, " can not be analyzed automatically", "\n", sep = ""))
      }
      
    } else {
      # Catches single-element clusters or other matrix issues
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Subcluster ", ClusterId, " can not be analyzed automatically", "\n", sep = ""))
    }
  } # End of for loop
}

# Start the process

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Reading data...","\n", sep=""))

data_list <- load_and_preprocess_data(min_pident = arguments$pident)

# Preprocess data
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Identifying if there are separate clusters","\n", sep=""))

clustering_data <- perform_initial_clustering()

cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Checking for nesting events ","\n", sep = ""))

check_nesting(ortho_counts = clustering_data$orthofinder_counts,
    seq_data = data_list$seqs,
    gene_data = data_list$genes,
    blast_links = data_list$blast_results)

# Process Clusters
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing ",clustering_data$Individual_clusters," individual cluster(s) identified","\n", sep=""))

# Process each cluster
Selected_Cluster_ID <- analyze_individual_clusters(
    Individual_clusters = clustering_data$Individual_clusters,
    transposase_OrthoFinder = clustering_data$transposed_counts,
    fit_Orthofinder = clustering_data$cluster_fit,
    seq_data = data_list$seqs,
    gene_data = data_list$genes,
    blast_links = data_list$blast_results)

if((clustering_data$Individual_clusters == 1) & (arguments$subclusters >= 2)){
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))
  
  core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl)

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Generating main plot","\n", sep=""))

  plot_cluster_synteny(Cluster_matrix=clustering_data$transposed_counts,orthogroup_table=data_list$orthogroups_table,core_genes=core_genes,seq_data=data_list$seqs,gene_data=data_list$genes,blast_links=data_list$links)

  if(arguments$captainRemoval <= 1) {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing ",arguments$subclusters," subclusters","\n", sep=""))

    fit <- cutree(clustering_data$Hierar_cl, k = arguments$subclusters) # k is the number of subclusters from the first approach
    New_dist <- as.matrix(clustering_data$distance_mat)

    plot_subcluster_synteny(
    subcluster_number = arguments$subclusters,
    orthogroup_table = data_list$orthogroups_table,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts, 
    core_genes = core_genes,   
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results)
  } else if(arguments$captainRemoval >= 2) {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))

    test_NbClust <- NbClust(clustering_data$distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(clustering_data$distance_mat)-1), index = "ball")
    test_fit <- cutree(clustering_data$Hierar_cl, h = 0.6)
    K_number <- round(mean(c(length(unique(test_NbClust$Best.partition)), length(unique(fit)))))
    fit <- cutree(clustering_data$Hierar_cl, k = K_number)

    New_dist <- as.matrix(clustering_data$distance_mat)

    plot_subcluster_synteny(
    subcluster_number = K_number,
    orthogroup_table = data_list$orthogroups_table,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts, 
    core_genes = core_genes,   
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results)
  }

} else if((clustering_data$Individual_clusters > 1) & (arguments$subclusters >= 2) & length(Selected_Cluster_ID) == 1) {

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))

  New_ortho_counts <- clustering_data$orthofinder_counts[,grep(Selected_Cluster_ID,clustering_data$cluster_fit)]
  Cluster_matrix <- clustering_data$transposed_counts[grep(Selected_Cluster_ID,clustering_data$cluster_fit),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  distance_mat <- dist(Cluster_matrix, method='binary')
  Hierar_cl <- hclust(distance_mat, method = 'average')
  
  core_genes <- core_genes_analysis(ortho_counts = New_ortho_counts, Cluster = Hierar_cl,
    Cluster_number = paste("_Cluster", Selected_Cluster_ID, sep=""))

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Generating main plot","\n", sep=""))

  plot_cluster_synteny(Cluster_matrix=Cluster_matrix,orthogroup_table=data_list$orthogroups_table,core_genes=core_genes,seq_data=data_list$seqs,gene_data=data_list$genes,blast_links=data_list$links, Cluster_number="1")

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))
  
  Cluster_matrix <- clustering_data$transposed_counts[grep(Selected_Cluster_ID,clustering_data$cluster_fit),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  distance_mat <- dist(Cluster_matrix, method='binary')
  Hierar_cl <- hclust(distance_mat, method = 'average')

  test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(distance_mat)-1), index = "ball")
  test_fit <- cutree(Hierar_cl, h = 0.6)
  K_number <- round(mean(c(length(unique(test_NbClust$Best.partition)), length(unique(test_fit)))))
  fit <- cutree(Hierar_cl, k = K_number)

  New_dist <- as.matrix(distance_mat)

  plot_subcluster_synteny(
    subcluster_number = K_number,
    orthogroup_table = data_list$orthogroups_table,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts,  
    core_genes = core_genes,  
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results,
    Cluster_number = paste("Cluster", Selected_Cluster_ID, sep=""))

} else if((clustering_data$Individual_clusters > 1) & (arguments$subclusters >= 2) & length(Selected_Cluster_ID) > 1) {

  for(MainClusterID in Selected_Cluster_ID ) {
    Cluster_matrix <- clustering_data$transposed_counts[grep(MainClusterID,clustering_data$cluster_fit),]

    if(! is.matrix(Cluster_matrix)) {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Cluster ", MainClusterID, " can not be analyzed automatically", "\n", sep = ""))
      next
    }

    Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]

    if(nrow(Cluster_matrix) < 4) {
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "Cluster ", MainClusterID, " can not be analyzed automatically", "\n", sep = ""))
      next
    }

    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes in Cluster ",MainClusterID,"\n", sep=""))

    New_ortho_counts <- clustering_data$orthofinder_counts[,grep(MainClusterID,clustering_data$cluster_fit)]
    distance_mat <- dist(Cluster_matrix, method='binary')
    Hierar_cl <- hclust(distance_mat, method = 'average')
  
    core_genes <- core_genes_analysis(ortho_counts = New_ortho_counts, Cluster = Hierar_cl,
      Cluster_number = paste("_Cluster", MainClusterID, sep=""))

    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Generating main plot","\n", sep=""))

    plot_cluster_synteny(Cluster_matrix=Cluster_matrix,orthogroup_table=data_list$orthogroups_table,core_genes=core_genes,seq_data=data_list$seqs,gene_data=data_list$genes,blast_links=data_list$links, Cluster_number=MainClusterID)

    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters in Cluster ",MainClusterID,"\n", sep=""))

    test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(Cluster_matrix)-1), index = "ball")
    test_fit <- cutree(Hierar_cl, h = 0.6)
    K_number <- round(mean(c(length(unique(test_NbClust$Best.partition)), length(unique(test_fit)))))
    fit <- cutree(Hierar_cl, k = K_number)

    New_dist <- as.matrix(distance_mat)

    plot_subcluster_synteny(
    subcluster_number = K_number,
    orthogroup_table = data_list$orthogroups_table,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts, 
    core_genes = core_genes,   
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results,
    Cluster_number = paste("Cluster", MainClusterID, "-", sep=""))

  } 
} else{
  cat(paste("    [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No subcluster to analyzed","\n", sep=""))
  
  if((clustering_data$Individual_clusters == 1) & (arguments$subclusters <= 1)) {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))
  
    core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl)

    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Generating plot","\n", sep=""))
    plot_cluster_synteny(Cluster_matrix=clustering_data$transposed_counts,orthogroup_table=data_list$orthogroups_table,core_genes=core_genes,seq_data=data_list$seqs,gene_data=data_list$genes,blast_links=data_list$links)
  }
} 

if(clustering_data$Individual_clusters == 1) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are differences between captain and cargo dendograms","\n", sep=""))

  phylo_tree <- read.tree("CaptainPhylogeny.nw")

  if(length(phylo_tree$tip.label) > 4){
    dist_matrix <- cophenetic.phylo(phylo_tree)
    min_dist <- min(dist_matrix)
    max_dist <- max(dist_matrix)
    normalized_dist_matrix <- (dist_matrix - min_dist) / (max_dist - min_dist)
    normalized_dist_matrix2 <- scale(normalized_dist_matrix)
    hclust_tree <- hclust(as.dist(normalized_dist_matrix), method = "average")
    hclust_tree$height <- sort(hclust_tree$height)

    test_NbClust <- NbClust(normalized_dist_matrix2, method = "average", min.nc = 1, max.nc = min(6, nrow(normalized_dist_matrix2)-1), index = "ball")
    test_fit <- cutree(hclust_tree, h = 0.5)
    K_number <- round(mean(c(length(unique(test_NbClust$Best.partition)),length(unique(test_fit)))))
    fit_tree <- cutree(hclust_tree, k = K_number)
    fit_tree <- fit_tree[order(names(fit_tree))]

    fit <- cutree(clustering_data$Hierar_cl, k = K_number) 
    fit <- fit[order(names(fit))]

    #First round of discordance - Cargo tree view
    contingency_table <- table(fit_tree, fit)
    best_match <- apply(contingency_table, 1, which.max)
    fit_aligned <- fit_tree
    for(i in 1:length(best_match)){
      original_label <- names(best_match)[i]
      new_label <- best_match[i]
      fit_aligned[fit_tree == original_label] <- as.numeric(new_label)
    }

    discordant_indices <- which(fit != fit_aligned)
    discordant_elements <- names(fit)[discordant_indices]

    #Second round of discordance - Captain tree view
    contingency_table <- table(fit, fit_tree)
    best_match <- apply(contingency_table, 1, which.max)
    fit_aligned <- fit
    for(i in 1:length(best_match)){
      original_label <- names(best_match)[i]
      new_label <- best_match[i]
      fit_aligned[fit == original_label] <- as.numeric(new_label)
    }
    discordant_indices <- which(fit_tree != fit_aligned)
    discordant_elements <- unique(c(discordant_elements,names(fit_tree)[discordant_indices]))

    if(length(discordant_elements) >= 1) {
      cat(paste("    [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","There ",length(discordant_elements)," discordant element(s) identified","\n", sep=""))
      write.table(discordant_elements, file = "Discordant_elements.txt", sep = '\t', row.names = F, col.names= F,quote = F)
    } else {
      cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No discordance detected","\n", sep=""))
    }

    dend_list <- dendlist(clustering_data$Hierar_cl, hclust_tree)
    svg(filename="CaptainVsCargo_tanglegram.svg", width = 16, height = min(49, length(fit)/2))
    dend_list %>% ladderize %>% 
    untangle(method = "step1side", k_seq = K_number:(length(fit)-1)) %>%
    set("branches_k_color", k=K_number) %>% 
    tanglegram(faster = TRUE, main_left = as.character("Cargo Cluster"),main_right = as.character("Captain tree"))
    invisible(dev.off())
  } else {
    dist_matrix <- cophenetic.phylo(phylo_tree)
    min_dist <- min(dist_matrix)
    max_dist <- max(dist_matrix)
    normalized_dist_matrix <- (dist_matrix - min_dist) / (max_dist - min_dist)
    normalized_dist_matrix2 <- scale(normalized_dist_matrix)
    hclust_tree <- hclust(as.dist(normalized_dist_matrix), method = "average")

    fit_tree <- cutree(hclust_tree, k = 2)
    fit_tree <- fit_tree[order(names(fit_tree))]

    fit <- cutree(clustering_data$Hierar_cl, k = 2) 
    fit <- fit[order(names(fit))]

    #First round of discordance - Cargo tree view
    contingency_table <- table(fit_tree, fit)
    best_match <- apply(contingency_table, 1, which.max)
    fit_aligned <- fit_tree
    for(i in 1:length(best_match)){
      original_label <- names(best_match)[i]
      new_label <- best_match[i]
      fit_aligned[fit_tree == original_label] <- as.numeric(new_label)
    }

    discordant_indices <- which(fit != fit_aligned)
    discordant_elements <- names(fit)[discordant_indices]

    #Second round of discordance - Captain tree view
    contingency_table <- table(fit, fit_tree)
    best_match <- apply(contingency_table, 1, which.max)
    fit_aligned <- fit
    for(i in 1:length(best_match)){
      original_label <- names(best_match)[i]
      new_label <- best_match[i]
      fit_aligned[fit == original_label] <- as.numeric(new_label)
    }

    discordant_indices <- which(fit_tree != fit_aligned)
    discordant_elements <- unique(c(discordant_elements,names(fit_tree)[discordant_indices]))

    if(length(discordant_elements) >= 1) {
      cat(paste("    [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","There ",length(discordant_elements)," discordant element(s) identified","\n", sep=""))
      write.table(discordant_elements, file = "Discordant_elements.txt", sep = '\t', row.names = F, col.names= F,quote = F)
    } else {
      cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No discordance detected","\n", sep=""))
    }

    dend_list <- dendlist(clustering_data$Hierar_cl, hclust_tree)
    svg(filename="CaptainVsCargo_tanglegram.svg", width = 16, height = min(49, length(fit)/2))
    dend_list %>% ladderize %>% 
    untangle(method = "step1side", k_seq = 2:(length(fit)-1)) %>%
    set("branches_k_color", k=2) %>% 
    tanglegram(faster = TRUE, main_left = as.character("Cargo Cluster"),main_right = as.character("Captain tree"))
    invisible(dev.off())
  } 
}
 
cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Finished","\n", sep=""))
