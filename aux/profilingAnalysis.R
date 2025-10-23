# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory [default %default]",metavar="PATH"),
  make_option(c("-s", "--subclusters"), type="integer", action = "store", default=1,
              help=" Number of subclusters get in the synteny analysis", metavar="number"),
  make_option(c("-c", "--captainRemoval"), type="integer", action = "store", default=0,
              help=" Number of elements removes from the original cluster", metavar="number")
)

# Parse the command-line arguments
arguments <- parse_args(OptionParser(option_list = option_list))

# Check for required arguments and flags
if (is.null(arguments$directory)) {
  stop("Error: directory must be provided.", call.=FALSE)
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

# Function block

load_and_preprocess_data <- function(
    blast_file = "clean_results.txt",
    gff_file = "Final_model.gff",
    min_pident = 70
) {
  
  # 1. Define BLAST Column Specifications
  blast_col_names <- c("qseqid", "sseqid", "qstart", "qend", "sstart", "send",
                       "pident", "length", "qlen", "slen")
  blast_col_classes <- c("character", "character", rep("numeric", 8))
  
  # 2. Read and Preprocess BLAST Results (blast_results)
  blast_results <- read.delim(
    blast_file,
    header = FALSE,
    sep = "\t",
    col.names = blast_col_names,
    colClasses = blast_col_classes
  ) %>%
    # Filter out self-hits and low-identity hits
    filter(qseqid != sseqid, pident > min_pident)
  
  # 3. Create Link Table (links)
  links <- blast_results %>%
    mutate(alignment_id = row_number()) %>%
    select(
      seq_id = qseqid, start = qstart, end = qend,
      seq_id2 = sseqid, start2 = sstart, end2 = send,
      pident, length
    )
  
  # 4. Create Sequence Length Table (seqs)
  seqs <- bind_rows(
    blast_results %>% select(seq_id = qseqid, length = qlen),
    blast_results %>% select(seq_id = sseqid, length = slen)
  ) %>%
    distinct() %>%
    group_by(seq_id) %>%
    summarize(length = max(length)) %>%
    ungroup()
  
  # 5. Define GFF Column Specifications
  gff_col_names <- c("seqname", "source", "feature", "start", "end", "score",
                     "strand", "frame", "attribute")
  gff_col_classes <- c(rep("character", 3), rep("numeric", 2), rep("character", 4))
  
  # 6. Read and Preprocess GFF Features (genes)
  genes <- read.delim(
    gff_file,
    header = FALSE,
    comment.char = "#",
    sep = "\t",
    col.names = gff_col_names,
    colClasses = gff_col_classes
  ) %>%
    filter(feature == "gene") %>%
    select(seq_id = seqname, start, end, strand) %>%
    mutate(
      seq_id = as.character(seq_id),
      type = "CDS"
    )
  
  # 7. Return all processed data frames as a list
  return(list(
    blast_results = blast_results,
    links = links,
    seqs = seqs,
    genes = genes
  ))
}

perform_initial_clustering <- function(
    orthogroups_file = "Orthogroups.GeneCount.tsv"
) {
  
  # 1. Read OrthoFinder Gene Counts
  OrthoFinder <- read.table(
    orthogroups_file,
    header = TRUE,
    check.names = FALSE,
    row.names = 1
  )
  
  # 2. Preprocess: Remove the "Total" column
  OrthoFinder2 <- OrthoFinder[, !names(OrthoFinder) %in% c("Total")]
  
  # 3. Transpose the matrix
  transposed_counts <- t(OrthoFinder2)
  
  # 4. Calculate Distance Matrix
  distance_mat <- dist(transposed_counts, method = 'binary')
  
  # 5. Perform Hierarchical Clustering
  Hierar_cl <- hclust(distance_mat, method = 'average')
  
  # 6. Identify Initial Number of Clusters
  fit_Orthofinder <- cutree(Hierar_cl, h = 0.99999999)
  Individual_clusters <- length(unique(fit_Orthofinder))
  
  # 7. Return all processed objects as a list
  return(list(
    transposed_counts = transposed_counts,
    orthofinder_counts = OrthoFinder2,
    distance_mat = distance_mat,
    Hierar_cl = Hierar_cl,
    cluster_fit = fit_Orthofinder,
    Individual_clusters = Individual_clusters
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
  
  # Initialize the vector to track clusters that were successfully plotted
  Selected_clusters <- numeric()
  
  for (ClusterId in 1:Individual_clusters) {
    
    # 1. Extract Cluster Matrix and Handle Single/Zero-Element Case
    Cluster_matrix <- transposase_OrthoFinder[grep(ClusterId, fit_Orthofinder), ]
    
    if (is.array(Cluster_matrix)) {
      # Remove orthogroups with zero counts in this cluster
      Cluster_matrix <- Cluster_matrix[, colSums(abs(Cluster_matrix)) > 0]
    } else {
      # Handle case where cluster might have one element by creating a false empty matrix
      Cluster_matrix <- matrix(, nrow = 1, ncol = 1)
      rownames(Cluster_matrix) <- names(fit_Orthofinder[grep(2,fit_Orthofinder)])
    }
    
    # 2. Check for sufficient orthogroups for analysis
    if (ncol(Cluster_matrix) > 2) {
      
      # 3. Create Hierarchical tree
      Cluster_distance <- dist(Cluster_matrix, method = 'binary')
      Cluster_Hier <- hclust(Cluster_distance, method = 'average')
      my_tree <- as.phylo(Cluster_Hier)
      
      # Save Tree and Cluster List
      if (Individual_clusters >= 2) {
        write.tree(phy = my_tree, file = paste("CargoHC_Cluster", ClusterId, ".nwk", sep = ""))
        write.table(Cluster_Hier$labels, file = paste("Cluster", ClusterId, ".txt", sep = ""), sep = '\t', row.names = F, col.names = F, quote = F)
      } else {
        write.tree(phy = my_tree, file = "CargoHC.nwk")
      }
      
      # 4. Create Profile Plot (Heatmap)
      Cluster_melt <- reshape2::melt(Cluster_matrix, as.is = T)
      Profile_plot <- ggplot(data = Cluster_melt, aes(x = reorder(Var2, value), y = Var1, fill = value)) +
        geom_tile() +
        scale_fill_viridis(option = "rocket", direction = -1) +
        labs(y = "element", x = "Orthogroup") +
        theme(axis.text.x = element_blank())
      
      # Save Profile Plot
      plot_width <- min(49, ncol(Cluster_matrix) * 0.2 + 2)
      plot_height <- min(49, nrow(Cluster_matrix))
      if (Individual_clusters >= 2) {
        ggsave(Profile_plot, filename = paste("CargoProfiling_Cluster", ClusterId, ".svg", sep = ""),
               width = plot_width, height = plot_height)
      } else {
        ggsave(Profile_plot, filename = "CargoProfiling.svg",
               width = plot_width, height = plot_height)
      }
      
      # 5. Data Preparation for Synteny Plot
      tree_sorted <- ladderize(my_tree, right = FALSE)
      selected_seqs <- row.names(Cluster_matrix)
      
      seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs)
      genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs)
      
      links_filtered <- blast_links %>%
        filter(qseqid %in% selected_seqs | sseqid %in% selected_seqs) %>%
        select(
          seq_id = qseqid, start = qstart, end = qend,
          seq_id2 = sseqid, start2 = sstart, end2 = send, pident
        )
      
      # 6. Align Sequence Order with Phylogenetic Tree
      p_tree_base <- ggtree::ggtree(tree_sorted, hang = 0.05) +
        theme(plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"))
      tree_data <- p_tree_base$data
      
      # Determine sequence order based on tree tip position
      tree_y_coords <- tree_data %>%
        filter(isTip) %>%
        select(seq_id = label, y) %>%
        # Invert and shift y-coordinates to match gggenomes' top-down plotting order
        mutate(y = max(y) - y + 1)
      
      ordered_seqs <- seqs_filtered %>%
        inner_join(tree_y_coords, by = "seq_id") %>%
        arrange(y)
      
      # 7. Create Phylogenetic and Genomic Plots
      p_tree <- p_tree_base +
        ggtree::geom_tiplab(align = TRUE, size = 3) +
        ggtree::theme_tree2()
      
      max_seq_len <- max(ordered_seqs$length, na.rm = TRUE)
      
      p_genome <- gggenomes::gggenomes(
        seqs = ordered_seqs,
        links = links_filtered,
        genes = genes_filtered
      ) +
        gggenomes::geom_seq(aes(y = y)) +
        gggenomes::geom_gene(aes(y = y)) +
        gggenomes::geom_link(aes(y = y, fill = pident), colour = NA) +
        gggenomes::geom_bin_label(aes(y = y), x = -10) +
        scale_x_continuous(
          labels = label_number(accuracy = 1),
          limits = c(0, max_seq_len)
        ) +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
      
      # Combine plots (Requires the 'patchwork' package)
      final_plot <- p_tree + p_genome
      
      # 8. Save Final Plot
      plot_width_final <- min(49, max(16, round(max_seq_len * 0.0001) + round(max(tree_sorted$edge.length, na.rm = TRUE) * 10)))
      plot_height_final <- min(49, nrow(Cluster_matrix))
      
      if (Individual_clusters >= 2) {
        ggsave(final_plot, filename = paste("CargoSynteny_Cluster", ClusterId, ".svg", sep = ""),
               width = plot_width_final, height = plot_height_final, limitsize = FALSE)
        # Record successful cluster
        Selected_clusters <- c(Selected_clusters, ClusterId)
      } else {
        ggsave(final_plot, filename = "CargoSynteny.svg",
               width = plot_width_final, height = plot_height_final, limitsize = FALSE)
      }
      
    } else {
      # 9. Handle Clusters with 0, 1, or 2 Orthogroups
      cat(paste("  [", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", "WARNING: There are issues with the following element: ", as.character(rownames(Cluster_matrix)), "\n", sep = ""))
    }
  }
  
  return(Selected_clusters)
}

core_genes_analysis <- function(    
    ortho_counts,
    Cluster,
    Cluster_number = ""
) {
  #core_genes <- vector()
  core_genes2 <- vector()
  #identify the possibility of general core genes
  Whole_core <- ortho_counts[rowSums(ortho_counts < 1 ) <= ncol(ortho_counts)*0.2,]
  core_genes <- rownames(Whole_core)

  if(length(core_genes) > 0) {

    write.table(core_genes, file = paste("Core_genes-General",Cluster_number, ".txt", sep = ""),
                        sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(Whole_core, file = paste("Core_genes-General",Cluster_number, ".csv", sep = ""),
                        row.names = T, quote = F)
  }

  fit <- cutree(Cluster, h = 0.8)

  if(length(unique(fit)) > 1) {
    for(SubClusterId in 1:length(unique(fit))) {
      OrthoFinder_subcluster <- subset(ortho_counts, select = names(fit[grep(SubClusterId,fit)]))
      OrthoFinder_subcluster <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster < 1 ) <= ncol(OrthoFinder_subcluster)*0.2,]
      #OrthoFinder_subcluster <- OrthoFinder_subcluster %>%  filter(!if_any(everything(), ~ .x == 0))
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

  return(list(
    general_core = core_genes,
    specific_core = core_genes2
    ))
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
  #elements_cluster <- vector()
  #elements_other <- vector()

  if(big){
    elements_cluster <- rownames(matrix)
    elements_other <- colnames(matrix)
  } else {
    elements_cluster <- rownames(matrix)
    elements_other <- tail(colnames(matrix),n=1)
  }

  selected_seqs2 <- c(elements_cluster,elements_other)

  OrthoFinder_subcluster_cluster <- subset(ortho_counts, select = elements_cluster)
  #OrthoFinder_subcluster_cluster <- OrthoFinder_subcluster_cluster[apply(OrthoFinder_subcluster_cluster!=0, 1, any),]
  OrthoFinder_subcluster_cluster <- OrthoFinder_subcluster_cluster %>%  filter(!if_all(everything(), ~ .x == 0))
  OrthoFinder_subcluster_other <- subset(ortho_counts, select = elements_other)
  OrthoFinder_subcluster_other <- OrthoFinder_subcluster_other %>%  filter(!if_all(everything(), ~ .x == 0))

  Gene_movement <- intersect(rownames(OrthoFinder_subcluster_cluster),rownames(OrthoFinder_subcluster_other))
  Nesting_test <- Gene_movement[ ! Gene_movement %in% core_genes$general_core]
  Gene_movement <- Gene_movement[ ! Gene_movement %in% c(core_genes$general_core, core_genes$specific_core)]
  

  if(length(Nesting_test) >= nrow(OrthoFinder_subcluster_cluster)*0.8 || length(Nesting_test) >= nrow(OrthoFinder_subcluster_other)*0.8) {
    single_movements <- TRUE
    multiple_movements <- FALSE
    write.table(selected_seqs2, file = paste(Cluster_number,"SubCluster", subcluster_number, "_nestingEvent.txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
  } else if(length(Gene_movement) > 0 && big) {
    single_movements <- TRUE
    OrthoFinder_subcluster <- ortho_counts[Gene_movement,selected_seqs2]
    OrthoFinder_subcluster <- OrthoFinder_subcluster[,colSums(OrthoFinder_subcluster) > 0]

    elements_cluster <- elements_cluster[elements_cluster %in% colnames(OrthoFinder_subcluster)]
    elements_other <- elements_other[elements_other %in% colnames(OrthoFinder_subcluster)]

    transposase_OrthoFinder <- t(OrthoFinder_subcluster)
    if(nrow(transposase_OrthoFinder) > 3){
      distance_mat <- dist(transposase_OrthoFinder, method='binary')
      Hierar_cl <- hclust(distance_mat, method = 'average')
      test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(transposase_OrthoFinder)-1), index = "ball")
 
      if(length(unique(test_NbClust$Best.partition)) > 1) {
        multiple_movements <- TRUE
      } else {
        multiple_movements <- FALSE
      }
    } else {
      multiple_movements <- FALSE
    }
    write.table(Gene_movement, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologs.txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(OrthoFinder_subcluster, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologsTable.csv", sep = ""),
                row.names = T, quote = F)
  } else if(length(Gene_movement) > 0) {
    OrthoFinder_subcluster <- ortho_counts[Gene_movement,selected_seqs2]
    single_movements <- TRUE
    multiple_movements <- FALSE
    write.table(Gene_movement, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologs.txt", sep = ""),
                sep = '\t', row.names = F, col.names = F, quote = F)
    write.csv(OrthoFinder_subcluster, file = paste(Cluster_number,"SubCluster", subcluster_number, "_moveOrthologsTable.csv", sep = ""),
                row.names = T, quote = F)
  } else {
    single_movements <- FALSE
    multiple_movements <- FALSE
  }
  return(list(
    movement = single_movements,
    multiple = multiple_movements,
    elements_in = elements_cluster,
    elements_out = elements_other,
    K_analysis = test_NbClust
    ))
}

plot_subcluster_synteny <- function(
    subcluster_number,  
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
    
    # --- Main Cluster Processing Logic ---
    if (is.matrix(matrix)) {
      # First reduction step
      reduce_matrix <- matrix[, colSums(matrix) < nrow(matrix) * 0.97]
      
      # Block 1: Enough sequences for detailed reduction (>= nrow + 2)
      if (ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
        reduce_matrix2 <- reduce_matrix[, (nrow(reduce_matrix) + 1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2) * 0.97, ]
        reduce_matrix3 <- reduce_matrix3[, colSums(reduce_matrix3 < 1) > 1]
        
        if (is.matrix(reduce_matrix3)) {
          
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
        
      # Block 2: Simple reduction (ncol > nrow)
      } else if (ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]

        Movement <- gene_movement_analysis(
          subcluster_number = ClusterId,  
          matrix = reduce_matrix,     
          ortho_counts = ortho_counts, 
          core_genes = core_genes,   
          Cluster_number = Cluster_number)

          selected_seqs2 <- c(Movement$elements_in,Movement$elements_out)
          selected_seqs2_defined <- Movement$movement
        
      # Block 3: Fallback using a 0.98 threshold
      } else {
        reduce_matrix <- matrix[, colSums(matrix) < nrow(matrix) * 0.98]
        
        if (ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
          reduce_matrix <- reduce_matrix[, names(sort(colSums(reduce_matrix), decreasing = F))]
          reduce_matrix2 <- reduce_matrix[, (nrow(reduce_matrix) + 1):ncol(reduce_matrix)]
          reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2) * 0.98, ]
          reduce_matrix3 <- reduce_matrix3[, colSums(reduce_matrix3 < 1) > 1]
          
          if (is.matrix(reduce_matrix3)) {
            selected_seqs2 <- c(rownames(reduce_matrix3), colnames(reduce_matrix3))
            
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
      
      # 4. Data preparation and Plotting (Runs only if sequences were selected)
      if (selected_seqs2_defined & ! multiple_movements) {
        
        # Filtering data
        seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs2)
        genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs2)
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
          geom_gene(aes(y = y)) +
          geom_link(aes(y = y, fill = pident), colour = NA) +
          geom_bin_label(aes(y = y), x = -10) +
          scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
        
        # Saving the plot
        ggsave(p_genome, filename = paste(Cluster_number,"CargoSynteny_SubCluster", ClusterId, ".svg", sep = ""),
               width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)),
               height = min(49, length(selected_seqs2)), limitsize = FALSE)
      } else if(selected_seqs2_defined & multiple_movements) {

        k_out_cluster <- unique(test_NbClust$Best.partition[elements_other])
        k_cluster <- unique(test_NbClust$Best.partition[elements_cluster])

        for(i in k_cluster) {
          for(j in k_out_cluster) {
              selected_seqs2 <- c(names(test_NbClust$Best.partition[grep(i,test_NbClust$Best.partition)][elements_cluster]),names(test_NbClust$Best.partition[grep(j,test_NbClust$Best.partition)][elements_other]))
              #selected_seqs2 <- selected_seqs2[!is.na(selected_seqs2)]
              selected_seqs2 <- selected_seqs2[!is.na(selected_seqs2)]

              # Filtering data
              seqs_filtered <- seq_data %>% filter(seq_id %in% selected_seqs2)
              genes_filtered <- gene_data %>% filter(seq_id %in% selected_seqs2)
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
              geom_gene(aes(y = y)) +
              geom_link(aes(y = y, fill = pident), colour = NA) +
              geom_bin_label(aes(y = y), x = -10) +
              scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
              theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm"))
        
              # Saving the plot
              ggsave(p_genome, filename = paste(Cluster_number,"CargoSynteny_SubCluster", ClusterId, "-", i,"vs",j,".svg", sep = ""),
                 width = min(49, max(16, round(max(ordered_seqs$length) * 0.0001) / 2)),
                 height = min(49, length(selected_seqs2)), limitsize = FALSE)
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
# Call directories in the Working directory
setwd(arguments$directory)

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Reading data...","\n", sep=""))

data_list <- load_and_preprocess_data()

# Preprocess data
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Identifying if there are separate clusters","\n", sep=""))

clustering_data <- perform_initial_clustering()

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

if((clustering_data$Individual_clusters == 1) & (arguments$subclusters >= 2) & (arguments$captainRemoval <= 1)){
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))
  
  core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl)

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing ",arguments$subclusters," subclusters","\n", sep=""))

  fit <- cutree(clustering_data$Hierar_cl, k = arguments$subclusters) # k is the number of subclusters from the first approach
  New_dist <- as.matrix(clustering_data$distance_mat)

  plot_subcluster_synteny(
    subcluster_number = arguments$subclusters,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts, 
    core_genes = core_genes,   
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results)

} else if((clustering_data$Individual_clusters == 1) & (arguments$subclusters >= 2) & (arguments$captainRemoval >= 2)) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))
  
  core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl)

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))

  test_NbClust <- NbClust(clustering_data$distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(clustering_data$distance_mat)-1), index = "ball")
  test_fit <- cutree(clustering_data$Hierar_cl, h = 0.6)
  K_number <- max(length(unique(test_NbClust$Best.partition)), length(unique(fit)))
  fit <- cutree(clustering_data$Hierar_cl, k = K_number)

  New_dist <- as.matrix(clustering_data$distance_mat)

  plot_subcluster_synteny(
    subcluster_number = K_number,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts, 
    core_genes = core_genes,   
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results)

} else if((clustering_data$Individual_clusters > 1) & (arguments$subclusters >= 2) & length(Selected_Cluster_ID) == 1) {

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))

  New_ortho_counts <- clustering_data$orthofinder_counts[,grep(Selected_Cluster_ID,clustering_data$cluster_fit)]
  Cluster_matrix <- clustering_data$transposed_counts[grep(Selected_Cluster_ID,clustering_data$cluster_fit),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  distance_mat <- dist(Cluster_matrix, method='binary')
  Hierar_cl <- hclust(distance_mat, method = 'average')
  
  core_genes <- core_genes_analysis(ortho_counts = New_ortho_counts, Cluster = Hierar_cl,
    Cluster_number = paste("_Cluster", Selected_Cluster_ID, sep=""))

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))
  
  Cluster_matrix <- clustering_data$transposed_counts[grep(Selected_Cluster_ID,clustering_data$cluster_fit),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  distance_mat <- dist(Cluster_matrix, method='binary')
  Hierar_cl <- hclust(distance_mat, method = 'average')

  test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(distance_mat)-1), index = "ball")
  test_fit <- cutree(Hierar_cl, h = 0.6)
  K_number <- max(length(unique(test_NbClust$Best.partition)), length(unique(test_fit)))
  fit <- cutree(Hierar_cl, k = K_number)

  New_dist <- as.matrix(distance_mat)

  plot_subcluster_synteny(
    subcluster_number = K_number,
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

    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes in",MainClusterID,"Cluster","\n", sep=""))

    New_ortho_counts <- clustering_data$orthofinder_counts[,grep(MainClusterID,clustering_data$cluster_fit)]
    distance_mat <- dist(Cluster_matrix, method='binary')
    Hierar_cl <- hclust(distance_mat, method = 'average')
  
    core_genes <- core_genes_analysis(ortho_counts = New_ortho_counts, Cluster = Hierar_cl,
      Cluster_number = paste("_Cluster", MainClusterID, sep=""))

    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters in",MainClusterID,"Cluster","\n", sep=""))

    

    test_NbClust <- NbClust(distance_mat, method = "average", min.nc = 1, max.nc = min(6, nrow(Cluster_matrix)-1), index = "ball")
    test_fit <- cutree(Hierar_cl, h = 0.6)
    K_number <- max(length(unique(test_NbClust$Best.partition)), length(unique(test_fit)))
    fit <- cutree(Hierar_cl, k = K_number)

    New_dist <- as.matrix(distance_mat)

    plot_subcluster_synteny(
    subcluster_number = K_number,
    dist_matrix = New_dist,    
    cluster_fit = fit,    
    ortho_counts = clustering_data$orthofinder_counts, 
    core_genes = core_genes,   
    seq_data = data_list$seqs,        
    gene_data = data_list$genes,      
    blast_links = data_list$blast_results,
    Cluster_number = paste("Cluster", MainClusterID, sep=""))

  } 
} else{
  cat(paste("    [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No subcluster to analyzed","\n", sep=""))
} 

if((clustering_data$Individual_clusters == 1) & (arguments$subclusters == 1)) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are possible core genes","\n", sep=""))
  
  core_genes <- core_genes_analysis(ortho_counts = clustering_data$orthofinder_counts, Cluster = clustering_data$Hierar_cl)

  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are differences between captain and cargo dendograms","\n", sep=""))

  phylo_tree <- read.tree("CaptainPhylogeny.nw")
  if(length(phylo_tree$tip.label) > 4){
    dist_matrix <- cophenetic.phylo(phylo_tree)
    min_dist <- min(dist_matrix)
    max_dist <- max(dist_matrix)
    normalized_dist_matrix <- (dist_matrix - min_dist) / (max_dist - min_dist)
    normalized_dist_matrix2 <- scale(normalized_dist_matrix)
    hclust_tree <- hclust(as.dist(normalized_dist_matrix), method = "average")

    test_NbClust <- NbClust(normalized_dist_matrix2, method = "average", min.nc = 1, max.nc = min(6, nrow(normalized_dist_matrix2)-1), index = "ball")
    test_fit <- cutree(hclust_tree, h = 0.5)
    K_number <- max(length(unique(test_NbClust$Best.partition)),length(unique(test_fit)))
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

      dend_list <- dendlist(clustering_data$Hierar_cl, hclust_tree)
      svg(filename="CaptainVsCargo_tanglegram.svg", width = 16, height = min(49, length(fit)/2))
      dend_list %>% ladderize %>% 
      untangle(method = "step1side", k_seq = K_number:(length(fit)-1)) %>%
      set("branches_k_color", k=K_number) %>% 
      tanglegram(faster = TRUE, main_left = as.character("Cargo Cluster"),main_right = as.character("Captain tree"))
      invisible(dev.off())
    } else {
      cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No discordance detected","\n", sep=""))
    }
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

      dend_list <- dendlist(clustering_data$Hierar_cl, hclust_tree)
      svg(filename="CaptainVsCargo_tanglegram.svg", width = 16, height = min(49, length(fit)/2))
      dend_list %>% ladderize %>% 
      untangle(method = "step1side", k_seq = 2:(length(fit)-1)) %>%
      set("branches_k_color", k=2) %>% 
      tanglegram(faster = TRUE, main_left = as.character("Cargo Cluster"),main_right = as.character("Captain tree"))
      invisible(dev.off())
    } else {
      cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No discordance detected","\n", sep=""))
    }
  } 
} 

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Finished","\n", sep=""))
