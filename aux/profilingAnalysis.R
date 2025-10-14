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
  make_option(c("-c", "--captainRemoval"), action = "store_true", default=FALSE,
              help=" Flag that identify if a captain have been remove from the original dataset. Store true if selected.")
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
suppressPackageStartupMessages(library(factoextra))
suppressPackageStartupMessages(library(dendextend))

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
if (!requireNamespace("factoextra", quietly = TRUE)) {
   stop("Package \"factoextra\" not installed. Please install it to run this script.", call. = FALSE)
}
if (!requireNamespace("dendextend", quietly = TRUE)) {
   stop("Package \"dendextend\" not installed. Please install it to run this script.", call. = FALSE)
}

  
# Call directories in the Working directory
setwd(arguments$directory)

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Reading data...","\n", sep=""))

col_names <- c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "pident", "length", "qlen", "slen")
col_classes <- c("character", "character", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric")

blast_results <- read.delim(
  "clean_results.txt",
  header = FALSE,
  sep = "\t",
  col.names = col_names,
  colClasses = col_classes
) %>%
  filter(qseqid != sseqid, pident > 70)

links <- blast_results %>%
  mutate(alignment_id = row_number()) %>%
  select(
    seq_id = qseqid,
    start = qstart,
    end = qend,
    seq_id2 = sseqid,
    start2 = sstart,
    end2 = send,
    pident,
    length
  )

seqs <- bind_rows( blast_results %>% select(seq_id = qseqid, length = qlen), blast_results %>% select(seq_id = sseqid, length = slen)) %>% distinct() %>% group_by(seq_id) %>% summarize(length = max(length)) %>% ungroup()

gff_col_names <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
gff_col_classes <- c("character", "character", "character", "numeric", "numeric", "character", "character", "character", "character")

genes <- read.delim("Final_model.gff",
                    header = FALSE,
                    comment.char = "#",
                    sep = "\t",
                    col.names = gff_col_names,
                    colClasses = gff_col_classes) %>%
  filter(feature == "gene") %>%
  select(seq_id = seqname, start, end, strand) %>%
  mutate(seq_id = as.character(seq_id),
         type = "CDS") 


# Preprocess data
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Identifying if there are separate clusters","\n", sep=""))
OrthoFinder <- read.table("Orthogroups.GeneCount.tsv", header = T, check.names = F, row.names = 1)
OrthoFinder2 <- OrthoFinder[ , ! names(OrthoFinder) %in% c("Total")]
transposase_OrthoFinder <- t(OrthoFinder2)
distance_mat <- dist(transposase_OrthoFinder, method='binary')
Hierar_cl <- hclust(distance_mat, method = 'average')

#Identify number of clusters
fit_Orthofinder <- cutree(Hierar_cl, h = 0.99999999)
Individual_clusters <- length(unique(fit_Orthofinder))

# Process Clusters
cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing ",Individual_clusters," individual cluster(s) identified","\n", sep=""))

# Process each cluster
Selected_clusters <- numeric()
for( ClusterId in 1:Individual_clusters ) {
  Cluster_matrix <- transposase_OrthoFinder[grep(ClusterId,fit_Orthofinder),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  if(ncol(Cluster_matrix) > 2){
    Cluster_distance <- dist(Cluster_matrix, method='binary')
    Cluster_Hier <- hclust(Cluster_distance, method = 'average')
    my_tree <- as.phylo(Cluster_Hier)

    if(Individual_clusters >= 2){
      write.tree(phy=my_tree, file=paste("CargoHC_Cluster",ClusterId,".nwk",sep=""))
      write.table(Cluster_Hier$labels, file =paste("Cluster",ClusterId,".txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
    } else {
      write.tree(phy=my_tree, file="CargoHC.nwk")
    }

    Cluster_melt <- melt(Cluster_matrix, as.is = T)
    Profile_plot <- ggplot(data = Cluster_melt, aes(x=reorder(Var2,value), y=Var1, fill=value)) + 
    geom_tile() + scale_fill_viridis(option="rocket", direction = -1) + labs(y="element",x="Orthogroup") + theme(axis.text.x=element_blank())

    if(Individual_clusters >= 2){
      ggsave(Profile_plot,filename=paste("CargoProfiling_Cluster",ClusterId,".svg",sep=""), width = min(49, ncol(Cluster_matrix)*0.2+2), height = min(49, nrow(Cluster_matrix)))
    } else {
      ggsave(Profile_plot,filename="CargoProfiling.svg", width = min(49, ncol(Cluster_matrix)*0.2+2), height = min(49, nrow(Cluster_matrix)))
    }
  
    tree_sorted <- ladderize(my_tree, right = FALSE)
    num_sequences <- length(tree_sorted$tip.label)

    selected_seqs <- row.names(Cluster_matrix)

    seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs)

    genes_filtered <- genes %>% filter(seq_id %in% selected_seqs)
    
    links_filtered <- blast_results %>%
      filter(qseqid %in% selected_seqs | sseqid %in% selected_seqs) %>%
      select(
      seq_id = qseqid,
      start = qstart,
      end = qend,
      seq_id2 = sseqid,
      start2 = sstart,
      end2 = send,
      pident)

    p_tree_base <- ggtree(tree_sorted, hang = 0.05) +
    theme(plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm")) # c(top, right, bottom, left)
    tree_data <- p_tree_base$data

    tree_y_coords <- tree_data %>%
      filter(isTip) %>%
      select(seq_id = label, y) %>%
      mutate(y = max(y) - y + 1)

    ordered_seqs <- seqs_filtered %>%
      inner_join(tree_y_coords, by = "seq_id") %>%
      arrange(y)

    p_tree <- p_tree_base +
      geom_tiplab(align = TRUE, size = 3) +
      theme_tree2() 

    p_genome <- gggenomes(
      seqs = ordered_seqs,
      links = links_filtered,
      genes = genes_filtered
    ) +
    geom_seq(aes(y = y)) +
    geom_gene(aes(y = y)) +
    geom_link(aes(y = y, fill = pident), colour = NA ) +
    geom_bin_label(aes(y = y), x = -10) +
    scale_x_continuous(
    labels = label_number(accuracy = 1),
    limits = c(0, max(ordered_seqs$length))
    ) +
    #scale_fill_gradient(low = "gray80", high = "gray60") +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) # c(top, right, bottom, left)

    final_plot <- p_tree + p_genome

    if(Individual_clusters >= 2){
      ggsave(final_plot,filename=paste("CargoSynteny_Cluster",ClusterId,".svg",sep=""), width = min(49, max(16,round(max(ordered_seqs$length)*0.0001)+round(max(tree_sorted$edge.length)*10))), height = min(49, nrow(Cluster_matrix)),limitsize = FALSE)
      Selected_clusters <- c(Selected_clusters,ClusterId)
    } else {
      ggsave(final_plot,filename="CargoSynteny.svg", width = min(49, max(16,round(max(ordered_seqs$length)*0.0001)+round(max(tree_sorted$edge.length)*10))), height = min(49, nrow(Cluster_matrix)),limitsize = FALSE)
    }
  } else {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","WARNING: There are issues with the following element: ",as.character(rownames(Cluster_matrix)),"\n", sep=""))
  }
  
}

if((Individual_clusters == 1) & (arguments$subclusters >= 2) & ! arguments$captainRemoval){
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing ",arguments$subclusters," subclusters","\n", sep=""))

  fit <- cutree(Hierar_cl, k = arguments$subclusters) # k is the number of subclusters from the first approach
  Table_fit <- as.data.frame(table(fit))
  New_dist <- as.matrix(distance_mat)

  for( ClusterId in 1:arguments$subclusters ) {

    matrix <- New_dist[grep(ClusterId,fit),] # The number is the number of the cluster that want to analyze
    Cluster_elements <- rownames(matrix)

    if(is.null(Cluster_elements)) {
      Cluster_elements <- names(grep(ClusterId,fit, value = TRUE))
    }

    write.table(Cluster_elements, file =paste("SubCluster",ClusterId,".txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)

    if(is.matrix(matrix)) {
      reduce_matrix <- matrix[,colSums(matrix) < nrow(matrix)*0.98]
      reduce_matrix <- reduce_matrix[,names(sort(colSums(reduce_matrix), decreasing = F))]
      if(ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2)*0.98,]
        selected_seqs2 <- c(rownames(reduce_matrix3),colnames(reduce_matrix3))

        OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
        OrthoFinder_subcluster2 <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster) > 0,]
        Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 > 0) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

        if(length(Gene_movement) > 0) {
          write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
        }

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>% filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
        filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
        select(
        seq_id = qseqid,
        start = qstart,
        end = qend,
        seq_id2 = sseqid,
        start2 = sstart,
        end2 = send,
        pident)

        ordered_seqs <- seqs_filtered %>%
        arrange(match(seq_id, selected_seqs2))
      } else if(ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        selected_seqs2 <- c(names(reduce_matrix2[reduce_matrix2 < 1]),tail(colnames(reduce_matrix), n = 1))

        OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
        OrthoFinder_subcluster2 <- OrthoFinder_subcluster[OrthoFinder_subcluster[, ncol(OrthoFinder_subcluster)] > 0,]
        Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 >=1) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

        if(length(Gene_movement) > 0) {
          write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
        }

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>% filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
        filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
        select(
        seq_id = qseqid,
        start = qstart,
        end = qend,
        seq_id2 = sseqid,
        start2 = sstart,
        end2 = send,
        pident)

        ordered_seqs <- seqs_filtered %>%
        arrange(match(seq_id, selected_seqs2))
      }
      if(exists("selected_seqs2")) {
        p_genome <- gggenomes(
        seqs = ordered_seqs,
        links = links_filtered,
        genes = genes_filtered
        ) +
        geom_seq(aes(y = y)) +
        geom_gene(aes(y = y)) +
        geom_link(aes(y = y, fill = pident), colour = NA ) +
        geom_bin_label(aes(y = y), x = -10) +
        scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
        #scale_fill_gradient(low = "gray80", high = "gray60") +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) # c(top, right, bottom, left)

        ggsave(p_genome,filename=paste("CargoSynteny_SubCluster",ClusterId,".svg",sep=""), width = min(49,max(16,round(max(ordered_seqs$length)*0.0001)/2)), height = min(49, length(selected_seqs2)), limitsize = FALSE)
      } else {
        cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Subcluster ",ClusterId," can not be analyzed automatically","\n", sep=""))
      }
    } else {
        cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Subcluster ",ClusterId," can not be analyzed automatically","\n", sep=""))
      }
  }
} else if((Individual_clusters == 1) & (arguments$subclusters >= 2) & arguments$captainRemoval) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))

  New_dist <- as.matrix(distance_mat)

  fviz_output <- fviz_nbclust(New_dist, FUN = hcut, method = "wss", k.max = min(10, nrow(New_dist)-1))
  wss_data <- fviz_output$data
  slopes <- diff(wss_data$y)
  slope_change <- abs(diff(slopes))
  elbow_point <- which.max(slope_change) + 1

  fit <- cutree(Hierar_cl, k = elbow_point) # k is the number of subclusters from the elbow approach

  Table_fit <- as.data.frame(table(fit))

  for( ClusterId in 1:elbow_point) {

    matrix <- New_dist[grep(ClusterId,fit),] # The number is the number of the cluster that want to analyze
    Cluster_elements <- rownames(matrix)

    if(is.null(Cluster_elements)) {
      Cluster_elements <- names(grep(ClusterId,fit, value = TRUE))
    }

    write.table(Cluster_elements, file =paste("SubCluster",ClusterId,".txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)

    if(is.matrix(matrix)) {
      reduce_matrix <- matrix[,colSums(matrix) < nrow(matrix)*0.98]
      reduce_matrix <- reduce_matrix[,names(sort(colSums(reduce_matrix), decreasing = F))]
      if(ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2)*0.98,]
        selected_seqs2 <- c(rownames(reduce_matrix3),colnames(reduce_matrix3))

        OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
        OrthoFinder_subcluster2 <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster) > 0,]
        Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 > 0) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

        if(length(Gene_movement) > 0) {
          write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
        }

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>%
        filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
        filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
        select(
        seq_id = qseqid,
        start = qstart,
        end = qend,
        seq_id2 = sseqid,
        start2 = sstart,
        end2 = send,
        pident)

        ordered_seqs <- seqs_filtered %>%
        arrange(match(seq_id, selected_seqs2))
      } else if(ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        selected_seqs2 <- c(names(reduce_matrix2[reduce_matrix2 < 1]),tail(colnames(reduce_matrix), n = 1))

        OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
        OrthoFinder_subcluster2 <- OrthoFinder_subcluster[OrthoFinder_subcluster[, ncol(OrthoFinder_subcluster)] > 0,]
        Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 >=1) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

        if(length(Gene_movement) > 0) {
          write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
        }

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>%
        filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
        filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
        select(
        seq_id = qseqid,
        start = qstart,
        end = qend,
        seq_id2 = sseqid,
        start2 = sstart,
        end2 = send,
        pident)

        ordered_seqs <- seqs_filtered %>%
        arrange(match(seq_id, selected_seqs2))
      }

      if(exists("selected_seqs2")) {
        p_genome <- gggenomes(
        seqs = ordered_seqs,
        links = links_filtered,
        genes = genes_filtered
        ) +
        geom_seq(aes(y = y)) +
        geom_gene(aes(y = y)) +
        geom_link(aes(y = y, fill = pident), colour = NA ) +
        geom_bin_label(aes(y = y), x = -10) +
        scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
        #scale_fill_gradient(low = "gray80", high = "gray60") +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) # c(top, right, bottom, left)

        ggsave(p_genome,filename=paste("CargoSynteny_SubCluster",ClusterId,".svg",sep=""), width = min(49,max(16,round(max(ordered_seqs$length)*0.0001)/2)), height = min(49, length(selected_seqs2)), limitsize = FALSE)
      } else{
        cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Subcluster ",ClusterId," can not be analyzed automatically","\n", sep=""))
      }
    }
  }
} else if((Individual_clusters > 1) & (arguments$subclusters >= 2) & length(Selected_clusters) == 1) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))

  Cluster_matrix <- transposase_OrthoFinder[grep(Selected_clusters,fit_Orthofinder),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  distance_mat <- dist(Cluster_matrix, method='binary')
  Hierar_cl <- hclust(distance_mat, method = 'average')

  New_dist <- as.matrix(distance_mat)

  fviz_output <- fviz_nbclust(New_dist, FUN = hcut, method = "wss", k.max = min(10, nrow(New_dist)-1))
  wss_data <- fviz_output$data
  slopes <- diff(wss_data$y)
  slope_change <- abs(diff(slopes))
  elbow_point <- which.max(slope_change) + 1

  fit <- cutree(Hierar_cl, k = elbow_point) # k is the number of subclusters from the elbow approach

  Table_fit <- as.data.frame(table(fit))

  for( ClusterId in 1:elbow_point) {

    matrix <- New_dist[grep(ClusterId,fit),] # The number is the number of the cluster that want to analyze
    Cluster_elements <- rownames(matrix)

    if(is.null(Cluster_elements)) {
      Cluster_elements <- names(grep(ClusterId,fit, value = TRUE))
    }

    write.table(Cluster_elements, file =paste("SubCluster",ClusterId,"_",Selected_clusters,".txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)

    if(is.matrix(matrix)) {
      reduce_matrix <- matrix[,colSums(matrix) < nrow(matrix)*0.98]
      reduce_matrix <- reduce_matrix[,names(sort(colSums(reduce_matrix), decreasing = F))]
      if(ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2)*0.98,]
        selected_seqs2 <- c(rownames(reduce_matrix3),colnames(reduce_matrix3))

        OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
        OrthoFinder_subcluster2 <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster) > 0,]
        Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 > 0) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

        if(length(Gene_movement) > 0) {
          write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_",Selected_clusters,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
        }

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>%
          filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
          filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(
          seq_id = qseqid,
          start = qstart,
          end = qend,
          seq_id2 = sseqid,
          start2 = sstart,
          end2 = send,
          pident)

        ordered_seqs <- seqs_filtered %>%
          arrange(match(seq_id, selected_seqs2))
      } else if(ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        selected_seqs2 <- c(names(reduce_matrix2[reduce_matrix2 < 1]),tail(colnames(reduce_matrix), n = 1))

        OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
        OrthoFinder_subcluster2 <- OrthoFinder_subcluster[OrthoFinder_subcluster[, ncol(OrthoFinder_subcluster)] > 0,]
        Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 >=1) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

        if(length(Gene_movement) > 0) {
          write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_",Selected_clusters,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
        }

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>%
          filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
          filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(
          seq_id = qseqid,
          start = qstart,
          end = qend,
          seq_id2 = sseqid,
          start2 = sstart,
          end2 = send,
          pident)

        ordered_seqs <- seqs_filtered %>%
        arrange(match(seq_id, selected_seqs2))
      }
      if(exists("selected_seqs2")) {
        p_genome <- gggenomes(
        seqs = ordered_seqs,
        links = links_filtered,
        genes = genes_filtered
        ) +
        geom_seq(aes(y = y)) +
        geom_gene(aes(y = y)) +
        geom_link(aes(y = y, fill = pident), colour = NA ) +
        geom_bin_label(aes(y = y), x = -10) +
        scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
        #scale_fill_gradient(low = "gray80", high = "gray60") +
        theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) # c(top, right, bottom, left)

        ggsave(p_genome,filename=paste("CargoSynteny_SubCluster",ClusterId,"_Cluster",Selected_clusters,".svg",sep=""), width = min(49,max(16,round(max(ordered_seqs$length)*0.0001)/2)), height = min(49, length(selected_seqs2)), limitsize = FALSE)
      } else{
        cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Subcluster ",ClusterId," from main cluster ",Selected_clusters," can not be analyzed automatically","\n", sep=""))
      }
    }
  }
} else if((Individual_clusters > 1) & (arguments$subclusters >= 2) & length(Selected_clusters) > 1) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing subclusters","\n", sep=""))

  for(MainClusterID in Selected_clusters ) {
    Cluster_matrix <- transposase_OrthoFinder[grep(MainClusterID,fit_Orthofinder),]
    Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]

    if(ncol(Cluster_matrix) < 4) {
      next
    }

    distance_mat <- dist(Cluster_matrix, method='binary')
    Hierar_cl <- hclust(distance_mat, method = 'average')
    New_dist <- as.matrix(distance_mat)

    fviz_output <- fviz_nbclust(New_dist, FUN = hcut, method = "wss", k.max = min(10, nrow(New_dist)-1))
    wss_data <- fviz_output$data
    slopes <- diff(wss_data$y)
    slope_change <- abs(diff(slopes))
    elbow_point <- which.max(slope_change) + 1

    fit <- cutree(Hierar_cl, k = elbow_point) # k is the number of subclusters from the elbow approach

    Table_fit <- as.data.frame(table(fit))

    for( ClusterId in 1:elbow_point) {

      matrix <- New_dist[grep(ClusterId,fit),] # The number is the number of the cluster that want to analyze
      Cluster_elements <- rownames(matrix)

      if(is.null(Cluster_elements)) {
        Cluster_elements <- names(grep(ClusterId,fit, value = TRUE))
      }

      write.table(Cluster_elements, file =paste("SubCluster",ClusterId,"_",MainClusterID,".txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)

      if(is.matrix(matrix)) {
        reduce_matrix <- matrix[,colSums(matrix) < nrow(matrix)*0.98]
        reduce_matrix <- reduce_matrix[,names(sort(colSums(reduce_matrix), decreasing = F))]
        if(ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
          reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
          reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2)*0.98,]
          selected_seqs2 <- c(rownames(reduce_matrix3),colnames(reduce_matrix3))

          OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
          OrthoFinder_subcluster2 <- OrthoFinder_subcluster[rowSums(OrthoFinder_subcluster) > 0,]
          Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 > 0) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

          if(length(Gene_movement) > 0) {
            write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_",MainClusterID,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
          }

          seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
          genes_filtered <- genes %>%
          filter(seq_id %in% selected_seqs2)
    
          links_filtered <- blast_results %>%
          filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(
          seq_id = qseqid,
          start = qstart,
          end = qend,
          seq_id2 = sseqid,
          start2 = sstart,
          end2 = send,
          pident)

          ordered_seqs <- seqs_filtered %>%
          arrange(match(seq_id, selected_seqs2))
        } else if(ncol(reduce_matrix) > nrow(reduce_matrix)) {
          reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
          selected_seqs2 <- c(names(reduce_matrix2[reduce_matrix2 < 1]),tail(colnames(reduce_matrix), n = 1))

          OrthoFinder_subcluster <- OrthoFinder2[,selected_seqs2]
          OrthoFinder_subcluster2 <- OrthoFinder_subcluster[OrthoFinder_subcluster[, ncol(OrthoFinder_subcluster)] > 0,]
          Gene_movement <- rownames(OrthoFinder_subcluster2[rowSums(OrthoFinder_subcluster2 >=1) >= round(ncol(OrthoFinder_subcluster2)*0.6),])

          if(length(Gene_movement) > 0) {
            write.table(Gene_movement, file =paste("SubCluster",ClusterId,"_",MainClusterID,"_moveOrthologs.txt",sep=""), sep = '\t', row.names = F, col.names= F,quote = F)
          }

          seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
          genes_filtered <- genes %>%
          filter(seq_id %in% selected_seqs2)
    
          links_filtered <- blast_results %>%
          filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(
          seq_id = qseqid,
          start = qstart,
          end = qend,
          seq_id2 = sseqid,
          start2 = sstart,
          end2 = send,
          pident)

          ordered_seqs <- seqs_filtered %>%
          arrange(match(seq_id, selected_seqs2))
        }
        if(exists("selected_seqs2")) {
          p_genome <- gggenomes(
          seqs = ordered_seqs,
          links = links_filtered,
          genes = genes_filtered
          ) +
          geom_seq(aes(y = y)) +
          geom_gene(aes(y = y)) +
          geom_link(aes(y = y, fill = pident), colour = NA ) +
          geom_bin_label(aes(y = y), x = -10) +
          scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
          #scale_fill_gradient(low = "gray80", high = "gray60") +
          theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) # c(top, right, bottom, left)

          ggsave(p_genome,filename=paste("CargoSynteny_SubCluster",ClusterId,"_",MainClusterID,".svg",sep=""), width = min(49,max(16,round(max(ordered_seqs$length)*0.0001)/2)), height = min(49, length(selected_seqs2)), limitsize = FALSE)
        } else{
          cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Subcluster ",ClusterId," from main cluster",MainClusterID," can not be analyzed automatically","\n", sep=""))
        }
      }
    }
  }
} else{
  cat(paste("    [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No subcluster to analyzed","\n", sep=""))
} 

if((Individual_clusters == 1) & (arguments$subclusters == 1)) {
  cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Analyzing if there are differences between captain and cargo dendograms","\n", sep=""))

  phylo_tree <- read.tree("CaptainPhylogeny.nw")
  dist_matrix <- cophenetic.phylo(phylo_tree)
  min_dist <- min(dist_matrix)
  max_dist <- max(dist_matrix)
  normalized_dist_matrix <- (dist_matrix - min_dist) / (max_dist - min_dist)
  normalized_dist_matrix2 <- scale(normalized_dist_matrix)
  fviz_output <- fviz_nbclust(normalized_dist_matrix2, FUN = hcut, method = "wss", k.max = min(10, nrow(normalized_dist_matrix2)-1))
  wss_data <- fviz_output$data
  slopes <- diff(wss_data$y)
  slope_change <- abs(diff(slopes))
  elbow_point <- which.max(slope_change) + 1
  hclust_tree <- hclust(as.dist(normalized_dist_matrix), method = "average")
  fit_tree <- cutree(hclust_tree, k = elbow_point)
  fit_tree <- fit_tree[order(names(fit_tree))]

  fit <- cutree(Hierar_cl, k = elbow_point) # k is the number of subclusters from the elbow method
  fit <- fit[order(names(fit))]

  contingency_table <- table(fit_tree, fit)
  best_match <- apply(contingency_table, 1, which.max)
  fit_aligned <- fit
  for(i in 1:length(best_match)){
    original_label <- names(best_match)[i]
    new_label <- best_match[i]
    fit_aligned[fit == original_label] <- as.numeric(new_label)
  }

  discordant_indices <- which(fit_tree != fit_aligned)
  discordant_elements <- names(fit_tree)[discordant_indices]

  if(length(discordant_elements) >= 1) {
    cat(paste("    [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","There ",length(discordant_elements)," discordant element(s) identified","\n", sep=""))
    write.table(discordant_elements, file = "Discordant_elements.txt", sep = '\t', row.names = F, col.names= F,quote = F)

    dend_list <- dendlist(Hierar_cl, hclust_tree)
    svg(filename="CaptainVsCargo_tanglegram.svg", width = 16, height = min(49, length(fit)/2))
    dend_list %>% ladderize %>% 
      untangle(method = "step1side", k_seq = elbow_point:(length(fit)-1)) %>%
      set("branches_k_color", k=elbow_point) %>% 
      tanglegram(faster = TRUE, main_left = as.character("Cargo Cluster"),main_right = as.character("Captain tree"))
    dev.off() 

    New_dist <- as.matrix(distance_mat)
    matrix <- New_dist[discordant_elements,]

    if(is.matrix(matrix)) {
      reduce_matrix <- matrix[,colSums(matrix) < nrow(matrix)*0.5]
      reduce_matrix <- reduce_matrix[,names(sort(colSums(reduce_matrix), decreasing = F))]
      if(ncol(reduce_matrix) >= nrow(reduce_matrix) + 2) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        reduce_matrix3 <- reduce_matrix2[rowSums(reduce_matrix2) < ncol(reduce_matrix2)*0.98,]
        selected_seqs2 <- c(rownames(reduce_matrix3),colnames(reduce_matrix3))

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>%
        filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
        filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
        select(
        seq_id = qseqid,
        start = qstart,
        end = qend,
        seq_id2 = sseqid,
        start2 = sstart,
        end2 = send,
        pident)

      } else if(ncol(reduce_matrix) > nrow(reduce_matrix)) {
        reduce_matrix2 <- reduce_matrix[,(nrow(reduce_matrix)+1):ncol(reduce_matrix)]
        selected_seqs2 <- c(names(reduce_matrix2[reduce_matrix2 < 1]),tail(colnames(reduce_matrix), n = 1))

        seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
        genes_filtered <- genes %>%
        filter(seq_id %in% selected_seqs2)
    
        links_filtered <- blast_results %>%
        filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
        select(
        seq_id = qseqid,
        start = qstart,
        end = qend,
        seq_id2 = sseqid,
        start2 = sstart,
        end2 = send,
        pident)

      }
      } else {
        reduce_matrix <- matrix[matrix < 0.5]

        if(length(reduce_matrix) > 1) {
          selected_seqs2 <- unique(c(discordant_elements, names(sort(matrix[matrix<0.5]))))
          seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs2)
    
          genes_filtered <- genes %>%
          filter(seq_id %in% selected_seqs2)
    
          links_filtered <- blast_results %>%
          filter(qseqid %in% selected_seqs2 | sseqid %in% selected_seqs2) %>%
          select(
          seq_id = qseqid,
          start = qstart,
          end = qend,
          seq_id2 = sseqid,
          start2 = sstart,
          end2 = send,
          pident)
        }
      }

    if(exists("selected_seqs2")) {
      p_genome <- gggenomes(
      seqs = seqs_filtered,
      links = links_filtered,
      genes = genes_filtered
      ) +
      geom_seq(aes(y = y)) +
      geom_gene(aes(y = y)) +
      geom_link(aes(y = y, fill = pident), colour = NA ) +
      geom_bin_label(aes(y = y), x = -10) +
      scale_x_continuous(labels = label_number(accuracy = 1), limits = c(0, max(ordered_seqs$length))) +
      #scale_fill_gradient(low = "gray80", high = "gray60") +
      theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) # c(top, right, bottom, left)

      ggsave(p_genome,filename="Captain_discordance.svg", width = min(49,max(16,round(max(ordered_seqs$length)*0.0001)/2)), height = min(49, length(selected_seqs2)), limitsize = FALSE)
    } else{
      cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","discordances can not be analyzed automatically","\n", sep=""))
    }
  } else {
    cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","No discordance detected","\n", sep=""))
  }
} 

cat(paste("  [",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"] ","Finished","\n", sep=""))
