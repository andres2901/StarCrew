# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory [default %default]",metavar="PATH")
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

  
# Call directories in the Working directory
setwd(arguments$directory)

cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]","Reading data...","\n"))

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
cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]","Identifying if there are separate clusters","\n"))
OrthoFinder <- read.table("Orthogroups.GeneCount.tsv", header = T, check.names = F, row.names = 1)
OrthoFinder2 <- OrthoFinder[ , ! names(OrthoFinder) %in% c("Total")]
transposase_OrthoFinder <- t(OrthoFinder2)
distance_mat <- dist(transposase_OrthoFinder, method='binary')
Hierar_cl <- hclust(distance_mat, method = 'average')

#Identify number of clusters
fit_Orthofinder <- cutree(Hierar_cl, h = 0.99999999)
Individual_clusters <- length(unique(fit_Orthofinder))

# Process Clusters
cat(paste("[",format(Sys.time(), "%Y-%m-%d %H:%M:%S"),"]","Analyzing",Individual_clusters,"individual cluster(s) identified","\n"))

# Process each cluster
for( ClusterId in 1:Individual_clusters ) {
  Cluster_matrix <- transposase_OrthoFinder[grep(ClusterId,fit_Orthofinder),]
  Cluster_matrix <- Cluster_matrix[,colSums(abs(Cluster_matrix)) > 0]
  Cluster_distance <- dist(Cluster_matrix, method='binary')
  Cluster_Hier <- hclust(Cluster_distance, method = 'average')
  my_tree <- as.phylo(Cluster_Hier)

  if(Individual_clusters >= 2){
    write.tree(phy=my_tree, file=paste("CargoHC_Cluster",ClusterId,".nwk",sep=""))
  } else {
    write.tree(phy=my_tree, file="CargoHC.nwk")
  }

  Cluster_melt <- melt(Cluster_matrix, as.is = T)
  Profile_plot <- ggplot(data = Cluster_melt, aes(x=reorder(Var2,value), y=Var1, fill=value)) + 
  geom_tile() + scale_fill_viridis(option="rocket", direction = -1) + labs(y="element",x="Orthogroup") + theme(axis.text.x=element_blank())



  if(Individual_clusters >= 2){
    ggsave(Profile_plot,filename=paste("CargoProfiling_Cluster",ClusterId,".svg",sep=""), width = ncol(Cluster_matrix)*0.2+2, height = nrow(Cluster_matrix))
  } else {
    ggsave(Profile_plot,filename="CargoProfiling.svg", width = ncol(Cluster_matrix)*0.2+2, height = nrow(Cluster_matrix))
  }
  

  tree_sorted <- ladderize(my_tree, right = FALSE)
  num_sequences <- length(tree_sorted$tip.label)

  selected_seqs <- row.names(Cluster_matrix)

  seqs_filtered <- seqs %>% filter(seq_id %in% selected_seqs)
  genes_filtered <- genes %>%
    filter(seq_id %in% selected_seqs)
  links_filtered <- blast_results %>%
    filter(qseqid %in% selected_seqs | sseqid %in% selected_seqs) %>%
    select(
      seq_id = qseqid,
      start = qstart,
      end = qend,
      seq_id2 = sseqid,
      start2 = sstart,
      end2 = send,
      pident
    )

  p_tree_base <- ggtree(tree_sorted, hang = 0.05) +
    # c(top, right, bottom, left)
    theme(plot.margin = unit(c(0.1, 0, 0.1, 0.1), "cm"))
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
    genes = genes_filtered,
  ) +
  geom_seq(aes(y = y)) +
  geom_gene(aes(y = y)) +
  geom_link(aes(y = y, fill = pident), colour = NA ) +
  geom_bin_label(aes(y = y), x = -10) +
  scale_x_continuous(
    labels = label_number(accuracy = 1),
    limits = c(0, max(ordered_seqs$length))
  ) +
  scale_fill_gradient(low = "gray80", high = "gray60") +
  # c(top, right, bottom, left)
  theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0), "cm")) 

  final_plot <- p_tree + p_genome

  if(Individual_clusters >= 2){
    ggsave(final_plot,filename=paste("CargoSynteny_Cluster",ClusterId,".svg",sep=""), width = max(16,round(max(ordered_seqs$length)*0.0001)+round(max(tree_sorted$edge.length)*10)), height = nrow(Cluster_matrix))
  } else {
    ggsave(final_plot,filename="CargoSynteny.svg", width = max(16,round(max(ordered_seqs$length)*0.0001)+round(max(tree_sorted$edge.length)*10)), height = nrow(Cluster_matrix))
  }
}
