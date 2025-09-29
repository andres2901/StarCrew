# Check and Call directory argument
suppressPackageStartupMessages(library(optparse))

if (!requireNamespace("optparse", quietly = TRUE)) {
  stop("Package \"optparse\" not installed. Please install it to run this script.", call. = FALSE)
}

option_list <- list(
  make_option(c("-d", "--directory"), type="character", action = "store", default=NULL, 
              help="path of the working directory [default %default]",metavar="PATH"),
  make_option(c("-c", "--clusters"), type="integer", action = "store", default=NULL,
              help=" Number of subclusters obtained in the syntenet approach [default %default]", metavar="number"),
  make_option(c("-g", "--gaps"), type="integer", action = "store", default=8,
              help="Number of maximum allowed gaps between anchor points for syntenet to call a collinear region [default %default]", metavar="number")
)

# Parse the command-line arguments
arguments <- parse_args(OptionParser(option_list = option_list))

# Check for required arguments and flags
if (is.null(arguments$directory)) {
  stop("Error: directory must be provided.", call.=FALSE)
}

if (is.null(arguments$clusters)) {
  stop("Error: number of clusters must be provided.", call.=FALSE)
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

OrthoFinder <- read.table("Orthogroups.GeneCount.tsv", header = T, check.names = F, row.names = 1)
OrthoFinder2 <- OrthoFinder[ , ! names(OrthoFinder) %in% c("Total")]
transposase_OrthoFinder <- t(OrthoFinder2)
distance_mat <- dist(transposase_OrthoFinder, method='binary')
Hierar_cl <- hclust(distance_mat, method = 'average')

#Identify number of clusters
fit_Orthofinder <- cutree(Hierar_cl, h = 0.99999999)
Individual_clusters <- length(unique(fit_Orthofinder))

if (Individual_clusters >= 2) {
  Clusters2 <- list()
  for( i in 1:Individual_clusters ) {
    Clusters2[i] <- "test"
  }
}

my_tree <- as.phylo(Hierar_cl)
write.tree(phy=my_tree, file="test_cluster.nw")

melted_test <- melt(transposase_test, as.is = T)
testing <- ggplot(data = melted_test, aes(x=reorder(Var2,value), y=Var1, fill=value)) + 
  geom_tile() + scale_fill_viridis(option="rocket", direction = -1) + labs(y="element",x="Orthogroup") + theme(axis.text.x=element_blank())

ggsave(testing,filename="testing_profiling.svg", width = 32, height = 16)

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

  tree_sorted <- ladderize(my_tree, right = FALSE)
num_sequences <- length(tree_sorted$tip.label)

p_tree_base <- ggtree(tree_sorted)
tree_data <- p_tree_base$data

tree_y_coords <- tree_data %>%
  filter(isTip) %>%
  select(seq_id = label, y) %>%
  mutate(y = max(y) - y + 1)

ordered_seqs <- seqs %>%
  inner_join(tree_y_coords, by = "seq_id") %>%
  arrange(y)

p_tree <- p_tree_base +
  geom_tiplab(align = TRUE, size = 3) +
  theme_tree2() +
  ylim(1, num_sequences)

p_genome <- gggenomes(
  seqs = ordered_seqs,
  links = links,
  genes = genes
) +
  geom_seq(aes(y = y)) +
  geom_gene(aes(y = y)) +
  geom_link(aes(y = y, fill = pident)) +
  geom_bin_label(aes(y = y), x = -10) +
  scale_x_continuous(
    labels = label_number(accuracy = 1),
    limits = c(0, max(ordered_seqs$length))
  )

final_plot <- p_tree + p_genome
ggsave(final_plot,filename="testing_syntenytree.svg", width = 32, height = 16)