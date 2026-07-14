# ===================================================================
# Script: 01_polb_tree_main.R
# Purpose: Build the main PolB phylogenetic tree figure (Figure 1)
# Requires: allseqs.nwk (from NuPhylo output)
# ===================================================================

library(ape)
library(tidytree)
library(treeio)
library(ggtree)
library(ggplot2)
library(dplyr)

tree <- read.tree("allseqs.nwk")
stopifnot(length(tree$tip.label) == 113)

our_candidates <- tree$tip.label[grepl("^ERR14872", tree$tip.label)]
stopifnot(length(our_candidates) == 5)

dist_matrix <- cophenetic(tree)
nearest_neighbors <- character(0)
for (cand in our_candidates) {
  dists <- dist_matrix[cand, ]
  dists <- dists[names(dists) != cand]
  dists <- dists[!(names(dists) %in% our_candidates)]
  nearest <- names(sort(dists))[1:3]
  nearest_neighbors <- c(nearest_neighbors, nearest)
}
nearest_neighbors <- unique(nearest_neighbors)

tip_data <- data.frame(label = tree$tip.label, stringsAsFactors = FALSE)
tip_data$is_ours <- tip_data$label %in% our_candidates
tip_data$display_label <- dplyr::case_when(
  tip_data$is_ours ~ gsub("_length.*", "", tip_data$label),
  tip_data$label %in% nearest_neighbors ~ tip_data$label,
  TRUE ~ ""
)
tip_data$label_type <- dplyr::case_when(
  tip_data$is_ours ~ "ours",
  tip_data$label %in% nearest_neighbors ~ "neighbor",
  TRUE ~ "other"
)

p_main <- ggtree(tree, layout = "circular", size = 0.3, color = "grey40") %<+% tip_data +
  geom_tippoint(aes(color = label_type, size = label_type)) +
  scale_color_manual(values = c(ours = "firebrick", neighbor = "black", other = "grey70"), guide = "none") +
  scale_size_manual(values = c(ours = 2.5, neighbor = 1.5, other = 0.6), guide = "none") +
  geom_tiplab(aes(label = display_label, color = label_type,
                  fontface = ifelse(label_type == "ours", "bold", "italic")),
              size = 2.5, offset = 0.05) +
  ggtitle("PolB tree: Oslo Fjord candidates (red) with 3 nearest reference neighbors each (black)") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggsave("Figure1_PolB_tree_main.pdf", plot = p_main, width = 14, height = 14, dpi = 300)
ggsave("Figure1_PolB_tree_main.png", plot = p_main, width = 14, height = 14, dpi = 300)
message("Figure 1 complete.")