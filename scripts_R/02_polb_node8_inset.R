# ===================================================================
# Script: 02_polb_node8_inset.R
# Purpose: Zoomed inset showing NODE_8's placement relative to
#          Asfarviridae reference isolates (Figure 2)
# Requires: allseqs.nwk
# ===================================================================

library(ape)
library(ggtree)
library(ggplot2)
library(dplyr)

tree <- read.tree("allseqs.nwk")

asfv_tips <- tree$tip.label[grepl("Asfarviridae", tree$tip.label)]
stopifnot(length(asfv_tips) == 3)

mrca_node <- getMRCA(tree, c("ERR14872777_NODE_8_length_176732_cov_16.217617", asfv_tips))
node8_clade <- extract.clade(tree, mrca_node)
stopifnot(length(node8_clade$tip.label) == 7)

sub_data <- data.frame(label = node8_clade$tip.label, stringsAsFactors = FALSE)
sub_data$is_ours <- grepl("^ERR14872", sub_data$label)
sub_data$is_asfv <- grepl("Asfarviridae", sub_data$label)
sub_data$display_label <- sub_data$label
sub_data$display_label[sub_data$is_ours] <- gsub("_length.*", "", sub_data$display_label[sub_data$is_ours])
sub_data$display_label[sub_data$is_asfv] <- gsub("Asfarviridae_", "ASFV_", sub_data$display_label[sub_data$is_asfv])
sub_data$point_color <- case_when(
  sub_data$is_ours ~ "ours",
  sub_data$is_asfv ~ "asfv",
  TRUE ~ "other"
)

p_inset <- ggtree(node8_clade, layout = "rectangular", size = 0.6) %<+% sub_data +
  geom_tippoint(aes(color = point_color), size = 3) +
  scale_color_manual(values = c(ours = "firebrick", asfv = "steelblue", other = "grey50"), guide = "none") +
  geom_tiplab(aes(label = display_label, color = point_color,
                  fontface = ifelse(point_color == "ours", "bold", "plain")),
              size = 3.5, offset = 0.01) +
  xlim(0, max(node8_clade$edge.length, na.rm = TRUE) * 8) +
  ggtitle("NODE_8's clade: placement alongside Asfarviridae (blue) reference isolates") +
  theme(plot.title = element_text(size = 11, face = "bold", hjust = 0.5))

ggsave("Figure2_NODE8_ASFV_inset.pdf", plot = p_inset, width = 8, height = 5, dpi = 300)
ggsave("Figure2_NODE8_ASFV_inset.png", plot = p_inset, width = 8, height = 5, dpi = 300)
message("Figure 2 complete.")