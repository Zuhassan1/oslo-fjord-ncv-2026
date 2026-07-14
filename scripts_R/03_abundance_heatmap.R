# ===================================================================
# Script: 03_abundance_heatmap.R
# Purpose: Heatmap of NCV candidate abundance across 9 samples (Figure 3)
# Requires: abundance_matrix_RPM.tsv
# ===================================================================

library(pheatmap)
library(dplyr)

abund <- read.table("abundance_matrix_RPM.tsv", header = TRUE, sep = "\t",
                    row.names = 1, check.names = FALSE)

# Candidates excluded from abundance interpretation per taxonomic assignment
# results (see manuscript "Taxonomic assignment" and "Abundance mapping" sections)
exclude <- c(
  "ERR14872777_NODE_249_length_32110_cov_20.287506",
  "ERR14872779_NODE_157_length_56008_cov_8.994245",
  "ERR14872777_NODE_117_length_50263_cov_15.779637",
  "ERR14872777_NODE_216_length_35230_cov_18.204805"
)
abund_22 <- abund[!(rownames(abund) %in% exclude), ]
stopifnot(nrow(abund_22) == 22)

mat_log <- log10(as.matrix(abund_22) + 1)
colnames(mat_log) <- gsub("_RPM", "", colnames(mat_log))
rownames(mat_log) <- gsub("_length.*", "", rownames(mat_log))

sample_meta <- data.frame(
  Site = c("Operastranda","Operastranda","Operastranda","Operastranda",
           "Operastranda","Operastranda","Operastranda","Huk","Huk"),
  Category = c("Flood","Non-flood","Non-flood","Non-flood",
               "Non-flood","Non-flood","Non-flood","Control","Control"),
  row.names = colnames(mat_log)
)

pheatmap(
  mat_log, cluster_rows = TRUE, cluster_cols = FALSE,
  annotation_col = sample_meta,
  color = colorRampPalette(c("white", "#4575b4", "#313695"))(100),
  main = "NCV candidate abundance across 9 samples (log10 RPM+1)",
  fontsize_row = 7, fontsize_col = 9,
  filename = "Figure3_abundance_heatmap.png", width = 10, height = 8
)
pheatmap(
  mat_log, cluster_rows = TRUE, cluster_cols = FALSE,
  annotation_col = sample_meta,
  color = colorRampPalette(c("white", "#4575b4", "#313695"))(100),
  main = "NCV candidate abundance across 9 samples (log10 RPM+1)",
  fontsize_row = 7, fontsize_col = 9,
  filename = "Figure3_abundance_heatmap.pdf", width = 10, height = 8
)
message("Figure 3 complete.")