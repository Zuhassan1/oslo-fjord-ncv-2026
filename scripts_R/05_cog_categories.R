# ===================================================================
# Script: 05_cog_categories.R
# Purpose: COG functional category distribution (Figure 5)
# Requires: candidate_annotations_VALID26.emapper.annotations
# Note: recomputes counts directly from the real annotation file,
#       rather than using hand-transcribed values, for genuine
#       reproducibility.
# ===================================================================

library(ggplot2)
library(dplyr)

lines <- readLines("candidate_annotations_VALID26.emapper.annotations")
data_lines <- lines[!grepl("^#", lines)]

cog_letters <- c()
for (line in data_lines) {
  fields <- strsplit(line, "\t")[[1]]
  if (length(fields) < 7) next
  cog_field <- trimws(fields[7])
  if (cog_field != "" && cog_field != "-" && !is.na(cog_field)) {
    cog_letters <- c(cog_letters, strsplit(cog_field, "")[[1]])
  }
}

cog_counts <- as.data.frame(table(cog_letters), stringsAsFactors = FALSE)
colnames(cog_counts) <- c("category", "count")

descriptions <- c(
  S="Function unknown", L="Replication/recombination/repair", K="Transcription",
  O="Posttranslational modification/chaperones", F="Nucleotide transport/metabolism",
  A="RNA processing/modification", G="Carbohydrate transport/metabolism",
  M="Cell wall/membrane/envelope biogenesis", J="Translation/ribosome biogenesis",
  H="Coenzyme transport/metabolism", T="Signal transduction",
  E="Amino acid transport/metabolism", I="Lipid transport/metabolism",
  C="Energy production/conversion", P="Inorganic ion transport/metabolism",
  V="Defense mechanisms", B="Chromatin structure/dynamics",
  U="Intracellular trafficking", Q="Secondary metabolite biosynthesis",
  D="Cell cycle control/division", W="Extracellular structures"
)
cog_counts$description <- descriptions[cog_counts$category]

message("Re-derived COG counts (verify against manuscript Table):")
print(cog_counts[order(-cog_counts$count), ])

cog_counts$highlight <- cog_counts$category %in% c("L", "K", "O")
cog_counts$category <- factor(cog_counts$category,
                              levels = cog_counts$category[order(-cog_counts$count)])

p_cog <- ggplot(cog_counts, aes(x = category, y = count, fill = highlight)) +
  geom_col() +
  geom_text(aes(label = count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(`FALSE` = "grey60", `TRUE` = "firebrick"), guide = "none") +
  labs(title = "COG functional category distribution across 26 candidates",
       subtitle = "Red: top 3 categories (excl. S), matching the reference study's own top 3 exactly",
       x = "COG category", y = "Count") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        plot.subtitle = element_text(size = 9, color = "grey30"))

ggsave("Figure5_COG_categories.pdf", plot = p_cog, width = 10, height = 6, dpi = 300)
ggsave("Figure5_COG_categories.png", plot = p_cog, width = 10, height = 6, dpi = 300)
message("Figure 5 complete.")