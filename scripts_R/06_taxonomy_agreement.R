# ===================================================================
# Script: 06_taxonomy_agreement_CORRECTED.R
# Purpose: TIGTOG vs. GVClass agreement summary (Figure 6) + detail table
# CORRECTED VERSION: uses GVClass's own taxonomy_majority column
# (its own designed consensus-confidence field) rather than extracting
# the first NCLDV-prefixed substring from the raw percentage breakdown,
# which incorrectly credited GVClass with an order-level call it did
# not actually make for several candidates (e.g. NODE_8, NODE_78).
# Requires: tigtog_results_v5.prediction_result.tsv, gvclass_summary.tsv
# ===================================================================

library(dplyr)
library(stringr)
library(ggplot2)

tigtog <- read.table("tigtog_results_v5.prediction_result.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
gvclass <- read.table("gvclass_summary.tsv", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
stopifnot(nrow(tigtog) == 26, nrow(gvclass) == 26)

tigtog_clean <- tigtog %>%
  transmute(candidate = str_replace(Sequence, "_length.*", ""), tigtog_order = Predicted_Order)

# Extract the order-level ('o_') component specifically from GVClass's own
# taxonomy_majority column. A blank value after 'o_' means GVClass itself
# did not reach a confident order-level majority call for that candidate.
extract_order <- function(taxonomy_majority_string) {
  parts <- str_split(taxonomy_majority_string, ";")[[1]]
  order_part <- parts[str_starts(parts, "o_")]
  order_val <- str_remove(order_part, "^o_")
  if (length(order_val) == 0 || order_val == "" || is.na(order_val)) {
    return("No confident order-level call")
  }
  return(order_val)
}

gvclass_clean <- gvclass %>%
  rowwise() %>%
  mutate(
    candidate = str_replace(query, "_length.*", ""),
    gvclass_order_raw = order,                       # raw percentage breakdown, kept for reference
    gvclass_majority_raw = taxonomy_majority,          # GVClass's own consensus field, kept for reference
    gvclass_order = extract_order(taxonomy_majority)
  ) %>%
  ungroup() %>%
  select(candidate, gvclass_order, gvclass_order_raw, gvclass_majority_raw)

comparison <- inner_join(tigtog_clean, gvclass_clean, by = "candidate")
stopifnot(nrow(comparison) == 26)

comparison$agreement <- case_when(
  comparison$tigtog_order == "Not_GV" & comparison$gvclass_order == "No confident order-level call" ~ "Both: no signal",
  comparison$gvclass_order == "No confident order-level call" ~ "GVClass: no confident order call",
  comparison$tigtog_order == comparison$gvclass_order ~ "Agree",
  TRUE ~ "Disagree"
)

message("CORRECTED agreement summary (using GVClass's own taxonomy_majority field):")
print(table(comparison$agreement))

message("\nFull candidate-by-candidate table:")
print(comparison %>% select(candidate, tigtog_order, gvclass_order, agreement))

comparison$agreement <- factor(comparison$agreement,
  levels = c("Agree", "Both: no signal", "GVClass: no confident order call", "Disagree"))
summary_counts <- comparison %>% count(agreement)

p_agree <- ggplot(summary_counts, aes(x = agreement, y = n, fill = agreement)) +
  geom_col() +
  geom_text(aes(label = n), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Agree" = "#2ca25f", "Both: no signal" = "#969696",
    "GVClass: no confident order call" = "#e6550d", "Disagree" = "#fdae61"
  ), guide = "none") +
  labs(title = "TIGTOG vs. GVClass taxonomic agreement across 26 candidates (corrected)",
       subtitle = "GVClass calls based on its own taxonomy_majority consensus field",
       x = NULL, y = "Number of candidates") +
  ylim(0, max(summary_counts$n) * 1.15) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("Figure6_taxonomy_agreement_CORRECTED.pdf", plot = p_agree, width = 8, height = 6, dpi = 300)
ggsave("Figure6_taxonomy_agreement_CORRECTED.png", plot = p_agree, width = 8, height = 6, dpi = 300)

write.csv(
  comparison %>% select(candidate, tigtog_order, gvclass_order, gvclass_order_raw, gvclass_majority_raw, agreement) %>% arrange(agreement, candidate),
  "Table_taxonomy_comparison_detail_CORRECTED.csv", row.names = FALSE
)
message("\nDone. Compare against the original Figure 6 (15 Agree / 2 no-signal / 1 non-NCLDV / 8 Disagree).")
