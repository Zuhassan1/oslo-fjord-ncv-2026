# ===================================================================
# Script: 04_huk_persistence_panel.R
# Purpose: Bar chart of the two ANI-confirmed Huk-persistent genome
#          pairs across all 9 samples (Figure 4)
# Requires: abundance_matrix_RPM.tsv
# ===================================================================

library(ggplot2)
library(dplyr)

abund <- read.table("abundance_matrix_RPM.tsv", header = TRUE, sep = "\t",
                    row.names = 1, check.names = FALSE)

huk_candidates <- c(
  "ERR14872781_NODE_137_length_51994_cov_10.210324",
  "ERR14872781_NODE_192_length_41414_cov_9.182886"
)
stopifnot(all(huk_candidates %in% rownames(abund)))
check <- abund[huk_candidates, ]

huk_pairs <- data.frame(
  candidate = rep(c("NODE_137 / NODE_68", "NODE_192 / NODE_108"), each = 9),
  sample = rep(c("T1(flood)","T2","T3","T4","T5","T6","T7","Huk-1","Huk-2"), times = 2),
  rpm = c(
    as.numeric(check["ERR14872781_NODE_137_length_51994_cov_10.210324", ]),
    as.numeric(check["ERR14872781_NODE_192_length_41414_cov_9.182886", ])
  )
)
huk_pairs$sample <- factor(huk_pairs$sample,
                           levels = c("T1(flood)","T2","T3","T4","T5","T6","T7","Huk-1","Huk-2"))
huk_pairs$is_huk <- huk_pairs$sample %in% c("Huk-1","Huk-2")

p_huk <- ggplot(huk_pairs, aes(x = sample, y = rpm, fill = is_huk)) +
  geom_col() +
  facet_wrap(~candidate, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "firebrick"), guide = "none") +
  labs(title = "Site-persistent NCV pairs: near-identical genomes (99.93% / 99.95% ANI)\nconsistently recovered at Huk across two sampling dates",
       x = NULL, y = "Abundance (RPM)") +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(size = 11, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figure4_Huk_persistence.pdf", plot = p_huk, width = 8, height = 6, dpi = 300)
ggsave("Figure4_Huk_persistence.png", plot = p_huk, width = 8, height = 6, dpi = 300)
message("Figure 4 complete.")