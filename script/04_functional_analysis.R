# PTB virome project
# 04_functional_analysis.R
# Functional enrichment of KO annotations.
#
# IMPORTANT: The current manuscript code performs the Fisher test on the
# number of vOTUs carrying each KO. The five-input GitHub dataset also
# contains votu_ko_abundance.tsv, but the exact row/column definition of
# that table has not been specified here. Therefore this script does not
# silently reinterpret that table as a vOTU-to-KO mapping.
#
# To reproduce the current analysis exactly, a vOTU -> KO mapping table is
# still required. Place it at data/annotations/votu_ko_mapping.tsv with:
#   vOTU    KO
#
# The Fisher p-values are intentionally UNADJUSTED here, matching the code
# supplied for this analysis. Do not add BH correction unless the analysis
# is explicitly changed.

rm(list = ls())

source("script/config.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)

metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)

metadata$Group <- recode(metadata$Group,
                          "Test-P" = "Patient",
                          "PTB" = "Patient",
                          "TB" = "Patient",
                          "HC" = "Health")

# Use the MaAsLin2 significant vOTUs from the differential analysis.
co_re <- read.delim("results/maaslin/vir/significant_results.tsv") %>%
  filter(metadata == "Group") %>%
  mutate(Group = ifelse(coef > 0, "HC", "TB"))

mapping_votu_ko <- read.delim(VOTU_KO_MAPPING_FILE,
                              check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(all(c("vOTU", "KO") %in% colnames(mapping_votu_ko)))

combined_data <- co_re %>%
  select(feature, Group) %>%
  inner_join(mapping_votu_ko, by = c("feature" = "vOTU")) %>%
  distinct(feature, Group, KO)

# The original analysis used the total number of differential vOTUs in each
# enrichment direction as the denominator.
total_votu <- co_re %>%
  distinct(feature, Group) %>%
  count(Group, name = "total_votu")

ko_counts <- combined_data %>%
  group_by(Group, KO) %>%
  summarise(count = n_distinct(feature), .groups = "drop")

merged_table <- ko_counts %>%
  pivot_wider(names_from = Group, values_from = count,
              values_fill = list(count = 0))

if (!all(c("HC", "TB") %in% colnames(merged_table))) {
  if (!"HC" %in% colnames(merged_table)) merged_table$HC <- 0L
  if (!"TB" %in% colnames(merged_table)) merged_table$TB <- 0L
}

n_hc <- total_votu$total_votu[match("HC", total_votu$Group)]
n_tb <- total_votu$total_votu[match("TB", total_votu$Group)]

merged_table <- merged_table %>%
  rename(count_in_HC = HC, count_in_TB = TB) %>%
  mutate(
    count_not_in_HC = n_hc - count_in_HC,
    count_not_in_TB = n_tb - count_in_TB
  )

fisher_test <- function(a, b, c, d) {
  fisher.test(matrix(c(a, b, c, d), nrow = 2))$p.value
}

merged_table <- merged_table %>%
  mutate(
    fisher_p_value = mapply(
      fisher_test,
      count_in_HC, count_in_TB,
      count_not_in_HC, count_not_in_TB
    ),
    occurrence_rate_in_HC = count_in_HC / n_hc,
    occurrence_rate_in_TB = count_in_TB / n_tb,
    diff = occurrence_rate_in_TB - occurrence_rate_in_HC,
    enriched_in = ifelse(diff > 0, "TB", "HC"),
    signif_label = case_when(
      fisher_p_value < 0.001 ~ "***",
      fisher_p_value < 0.01 ~ "**",
      fisher_p_value < 0.05 ~ "*",
      TRUE ~ ""
    )
  )

# KEGG descriptions are intentionally not hard-coded here. Add a mapping
# table under data/annotations/ko_description.tsv if available, with columns
# KO and Description.
if (file.exists("data/annotations/ko_description.tsv")) {
  ko_description <- read.delim("data/annotations/ko_description.tsv",
                                check.names = FALSE)
  merged_table <- merged_table %>% left_join(ko_description, by = "KO")
} else {
  merged_table$Description <- NA_character_
}

final_table <- merged_table %>%
  filter(fisher_p_value < 0.05) %>%
  mutate(Gene_Description = ifelse(
    is.na(Description), KO, paste(KO, Description, sep = ": ")
  ))

write.csv(final_table, "results/tables/KO_enriched_results.csv",
          row.names = FALSE, quote = FALSE)

# Plot 1: occurrence rate.
diff_data <- final_table %>% arrange(desc(diff))
sorted_genes <- unique(diff_data$Gene_Description)

long_data <- pivot_longer(
  diff_data,
  cols = c(occurrence_rate_in_HC, occurrence_rate_in_TB),
  names_to = "group",
  values_to = "occurrence_rate"
)
long_data$Gene_Description <- factor(long_data$Gene_Description,
                                      levels = sorted_genes)

custom_colors <- c(
  "occurrence_rate_in_HC" = "#2e8153",
  "occurrence_rate_in_TB" = "#bc5858"
)

num_levels <- length(levels(long_data$Gene_Description))
x_range <- c(0, max(long_data$occurrence_rate, na.rm = TRUE) * 1.05)

p_occurrence <- ggplot(
  long_data,
  aes(x = occurrence_rate,
      y = as.numeric(Gene_Description) + 0.5,
      color = group)
) +
  geom_point(size = 2) +
  scale_x_continuous(labels = percent, limits = x_range) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    breaks = seq_along(levels(long_data$Gene_Description)) + 0.5,
    labels = levels(long_data$Gene_Description),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  labs(x = "Occurrence Rate", y = "Gene: Description", color = "Group") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    panel.grid = element_blank(),
    legend.position = "top"
  )

ggsave("results/figures/kegg_occurrence.pdf", p_occurrence,
       width = 10, height = max(3, num_levels * 0.25 + 2), dpi = 300)

# Plot 2: volcano plot. Fisher p-values remain unadjusted.
volcano_data <- merged_table %>%
  mutate(
    log2FC = log2((occurrence_rate_in_TB + 1e-6) /
                  (occurrence_rate_in_HC + 1e-6)),
    neg_log10_p = -log10(fisher_p_value),
    enrichment = case_when(
      fisher_p_value < 0.05 & log2FC > 0 ~ "Enriched in TB",
      fisher_p_value < 0.05 & log2FC < 0 ~ "Enriched in HC",
      TRUE ~ "Not Significant"
    ),
    significant = fisher_p_value < 0.05
  )

p_kegg_volcano <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_p)) +
  geom_point(aes(color = enrichment), size = 1.5, alpha = 0.7) +
  geom_point(
    data = subset(volcano_data, significant),
    aes(fill = enrichment), shape = 21, color = "black", size = 3
  ) +
  geom_text_repel(
    data = subset(volcano_data, significant),
    aes(label = KO), size = 3, box.padding = 0.3, max.overlaps = 100
  ) +
  scale_color_manual(values = c(
    "Enriched in TB" = "#bc5858",
    "Enriched in HC" = "#2e8153",
    "Not Significant" = "grey70"
  )) +
  scale_fill_manual(values = c(
    "Enriched in TB" = "#bc5858",
    "Enriched in HC" = "#2e8153",
    "Not Significant" = "grey70"
  )) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "log2(Fold Change)", y = "-log10(p-value)",
       color = "Enrichment", fill = "Enrichment") +
  theme_minimal(base_size = 13)

ggsave("results/figures/kegg_volcano.pdf", p_kegg_volcano,
       height = 5, width = 7)
