#!/usr/bin/env Rscript

rm(list = ls())

library(tidyverse)
library(ggpubr)
library(vegan)
library(gridExtra)

source("../config/config.R")

# ============================================================
# Read input data
# ============================================================

sample_data <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

votu_abundance <- read.table(
  votu_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

virus_family_abundance <- read.table(
  virus_family_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ============================================================
# Check sample matching
# ============================================================

sample_data$Group <- factor(
  sample_data$Group,
  levels = c("Health", "Patient")
)

if (!all(colnames(votu_abundance) %in% sample_data$Sample)) {
  stop("Some vOTU abundance samples are absent from metadata.")
}

if (!all(colnames(virus_family_abundance) %in% sample_data$Sample)) {
  stop("Some virus-family abundance samples are absent from metadata.")
}

sample_data <- sample_data[
  match(colnames(votu_abundance), sample_data$Sample),
]

if (any(is.na(sample_data$Sample))) {
  stop("Sample matching failed.")
}

# ============================================================
# Colors
# ============================================================

custom_colors <- c(
  "Patient" = "#D44D5C",
  "Health" = "#33A02C"
)

# ============================================================
# Alpha diversity
# ============================================================

set.seed(123)

calculate_alpha_diversity(
  abundance_data = votu_abundance,
  sample_col = "Sample",
  group_col = "Group",
  group_data = sample_data,
  result_name = "alpha_re"
)

p_simpson <- ggplot(
  alpha_re,
  aes(x = Group, y = Simpson, fill = Group)
) +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  geom_signif(
    comparisons = list(c("Health", "Patient"))
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_bw()

ggsave(
  file.path(figure_dir, "otus_simpson.pdf"),
  p_simpson,
  width = 5,
  height = 3
)

p_shannon <- ggplot(
  alpha_re,
  aes(x = Group, y = Shannon, fill = Group)
) +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  geom_signif(
    comparisons = list(c("Health", "Patient"))
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_bw()

ggsave(
  file.path(figure_dir, "otus_shannon.pdf"),
  p_shannon,
  width = 5,
  height = 3
)

p_richness <- ggplot(
  alpha_re,
  aes(x = Group, y = Richness, fill = Group)
) +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(
    method = "wilcox.test",
    label = "p.signif"
  ) +
  geom_signif(
    comparisons = list(c("Health", "Patient"))
  ) +
  scale_fill_manual(values = custom_colors) +
  theme_bw()

ggsave(
  file.path(figure_dir, "otus_richness.pdf"),
  p_richness,
  width = 5,
  height = 3
)

write.table(
  alpha_re,
  file.path(table_dir, "alpha_diversity.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# ============================================================
# Beta diversity: vOTU level
# ============================================================

plot_pcoa_pro(
  abd_dt = votu_abundance,
  sample_dt = sample_data,
  colors = custom_colors,
  adonis_out_name = file.path(table_dir, "adonis_votu")
)

# ============================================================
# Beta diversity: virus family level
# ============================================================

plot_pcoa_pro(
  abd_dt = virus_family_abundance,
  sample_dt = sample_data,
  colors = custom_colors,
  adonis_out_name = file.path(table_dir, "adonis_virus_family")
)

# ============================================================
# Betadisper
# ============================================================

abd_matrix <- t(votu_abundance)

sample_data_aligned <- sample_data[
  match(rownames(abd_matrix), sample_data$Sample),
]

if (any(is.na(sample_data_aligned$Sample))) {
  stop("Sample matching failed before betadisper.")
}

dist_matrix <- vegdist(
  abd_matrix,
  method = "bray"
)

dispersion_model <- betadisper(
  dist_matrix,
  sample_data_aligned$Group
)

dispersion_test <- permutest(
  dispersion_model,
  permutations = 999
)

print(dispersion_test)

dispersion_p <- dispersion_test$tab$`Pr(>F)`[1]

cat(
  "betadisper p-value =",
  dispersion_p,
  "\n"
)

write.table(
  data.frame(
    test = "betadisper",
    p_value = dispersion_p
  ),
  file.path(table_dir, "betadisper.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)
