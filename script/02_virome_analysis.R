# PTB virome project
# 02_virome_analysis.R
# Alpha diversity, beta diversity, and dispersion testing.

rm(list = ls())

source("script/config.R")

library(tidyverse)
library(ggpubr)
library(vegan)
library(gridExtra)

source("script/calculate_alpha.R")
source("script/plot_aplha.R")
source("script/plot_beta.R")

metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)
virus_family <- read.delim(VIRUS_FAMILY_FILE, check.names = FALSE, row.names = 1)

metadata$Group <- recode(metadata$Group,
                          "Test-P" = "Patient",
                          "PTB" = "Patient",
                          "TB" = "Patient",
                          "HC" = "Health")

custom_colors <- c("Patient" = "#D44D5C", "Health" = "#33a02c")

# Alpha diversity
set.seed(123)
calculate_alpha_diversity(
  abundance_data = votu,
  sample_col = "Sample",
  group_col = "Group",
  group_data = metadata,
  result_name = "alpha_re"
)

p_simpson <- ggplot(alpha_re, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  geom_signif(comparisons = list(c("Health", "Patient"))) +
  scale_fill_manual(values = custom_colors) +
  theme_bw()
ggsave("results/figures/otus_simpson.pdf", p_simpson, width = 5, height = 3)

p_shannon <- ggplot(alpha_re, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  geom_signif(comparisons = list(c("Health", "Patient"))) +
  scale_fill_manual(values = custom_colors) +
  theme_bw()
ggsave("results/figures/otus_shannon.pdf", p_shannon, width = 5, height = 3)

p_obs <- ggplot(alpha_re, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot(alpha = 0.6) +
  stat_compare_means(method = "wilcox.test", label = "p.signif") +
  geom_signif(comparisons = list(c("Health", "Patient"))) +
  scale_fill_manual(values = custom_colors) +
  theme_bw()
ggsave("results/figures/otus_richness.pdf", p_obs, width = 5, height = 3)

# Beta diversity
plot_pcoa_pro(votu, metadata, colors = custom_colors, adonis_out_name = "adonis_votu")
plot_pcoa_pro(virus_family, metadata, colors = custom_colors, adonis_out_name = "adonis_family")

# Betadisper for vOTU-level Bray-Curtis distances.
abd_matrix <- if (grepl("^v\\d+", rownames(votu)[1])) t(votu) else votu
metadata_aligned <- metadata[match(rownames(abd_matrix), metadata$Sample), , drop = FALSE]
if (any(is.na(metadata_aligned$Sample))) stop("Sample alignment failed for betadisper.")

dist_matrix <- vegdist(abd_matrix, method = "bray")
mod <- betadisper(dist_matrix, metadata_aligned$Group)
disp_test <- permutest(mod, permutations = 999)
print(disp_test)
p_val <- disp_test$tab$`Pr(>F)`[1]
cat("betadisper p-value =", p_val, "\n")
