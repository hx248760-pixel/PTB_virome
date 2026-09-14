# PTB virome project
# 03_differential_analysis.R
# Differential viral and bacterial features, volcano plots, family composition,
# and predicted viral-host genus plots.

rm(list = ls())

source("script/config.R")

library(tidyverse)
library(Maaslin2)
library(ggplot2)

metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)
virus_family <- read.delim(VIRUS_FAMILY_FILE, check.names = FALSE, row.names = 1)
bacteria <- read.delim(BACTERIAL_FILE, check.names = FALSE, row.names = 1)

metadata$Group <- recode(metadata$Group,
                          "Test-P" = "Patient",
                          "PTB" = "Patient",
                          "TB" = "Patient",
                          "HC" = "Health")
rownames(metadata) <- metadata$Sample

# MaAsLin2 settings are retained from the original analysis.
Maaslin2(
  votu, metadata, "results/maaslin/vir",
  min_abundance = 0.01, min_prevalence = 0.01,
  normalization = "NONE",
  fixed_effects = c("Group", "Sex", "Age"),
  reference = c("Group,Patient"),
  plot_heatmap = TRUE, plot_scatter = TRUE
)

vir_results <- read.delim("results/maaslin/vir/all_results.tsv")
write.table(vir_results, "results/tables/sig_vir_otu.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

Maaslin2(
  virus_family, metadata, "results/maaslin/vir_family",
  min_abundance = 0.01, min_prevalence = 0.01,
  normalization = "NONE",
  fixed_effects = c("Group", "Sex", "Age"),
  reference = c("Group,Patient"),
  plot_heatmap = TRUE, plot_scatter = TRUE
)

family_results <- read.delim("results/maaslin/vir_family/all_results.tsv")
write.table(family_results, "results/tables/sig_vir_family.txt",
            quote = FALSE, row.names = FALSE, sep = "\t")

# Volcano plot: exploratory q < 0.25 threshold.
volcano_data <- vir_results %>%
  filter(metadata == "Group") %>%
  mutate(Group = case_when(
    coef > 0 & qval < 0.25 ~ "Health",
    coef < 0 & qval < 0.25 ~ "Patient",
    TRUE ~ "Not significant"
  ))

p_volcano <- ggplot(volcano_data, aes(x = coef, y = -log10(qval), color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Health" = "#33a02c",
    "Patient" = "#D44D5C",
    "Not significant" = "#AAAAAA"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.25), linetype = "dashed") +
  labs(x = "Coefficient (coef)", y = "-log10(q value)",
       title = "Volcano Plot of vOTU Features", color = "Group") +
  theme_bw()

ggsave("results/figures/votu_volcano_q0.25.pdf", p_volcano, width = 7, height = 5)

# The following host/family plotting requires annotation tables that are not
# among the five primary abundance tables and the user-supplied helper
# plot_virf_host.R. These files are expected under data/annotations/.
# Required examples:
#   data/annotations/vOTU.family.tsv
#   data/annotations/vOTU.host.tax.genus.tsv
#
# plot_virf_host.R will be added after its plotting function is finalized.
