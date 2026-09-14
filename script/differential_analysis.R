#!/usr/bin/env Rscript

rm(list = ls())

library(tidyverse)
library(Maaslin2)

source("../config/config.R")

# ============================================================
# Input
# ============================================================

sample_data <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

votu_abundance <- read.table(
  votu_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

virus_family_abundance <- read.table(
  virus_family_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

bacterial_species_abundance <- read.table(
  bacterial_species_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE
)

rownames(sample_data) <- sample_data$Sample

sample_data$Group <- factor(
  sample_data$Group,
  levels = c("Patient", "Health")
)

# ============================================================
# Virus vOTU differential abundance
# ============================================================

virus_output <- file.path(
  result_dir,
  "maaslin_virus"
)

dir.create(
  virus_output,
  recursive = TRUE,
  showWarnings = FALSE
)

Maaslin2(
  input_data = votu_abundance,
  input_metadata = sample_data,
  output = virus_output,
  min_abundance = 0.01,
  min_prevalence = 0.01,
  normalization = "NONE",
  fixed_effects = c("Group", "Sex", "Age"),
  reference = c("Group,Patient"),
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

# ============================================================
# Virus family differential abundance
# ============================================================

family_output <- file.path(
  result_dir,
  "maaslin_virus_family"
)

dir.create(
  family_output,
  recursive = TRUE,
  showWarnings = FALSE
)

Maaslin2(
  input_data = virus_family_abundance,
  input_metadata = sample_data,
  output = family_output,
  min_abundance = 0.01,
  min_prevalence = 0.01,
  normalization = "NONE",
  fixed_effects = c("Group", "Sex", "Age"),
  reference = c("Group,Patient"),
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

# ============================================================
# Extract Group-associated bacterial features
# ============================================================

bac_result <- read.table(
  file.path(
    result_dir,
    "maaslin_bacteria",
    "significant_results.tsv"
  ),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

vir_result <- read.table(
  file.path(
    result_dir,
    "maaslin_virus",
    "significant_results.tsv"
  ),
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

sig_bac <- bac_result %>%
  filter(metadata == "Group")

sig_vir <- vir_result %>%
  filter(metadata == "Group")

# ============================================================
# Extract abundance matrices
# ============================================================

sig_bac_features <- intersect(
  sig_bac$feature,
  rownames(bacterial_species_abundance)
)

sig_vir_features <- intersect(
  sig_vir$feature,
  rownames(votu_abundance)
)

sig_bac_abundance <- bacterial_species_abundance[
  sig_bac_features,
  ,
  drop = FALSE
]

sig_vir_abundance <- votu_abundance[
  sig_vir_features,
  ,
  drop = FALSE
]

write.table(
  sig_bac_abundance,
  file.path(data_dir, "sig_bacterial_species.tsv"),
  quote = FALSE,
  sep = "\t"
)

write.table(
  sig_vir_abundance,
  file.path(data_dir, "sig_votu.tsv"),
  quote = FALSE,
  sep = "\t"
)

# ============================================================
# Export differential results
# ============================================================

write.table(
  bac_result,
  file.path(table_dir, "bacterial_maaslin_all_results.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)

write.table(
  vir_result,
  file.path(table_dir, "virus_maaslin_all_results.tsv"),
  quote = FALSE,
  sep = "\t",
  row.names = FALSE
)
