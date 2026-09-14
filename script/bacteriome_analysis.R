#!/usr/bin/env Rscript

rm(list = ls())

library(tidyverse)
library(Maaslin2)

source("../config/config.R")

# ============================================================
# Read data
# ============================================================

sample_data <- read.table(
  metadata_file,
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

bacterial_species_abundance <- read.table(
  bacterial_species_file,
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(sample_data) <- sample_data$Sample

# ============================================================
# Check samples
# ============================================================

missing_samples <- setdiff(
  colnames(bacterial_species_abundance),
  rownames(sample_data)
)

if (length(missing_samples) > 0) {
  stop(
    "Samples in bacterial abundance table are missing from metadata: ",
    paste(missing_samples, collapse = ", ")
  )
}

sample_data <- sample_data[
  colnames(bacterial_species_abundance),
  ,
  drop = FALSE
]

# ============================================================
# MaAsLin2
# ============================================================

dir.create(
  file.path(result_dir, "maaslin_bacteria"),
  recursive = TRUE,
  showWarnings = FALSE
)

fit_bacteria <- Maaslin2(
  input_data = bacterial_species_abundance,
  input_metadata = sample_data,
  output = file.path(result_dir, "maaslin_bacteria"),
  min_abundance = 0.01,
  min_prevalence = 0.01,
  normalization = "NONE",
  fixed_effects = c("Group", "Sex", "Age"),
  reference = c("Group,Patient"),
  plot_heatmap = TRUE,
  plot_scatter = TRUE
)

cat("Bacterial MaAsLin2 analysis completed.\n")
