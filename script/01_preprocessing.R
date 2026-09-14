# PTB virome project
# 01_preprocessing.R
#
# Final GitHub input tables:
#   data/metadata.tsv
#   data/votu_abundance.tsv
#   data/virus_family_abundance.tsv
#   data/bacterial_species_abundance.tsv
#   data/votu_ko_abundance.tsv
#
# This script validates the common sample identifiers and standardized metadata.

rm(list = ls())

source("script/config.R")

library(tidyverse)

metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)
virus_family <- read.delim(VIRUS_FAMILY_FILE, check.names = FALSE, row.names = 1)
bacteria <- read.delim(BACTERIAL_FILE, check.names = FALSE, row.names = 1)
ko_abundance <- read.delim(KO_ABUNDANCE_FILE, check.names = FALSE, row.names = 1)

required_metadata <- c("Sample", "Group", "Sex", "Age")
stopifnot(all(required_metadata %in% colnames(metadata)))

metadata$Group <- recode(metadata$Group,
                          "Test-P" = "Patient",
                          "PTB" = "Patient",
                          "TB" = "Patient",
                          "HC" = "Health")
metadata$Group <- factor(metadata$Group, levels = c("Health", "Patient"))

if (anyDuplicated(metadata$Sample)) {
  stop("Duplicate sample IDs detected in data/metadata.tsv")
}

abundance_tables <- list(
  vOTU = votu,
  virus_family = virus_family,
  bacteria = bacteria
)

for (nm in names(abundance_tables)) {
  missing_samples <- setdiff(colnames(abundance_tables[[nm]]), metadata$Sample)
  if (length(missing_samples) > 0) {
    stop(sprintf("Samples in %s abundance table are missing from metadata.tsv: %s",
                nm, paste(missing_samples, collapse = ", ")))
  }
}

# Keep the sample order defined by metadata.
common_samples <- metadata$Sample
votu <- votu[, intersect(common_samples, colnames(votu)), drop = FALSE]
virus_family <- virus_family[, intersect(common_samples, colnames(virus_family)), drop = FALSE]
bacteria <- bacteria[, intersect(common_samples, colnames(bacteria)), drop = FALSE]

write.table(metadata, "results/tables/metadata_standardized.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("Input validation completed.\n")
cat("Samples:", nrow(metadata), "\n")
cat("Health:", sum(metadata$Group == "Health"), "\n")
cat("Patient:", sum(metadata$Group == "Patient"), "\n")
