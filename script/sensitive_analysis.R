# ==============================================================================
# PTB Virome Project
# MaAsLin2 sensitivity analysis: alternative abundance / prevalence thresholds
#
# Primary analysis : min_abundance = 0.01, min_prevalence = 0.01
# Sensitivity      : one parameter changed at a time, all other MaAsLin2
#                    parameters identical to the primary analysis (config.R).
#
# Input  : data/metadata.tsv, data/votu_abundance.tsv
# Output : results/tables/maaslin_sensitivity_summary.tsv
#          results/tables/maaslin_sensitivity_overlap.tsv
#          results/tables/maaslin_sensitivity_core_vOTUs.tsv
#          results/maaslin/maaslin_output_<label>/  (full MaAsLin2 output)
#
# No figures are produced by this script.
# ==============================================================================

rm(list = ls())

source("script/config.R")

suppressPackageStartupMessages({
  library(Maaslin2)
  library(dplyr)
  library(tidyr)
})

dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. Read input data
# ------------------------------------------------------------------------------
metadata <- read.delim(
  METADATA_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

votu <- read.delim(
  VOTU_FILE,
  check.names = FALSE,
  row.names = 1
)

stopifnot(all(metadata$Group %in% GROUP_LEVELS))
rownames(metadata) <- metadata$Sample

# ------------------------------------------------------------------------------
# 2. Sensitivity-analysis parameter grid
#
#    Primary analysis: ma = 0.01, mp = 0.01
#    Each sensitivity run changes only ONE parameter relative to the primary.
# ------------------------------------------------------------------------------
param_grid <- data.frame(
  label          = c("Primary_ma001_mp001",
                     "A_ma0005_mp001",
                     "B_ma002_mp001",
                     "C_ma001_mp005",
                     "D_ma001_mp010",
                     "E_ma001_mp020"),
  min_abundance  = c(0.01, 0.005, 0.02, 0.01, 0.01, 0.01),
  min_prevalence = c(0.01, 0.01,  0.01, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------------------
# 3. Run MaAsLin2 for each parameter combination
# ------------------------------------------------------------------------------
for (i in seq_len(nrow(param_grid))) {
  label <- param_grid$label[i]
  ma    <- param_grid$min_abundance[i]
  mp    <- param_grid$min_prevalence[i]
  
  output_dir <- file.path(MAASLIN_DIR, paste0("maaslin_output_", label))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  message("============================================================")
  message("Running MaAsLin2: ", label)
  message("  min_abundance  = ", ma)
  message("  min_prevalence = ", mp)
  message("  output         = ", output_dir)
  message("============================================================")
  
  Maaslin2(
    input_data     = votu,
    input_metadata = metadata,
    output         = output_dir,
    min_abundance  = ma,
    min_prevalence = mp,
    normalization  = MAASLIN_NORMALIZATION,
    fixed_effects  = MAASLIN_FIXED_EFFECTS,
    reference      = MAASLIN_REFERENCE,
    plot_heatmap   = FALSE,
    plot_scatter   = FALSE
  )
}

# ------------------------------------------------------------------------------
# 4. Collect significant vOTUs (metadata == "Group") from each run
# ------------------------------------------------------------------------------
extract_sig_group <- function(result_file) {
  if (!file.exists(result_file)) return(character(0))
  df <- read.delim(result_file, stringsAsFactors = FALSE, check.names = FALSE)
  df %>%
    filter(metadata == "Group") %>%
    filter(qval < MAASLIN_Q_THRESHOLD) %>%
    pull(feature) %>%
    unique()
}

sig_list <- lapply(param_grid$label, function(label) {
  rf <- file.path(MAASLIN_DIR, paste0("maaslin_output_", label),
                  "significant_results.tsv")
  extract_sig_group(rf)
})
names(sig_list) <- param_grid$label

# ------------------------------------------------------------------------------
# 5. Summary table: number of significant vOTUs per threshold
# ------------------------------------------------------------------------------
primary_label <- "Primary_ma001_mp001"
primary_set   <- sig_list[[primary_label]]

summary_df <- data.frame(
  Analysis       = param_grid$label,
  min_abundance  = param_grid$min_abundance,
  min_prevalence = param_grid$min_prevalence,
  n_significant  = sapply(sig_list, length),
  n_overlap_with_primary = sapply(sig_list, function(x) length(intersect(x, primary_set))),
  pct_overlap_with_primary = round(
    100 * sapply(sig_list, function(x) length(intersect(x, primary_set))) /
      length(primary_set), 1
  ),
  stringsAsFactors = FALSE
)

write.table(
  summary_df,
  file = file.path(TABLE_DIR, "maaslin_sensitivity_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ------------------------------------------------------------------------------
# 6. Pairwise overlap matrix (Jaccard / shared counts)
# ------------------------------------------------------------------------------
all_labels <- names(sig_list)
n <- length(all_labels)

shared_mat <- matrix(0L, nrow = n, ncol = n,
                     dimnames = list(all_labels, all_labels))
for (i in seq_len(n)) {
  for (j in seq_len(n)) {
    shared_mat[i, j] <- length(intersect(sig_list[[i]], sig_list[[j]]))
  }
}

shared_df <- data.frame(
  Analysis = rownames(shared_mat),
  as.data.frame(shared_mat),
  check.names = FALSE
)

write.table(
  shared_df,
  file = file.path(TABLE_DIR, "maaslin_sensitivity_overlap.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ------------------------------------------------------------------------------
# 7. Core vOTUs significant across ALL thresholds
# ------------------------------------------------------------------------------
core_features <- Reduce(intersect, sig_list)

core_df <- data.frame(
  feature = core_features,
  stringsAsFactors = FALSE
)

write.table(
  core_df,
  file = file.path(TABLE_DIR, "maaslin_sensitivity_core_vOTUs.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ------------------------------------------------------------------------------
# 8. Console summary
# ------------------------------------------------------------------------------
message("============================================================")
message("Sensitivity analysis summary")
message("============================================================")
print(summary_df)
message("------------------------------------------------------------")
message("Core vOTUs significant across ALL thresholds: ",
        length(core_features))
message("------------------------------------------------------------")
message("Tables written to: ", TABLE_DIR)
message("============================================================")
