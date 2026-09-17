# ==============================================================================
# PTB Virome Project
# 01_bac_analysis.R
# Differential abundance analysis for bacterial species using MaAsLin2.
# ==============================================================================

rm(list = ls())

source("script/config.R")

library(Maaslin2)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
# 读取元数据
metadata <- read.delim(METADATA_FILE, row.names = 1, check.names = FALSE)

# 读取细菌丰度表（首列为 name，自动读取并转置）
bac_raw <- read.delim(BACTERIAL_FILE, check.names = FALSE, stringsAsFactors = FALSE)

# # 提取物种名称 (格式如 name 列所示，剔除可能存在的分类层级前缀)
# bac_raw[[1]] <- sub(".*s__", "", bac_raw[[1]])

# 转换为以 Sample 为行、Feature (Species) 为列的数据框
df_data <- bac_raw %>%
  column_to_rownames(var = colnames(bac_raw)[1]) %>%
  t() %>%
  as.data.frame()

# 确保样本 ID 对齐
common_samples <- intersect(rownames(metadata), rownames(df_data))
metadata <- metadata[common_samples, , drop = FALSE]
df_data  <- df_data[common_samples, , drop = FALSE]

# ------------------------------------------------------------------------------
# 2. Run MaAsLin2 Analysis
# ------------------------------------------------------------------------------
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

Maaslin2(
  input_data     = df_data,
  input_metadata = metadata,
  output         = file.path(MAASLIN_DIR, "bac"),
  min_abundance  = MAASLIN_MIN_ABUNDANCE,
  min_prevalence = MAASLIN_MIN_PREVALENCE,
  normalization  = MAASLIN_NORMALIZATION,
  fixed_effects  = MAASLIN_FIXED_EFFECTS,
  reference      = MAASLIN_REFERENCE,
  plot_heatmap   = FALSE,
  plot_scatter   = FALSE
)
