# ==============================================================================
# PTB Virome Project
# 05_prepare_inputf.R
# Prepare input files for correlation analysis, run calc_corr.py, and filter results.
# ==============================================================================

rm(list = ls())
source("script/config.R")

library(dplyr)
library(tidyr)

# ------------------------------------------------------------------------------
# 1. 读取 MaAsLin2 显著结果和丰度表
# ------------------------------------------------------------------------------
# 读取显著结果
sig_bac <- read.delim(file.path(MAASLIN_DIR, "bac", "significant_results.tsv"), 
                      check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(metadata == "Group")

sig_vir <- read.delim(file.path(MAASLIN_DIR, "vir", "significant_results.tsv"), 
                      check.names = FALSE, stringsAsFactors = FALSE) %>%
  filter(metadata == "Group")

# 读取丰度表 (行名为物种/vOTU，列名为样本)
vir_spe <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)
bac_spe <- read.delim(BACTERIAL_FILE, check.names = FALSE, row.names = 1)

# 筛选显著特征
sig_vir_spe <- vir_spe[rownames(vir_spe) %in% sig_vir$feature, , drop = FALSE]
sig_bac_spe <- bac_spe[rownames(bac_spe) %in% sig_bac$feature, , drop = FALSE]

# ------------------------------------------------------------------------------
# 2. 按组划分样本并筛选出现率 > OCCURRENCE_THRESHOLD 的物种
# ------------------------------------------------------------------------------
metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
rownames(metadata) <- metadata$Sample

# 根据 config 中的分组名称：Health 和 Patient
health_samples <- rownames(metadata)[metadata$Group == "Health"]
patient_samples <- rownames(metadata)[metadata$Group == "Patient"]

# 提取各组样本的丰度子集
hc_sig_vir <- sig_vir_spe[, colnames(sig_vir_spe) %in% health_samples, drop = FALSE]
hc_sig_bac <- sig_bac_spe[, colnames(sig_bac_spe) %in% health_samples, drop = FALSE]

tb_sig_vir <- sig_vir_spe[, colnames(sig_vir_spe) %in% patient_samples, drop = FALSE]
tb_sig_bac <- sig_bac_spe[, colnames(sig_bac_spe) %in% patient_samples, drop = FALSE]

# 出现率筛选函数
filter_by_prevalence <- function(df, threshold = OCCURRENCE_THRESHOLD) {
  if (ncol(df) == 0) return(df)
  prevalence <- rowSums(df > 0) / ncol(df)
  return(df[prevalence > threshold, , drop = FALSE])
}

hc_selected_vir <- filter_by_prevalence(hc_sig_vir)
hc_selected_bac <- filter_by_prevalence(hc_sig_bac)
tb_selected_vir <- filter_by_prevalence(tb_sig_vir)
tb_selected_bac <- filter_by_prevalence(tb_sig_bac)

# 保存为 calc_corr.py 的输入文件
corr_input_dir <- file.path(TABLE_DIR, "correlation_input")
dir.create(corr_input_dir, recursive = TRUE, showWarnings = FALSE)

write.table(hc_selected_vir, file.path(corr_input_dir, "hc_sig_vir.txt"), 
            quote = FALSE, sep = "\t")
write.table(hc_selected_bac, file.path(corr_input_dir, "hc_sig_bac.txt"), 
            quote = FALSE, sep = "\t")
write.table(tb_selected_vir, file.path(corr_input_dir, "tb_sig_vir.txt"), 
            quote = FALSE, sep = "\t")
write.table(tb_selected_bac, file.path(corr_input_dir, "tb_sig_bac.txt"), 
            quote = FALSE, sep = "\t")

# ------------------------------------------------------------------------------
# 3. 调用 calc_corr.py 计算相关性
# ------------------------------------------------------------------------------
# calc_corr.py 位于 script/ 目录下
calc_corr_script <- file.path("script", "calc_corr.py")

# 定义输出文件
hc_re_out <- file.path(corr_input_dir, "hc_re.txt")
tb_re_out <- file.path(corr_input_dir, "tb_re.txt")

# 构建命令 (使用 python3 调用)
python_cmd <- "python3"
args_hc <- c(calc_corr_script, "-i", file.path(corr_input_dir, "hc_sig_vir.txt"),
             "-I", file.path(corr_input_dir, "hc_sig_bac.txt"),
             "-o", hc_re_out, "-s", "0")
args_tb <- c(calc_corr_script, "-i", file.path(corr_input_dir, "tb_sig_vir.txt"),
             "-I", file.path(corr_input_dir, "tb_sig_bac.txt"),
             "-o", tb_re_out, "-s", "0")

message("Running correlation analysis for HC...")
system2(python_cmd, args = args_hc)
message("Running correlation analysis for TB...")
system2(python_cmd, args = args_tb)

# ------------------------------------------------------------------------------
# 4. 读取相关性结果，计算 FDR，过滤并保存
# ------------------------------------------------------------------------------
re_hc <- read.delim(hc_re_out, check.names = FALSE, stringsAsFactors = FALSE)
re_hc$fdr <- p.adjust(re_hc$pvalue, method = "fdr")

re_tb <- read.delim(tb_re_out, check.names = FALSE, stringsAsFactors = FALSE)
re_tb$fdr <- p.adjust(re_tb$pvalue, method = "fdr")

re_hc_filter <- re_hc %>% 
  filter(abs(corr) >= CORRELATION_THRESHOLD & fdr <= FDR_THRESHOLD)

re_tb_filter <- re_tb %>% 
  filter(abs(corr) >= CORRELATION_THRESHOLD & fdr <= FDR_THRESHOLD)

# 保存到 TABLE_DIR，供 05_network_analysis.R 使用
write.table(re_hc_filter, file.path(TABLE_DIR, "hc_sig.re"), 
            quote = FALSE, sep = "\t", row.names = FALSE)
write.table(re_tb_filter, file.path(TABLE_DIR, "tb_sig.re"), 
            quote = FALSE, sep = "\t", row.names = FALSE)

message("Correlation analysis completed. Results saved to ", TABLE_DIR)
