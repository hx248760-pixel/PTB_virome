# ==============================================================================
# PTB Virome Project
# 02_virome_analysis.R
# Alpha diversity, beta diversity, dispersion testing, and viral family composition.
# ==============================================================================

rm(list = ls())
source("script/config.R")

library(tidyverse)
library(ggpubr)
library(vegan)
library(gridExtra)

# 加载自定义工具脚本
source("script/calculate_alpha.R")
source("script/plot_aplha.R")
source("script/plot_beta.R")
source("script/plot_family_barplot.R")

# ------------------------------------------------------------------------------
# 1. 读取输入文件 (Data Loading)
# ------------------------------------------------------------------------------
metadata     <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu         <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)
virus_family <- read.delim(VIRUS_FAMILY_FILE, check.names = FALSE, row.names = 1)

# 校验分组列
stopifnot("Group" %in% colnames(metadata))
stopifnot(all(metadata$Group %in% GROUP_LEVELS))
custom_colors <- GROUP_COLORS

# ------------------------------------------------------------------------------
# 2. Alpha 多样性计算与绘图 (Alpha Diversity)
# ------------------------------------------------------------------------------
set.seed(123)

# 修正 1：显式接收返回的数据框 alpha_re，彻底避免 object not found 报错
alpha_re <- calculate_alpha_diversity(
  abundance_data = votu,
  sample_col     = "Sample",
  group_col      = "Group",
  group_data     = metadata
)

# 导出 Alpha 多样性指标表格供审查
write.table(alpha_re, file.path(TABLE_DIR, "alpha_diversity_metrics.tsv"), 
            sep = "\t", quote = FALSE, row.names = FALSE)

# 绘制与保存 Alpha 多样性箱线图
p_simpson <- ggplot(alpha_re, aes(x = Group, y = Simpson, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = custom_colors) +
  theme_bw() + theme(legend.position = "none") + labs(x = "", y = "Simpson Index")

p_shannon <- ggplot(alpha_re, aes(x = Group, y = Shannon, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = custom_colors) +
  theme_bw() + theme(legend.position = "none") + labs(x = "", y = "Shannon Index")

p_obs <- ggplot(alpha_re, aes(x = Group, y = Richness, fill = Group)) +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.4, size = 1.5) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = custom_colors) +
  theme_bw() + theme(legend.position = "none") + labs(x = "", y = "Observed vOTU Richness")

ggsave(file.path(FIGURE_DIR, "otus_simpson.pdf"), p_simpson, width = 4, height = 3.5)
ggsave(file.path(FIGURE_DIR, "otus_shannon.pdf"), p_shannon, width = 4, height = 3.5)
ggsave(file.path(FIGURE_DIR, "otus_richness.pdf"), p_obs, width = 4, height = 3.5)

# ------------------------------------------------------------------------------
# 3. Beta 多样性 PCoA (Beta Diversity)
# ------------------------------------------------------------------------------
plot_pcoa_pro(
  abd_dt          = votu,
  sample_dt       = metadata,
  colors          = custom_colors,
  adonis_out_name = file.path(TABLE_DIR, "adonis_votu")
)

plot_pcoa_pro(
  abd_dt          = virus_family,
  sample_dt       = metadata,
  colors          = custom_colors,
  adonis_out_name = file.path(TABLE_DIR, "adonis_fa")
)

# ------------------------------------------------------------------------------
# 4. Bray-Curtis 离散度齐性检验 (Betadisper Test)
# ------------------------------------------------------------------------------
# 修正 2：改用更加稳健的维度对齐方式 (确保行为 Sample，列为 vOTUs)
if (ncol(votu) == nrow(metadata)) {
  abd_matrix <- t(votu)
} else {
  abd_matrix <- as.matrix(votu)
}

# 按照 Sample 精确对齐 metadata
common_samples <- intersect(rownames(abd_matrix), metadata$Sample)
abd_matrix_aligned <- abd_matrix[common_samples, ]
metadata_aligned   <- metadata[match(common_samples, metadata$Sample), ]

dist_matrix <- vegdist(abd_matrix_aligned, method = "bray")
mod         <- betadisper(dist_matrix, metadata_aligned$Group)
disp_test   <- permutest(mod, permutations = BETADISPER_PERMUTATIONS)

# 输出组间离散度检验结果表
write.table(as.data.frame(disp_test$tab), file.path(TABLE_DIR, "betadisper_votu.tsv"),
            sep = "\t", quote = FALSE)
print(disp_test)

# ------------------------------------------------------------------------------
# 5. 病毒科水平组成堆叠图 (Viral Family Barplot)
# ------------------------------------------------------------------------------
# 修正 3：复用开头读取好的 virus_family，无需重复读盘
plot_family_abundance(
  abundance_df    = virus_family,
  meta_df         = metadata,
  group_col       = "Group",
  bar_colors      = c("#D44D5C", "#33A02C", "#1F78B4", "#FF7F00",
                      "#6A3D9A", "#B15928", "#E31A1C", "#A6CEE3",
                      "#B2DF8A", "#FDBF6F", "#CAB2D6", "#FFFF99", "grey"),
  top_n           = 12,
  y_range         = c(0, 50),
  cluster_samples = TRUE,
  out_pdf         = file.path(FIGURE_DIR, "family.pdf")
)
