# ==============================================================================
# PTB Virome Project
# 04_functional_analysis.R
# Functional enrichment of KO annotations (Fisher's exact test) & Visualization.
# ==============================================================================

rm(list = ls())

source("script/config.R")

library(tidyverse)
library(ggplot2)
library(ggrepel)
library(scales)

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
metadata <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu     <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)

# Use the MaAsLin2 significant vOTUs from the differential analysis
co_re <- read.delim("results/maaslin/vir/significant_results.tsv") %>%
  filter(metadata == "Group") %>%
  mutate(Group = ifelse(coef > 0, "HC", "TB"))

mapping_votu_ko <- read.delim(VOTU_KO_MAPPING_FILE,header = F,
                              check.names = FALSE, stringsAsFactors = FALSE)

names(mapping_votu_ko)[1:2] <- c("vOTU", "KO")

combined_data <- co_re %>%
  select(feature, Group) %>%
  inner_join(mapping_votu_ko, by = c("feature" = "vOTU")) %>%
  distinct(feature, Group, KO)

# === 预设 vOTU 总数（作为分母） ===
total_votu <- co_re %>%
  distinct(feature, Group) %>%
  count(Group, name = "total_votu")

n_hc <- total_votu$total_votu[match("HC", total_votu$Group)]
n_tb <- total_votu$total_votu[match("TB", total_votu$Group)]

# === 统计每个 KO 在两个组中出现的 vOTU 数量 ===
ko_counts <- combined_data %>%
  group_by(Group, KO) %>%
  summarise(count = n_distinct(feature), .groups = "drop")

# === 构建完整频数表，补全 0 ===
merged_table <- ko_counts %>%
  pivot_wider(names_from = Group, values_from = count, values_fill = list(count = 0))

if (!"HC" %in% colnames(merged_table)) merged_table$HC <- 0L
if (!"TB" %in% colnames(merged_table)) merged_table$TB <- 0L

merged_table <- merged_table %>%
  rename(count_in_HC = HC, count_in_TB = TB) %>%
  mutate(
    count_not_in_HC = n_hc - count_in_HC,
    count_not_in_TB = n_tb - count_in_TB
  )

# ------------------------------------------------------------------------------
# 2. Fisher Exact Test
# ------------------------------------------------------------------------------
fisher_test <- function(a, b, c, d) {
  mat <- matrix(c(a, b, c, d), nrow = 2)
  fisher.test(mat)$p.value
}

merged_table <- merged_table %>%
  mutate(
    fisher_p_value = mapply(
      fisher_test,
      count_in_HC, count_in_TB,
      count_not_in_HC, count_not_in_TB
    ),
    
    # 出现率计算
    occurrence_rate_in_HC = count_in_HC / n_hc,
    occurrence_rate_in_TB = count_in_TB / n_tb,
    
    # 差值 & 富集方向
    diff = occurrence_rate_in_TB - occurrence_rate_in_HC,
    enriched_in = ifelse(diff > 0, "TB", "HC"),
    
    # 添加显著性标记
    signif_label = case_when(
      fisher_p_value < 0.001 ~ "***",
      fisher_p_value < 0.01  ~ "**",
      fisher_p_value < 0.05  ~ "*",
      TRUE ~ ""
    )
  )

# === 合并注释 ===
# 读取 KO 注释文件 (无表头，自动命名为 V1, V2, ...)
ko_description <- read.delim(
  KO_DESCRIPTION_FILE,
  header = FALSE,
  sep = "\t",
  quote = "",
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  # V5 为 KO 号 (如 K00844)，V6 为基因功能描述 (如 HK; hexokinase...)
  select(KO = V7, Description = V8) %>%
  distinct(KO, .keep_all = TRUE) # 一个 KO 只保留一条描述，防止多对多关联膨胀

# 合并注释并筛选显著结果
final_table <- merged_table %>%
  left_join(ko_description, by = "KO") %>%
  filter(!is.na(Description)) %>% # 过滤掉未匹配到描述的 KO
  mutate(Gene_Description = paste(KO, Description, sep = ": "))%>%
  filter(fisher_p_value<=0.05)

write.csv(final_table, file = "results/tables/KO_enriched_results.csv", row.names = FALSE, quote = FALSE)

# ------------------------------------------------------------------------------
# 3. Plot 1: Occurrence Rate (原版风格：黄条带 + 显著性星号)
# ------------------------------------------------------------------------------

diff_data <- final_table %>% arrange(desc(diff))
sorted_genes <- unique(diff_data$Gene_Description)

# 转为长格式
long_data <- pivot_longer(
  diff_data,
  cols = c(occurrence_rate_in_HC, occurrence_rate_in_TB),
  names_to = "group",
  values_to = "occurrence_rate"
)

# 设置 KO 排序
long_data$Gene_Description <- factor(long_data$Gene_Description, levels = sorted_genes)

custom_colors <- c(
  "occurrence_rate_in_HC" = "#2e8153", 
  "occurrence_rate_in_TB" = "#bc5858"
)

num_levels <- length(levels(long_data$Gene_Description))
x_range <- c(0, max(long_data$occurrence_rate, na.rm = TRUE) * 1.05)

# 背景色条带（偶数行）
bg_rects <- data.frame(
  ymin = seq(1, num_levels, by = 2),
  ymax = seq(2, num_levels + 1, by = 2)
)

# 星号位置数据（标在 TB 的点右侧）
star_data <- diff_data %>%
  mutate(
    y_pos = seq_along(Gene_Description) + 0.5,
    x_pos = occurrence_rate_in_TB + 0.01  # 星号稍偏右
  ) %>%
  filter(signif_label != "")

# 画图
p_occurrence <- ggplot(long_data, aes(x = occurrence_rate, y = as.numeric(Gene_Description) + 0.5, color = group)) +
  # 背景色条带
  geom_rect(
    data = bg_rects,
    aes(ymin = ymin, ymax = ymax),
    xmin = x_range[1], xmax = x_range[2],
    fill = "#FFFACD", alpha = 0.3,
    inherit.aes = FALSE
  ) +
  # 散点
  geom_point(size = 2) +
  # 添加显著性星号
  geom_text(
    data = star_data, 
    aes(x = x_pos, y = y_pos, label = signif_label),
    inherit.aes = FALSE, size = 4.5, hjust = 0
  ) +
  scale_x_continuous(labels = percent, limits = x_range) +
  scale_color_manual(values = custom_colors) +
  scale_y_continuous(
    breaks = seq_along(levels(long_data$Gene_Description)) + 0.5,
    labels = levels(long_data$Gene_Description),
    expand = expansion(add = c(0.5, 0.5))
  ) +
  labs(x = "Occurrence Rate", y = "Gene: Description", color = "Group") +
  theme_minimal(base_size = 12) +
  theme(
    axis.title       = element_text(size = 12, face = "bold"),
    axis.text        = element_text(size = 10),
    axis.line        = element_line(color = "black"),
    axis.ticks       = element_line(color = "black"),
    panel.grid       = element_blank(),
    legend.position  = "top",
    legend.direction = "horizontal"
  )

# 保存图像
ggsave("results/figures/kegg_occurrence.pdf", plot = p_occurrence, 
       width = 10, height = num_levels * 0.25 + 2, units = "in", dpi = 300)

# ------------------------------------------------------------------------------
# 4. Plot 2: Volcano Plot (原版风格)
# ------------------------------------------------------------------------------
volcano_data <- merged_table %>%
  mutate(
    log2FC = log2((occurrence_rate_in_TB + 1e-6) / (occurrence_rate_in_HC + 1e-6)),
    neg_log10_p = -log10(fisher_p_value),
    enrichment = case_when(
      fisher_p_value < 0.05 & log2FC > 0 ~ "Enriched in TB",
      fisher_p_value < 0.05 & log2FC < 0 ~ "Enriched in HC",
      TRUE ~ "Not Significant"
    ),
    significant = fisher_p_value < 0.05
  )

p_kegg_volcano <- ggplot(volcano_data, aes(x = log2FC, y = neg_log10_p)) +
  # 所有点，按 enrichment 上色
  geom_point(aes(color = enrichment), size = 1.5, alpha = 0.7) +
  
  # 显著点，高亮为中空圆圈，填充颜色也根据 enrichment 分类
  geom_point(
    data = subset(volcano_data, significant),
    aes(x = log2FC, y = neg_log10_p, fill = enrichment),
    shape = 21, color = "black", size = 3, stroke = 1
  ) +
  
  # 添加显著 KO 编号标签
  geom_text_repel(
    data = subset(volcano_data, significant),
    aes(label = KO),
    size = 3,
    box.padding = 0.3,
    max.overlaps = 100
  ) +
  
  # 自定义颜色
  scale_color_manual(values = c(
    "Enriched in TB"  = "#bc5858",
    "Enriched in HC"  = "#2e8153",
    "Not Significant" = "grey70"
  )) +
  scale_fill_manual(values = c(
    "Enriched in TB"  = "#bc5858",
    "Enriched in HC"  = "#2e8153",
    "Not Significant" = "grey70"
  )) +
  
  # 参考线
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  
  # 标题与标签
  labs(
    title = "Volcano Plot of Functional Enrichment with KO Labels",
    x = "log2(Fold Change)",
    y = "-log10(p-value)",
    color = "Enrichment",
    fill = "Enrichment"
  ) +
  
  # 美化主题
  theme_minimal(base_size = 13)

# 保存图像
ggsave("results/figures/kegg_volcano.pdf", plot = p_kegg_volcano, height = 5, width = 7)
