# ==============================================================================
# PTB Virome Project
# 03_differential_analysis.R
# Differential viral/family features, volcano plot, family pie charts, and viral-host plots.
# ==============================================================================

rm(list = ls())
source("script/config.R")

library(tidyverse)
library(Maaslin2) # 修正 1：规范包名大小写
library(ggplot2)

dir.create(MAASLIN_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. 读取输入文件 (Data Loading)
# ------------------------------------------------------------------------------
metadata     <- read.delim(METADATA_FILE, check.names = FALSE, stringsAsFactors = FALSE)
votu         <- read.delim(VOTU_FILE, check.names = FALSE, row.names = 1)
virus_family <- read.delim(VIRUS_FAMILY_FILE, check.names = FALSE, row.names = 1)

stopifnot(all(metadata$Group %in% GROUP_LEVELS))
rownames(metadata) <- metadata$Sample

# ------------------------------------------------------------------------------
# 2. MaAsLin2 vOTU 级别差异分析 (vOTU Differential Abundance)
# ------------------------------------------------------------------------------
Maaslin2(
  input_data     = votu,
  input_metadata = metadata,
  output         = file.path(MAASLIN_DIR, "vir"),
  min_abundance  = MAASLIN_MIN_ABUNDANCE,
  min_prevalence = MAASLIN_MIN_PREVALENCE,
  normalization  = MAASLIN_NORMALIZATION,
  fixed_effects  = MAASLIN_FIXED_EFFECTS,
  reference      = MAASLIN_REFERENCE,
  plot_heatmap   = FALSE,
  plot_scatter   = FALSE
)

# 提取结果并增加 95% CI 计算
vir_results <- read.delim(file.path(MAASLIN_DIR, "vir", "significant_results.tsv"), check.names = FALSE) %>%
  filter(metadata == "Group") %>%
  mutate(
    CI_95_Low  = coef - 1.96 * stderr,
    CI_95_High = coef + 1.96 * stderr,
    CI_95      = sprintf("[%.4f, %.4f]", CI_95_Low, CI_95_High)
  )

write.table(vir_results, file.path(TABLE_DIR, "sig_vir_otu.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# ------------------------------------------------------------------------------
# 3. MaAsLin2 病毒科级别差异分析 (Family Differential Abundance)
# ------------------------------------------------------------------------------
Maaslin2(
  input_data     = virus_family,
  input_metadata = metadata,
  output         = file.path(MAASLIN_DIR, "vir_family"),
  min_abundance  = MAASLIN_MIN_ABUNDANCE,
  min_prevalence = MAASLIN_MIN_PREVALENCE,
  normalization  = MAASLIN_NORMALIZATION,
  fixed_effects  = MAASLIN_FIXED_EFFECTS,
  reference      = MAASLIN_REFERENCE,
  plot_heatmap   = FALSE,
  plot_scatter   = FALSE
)

family_results <- read.delim(file.path(MAASLIN_DIR, "vir_family", "significant_results.tsv"), check.names = FALSE) %>%
  filter(metadata == "Group") %>%
  mutate(
    CI_95_Low  = coef - 1.96 * stderr,
    CI_95_High = coef + 1.96 * stderr,
    CI_95      = sprintf("[%.4f, %.4f]", CI_95_Low, CI_95_High)
  )

write.table(family_results, file.path(TABLE_DIR, "sig_vir_family.txt"), quote = FALSE, row.names = FALSE, sep = "\t")

# ------------------------------------------------------------------------------
# 4. 火山图绘制 (Volcano Plot: exploratory q < 0.25 threshold)
# ------------------------------------------------------------------------------
# 注意：MAASLIN_REFERENCE 为 "Group,Patient"，因此 coef > 0 代表 Health 富集
votu_all <- read.delim(file.path(MAASLIN_DIR, "vir", "all_results.tsv"), check.names = FALSE) %>%
  filter(metadata == "Group")

re_ma <- votu_all %>%
  mutate(Group = case_when(
    coef > 0 & qval < MAASLIN_Q_THRESHOLD ~ "Health",
    coef < 0 & qval < MAASLIN_Q_THRESHOLD ~ "Patient",
    TRUE ~ "Not significant"
  ))

p_volcano <- ggplot(re_ma, aes(x = coef, y = -log10(qval), color = Group)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c(
    "Health"          = "#33a02c",
    "Patient"         = "#D44D5C",
    "Not significant" = "#AAAAAA"
  )) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = -log10(MAASLIN_Q_THRESHOLD), linetype = "dashed", color = "darkgray") +
  labs(x = "Coefficient (coef)", y = "-log10(q value)",
       title = "Volcano Plot of Features", color = "Group") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        axis.title = element_text(size = 12), 
        legend.title = element_text(size = 10),
        legend.position = "bottom")

ggsave(file.path(FIGURE_DIR, "votu_volcano_q0.25.pdf"), p_volcano, width = 7, height = 5)

# ------------------------------------------------------------------------------
# 5. 显著 vOTU 病毒科饼图绘制 (Pie Charts)
# ------------------------------------------------------------------------------
family_color_map <- c(
  "unknow"           = "#D44D5C", 
  "Aliceevansviridae"= "#af5727",
  "Herelleviridae"   = "#fe7d00", 
  "Inoviridae"       = "#b1de88",
  "Fervensviridae"   = "#D56AA0", 
  "Microviridae"     = "#329e2b",
  "Winoviridae"      = "#1e77b3"
)

draw_pie_chart <- function(data, title_text, output_file, color_map) {
  if (nrow(data) == 0) {
    warning("No features available for: ", title_text)
    return(invisible(NULL))
  }
  plot_df <- data %>%
    count(Family, name = "n") %>%
    mutate(Family = ifelse(is.na(Family) | Family == "", "unknow", Family))
  
  p <- ggplot(plot_df, aes(x = "", y = n, fill = Family)) +
    geom_col(width = 1, color = "white") +
    coord_polar(theta = "y") +
    scale_fill_manual(values = color_map, na.value = "#AAAAAA") +
    labs(title = title_text, fill = "Virus family") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  
  ggsave(output_file, p, width = 5, height = 5)
  invisible(p)
}

sig_otu <- read.delim(file.path(MAASLIN_DIR, "vir", "significant_results.tsv"), check.names = FALSE) %>%
  filter(metadata == "Group")

vir_tax <- read.delim(VOTU_FAMILY_ANNOTATION, header = F, check.names = FALSE)
sig_otu$Family <- vir_tax$V4[match(sig_otu$feature, vir_tax$V1)]
sig_otu$Family[is.na(sig_otu$Family)] <- "unknow"
sig_otu$Group <- case_when(sig_otu$coef < 0 ~ "Patient", sig_otu$coef > 0 ~ "Health")

Health_data  <- sig_otu %>% filter(Group == "Health")
Patient_data <- sig_otu %>% filter(Group == "Patient")

draw_pie_chart(Patient_data,
               paste0("TB-enriched (n=", dplyr::n_distinct(Patient_data$feature), ")"),
               file.path(FIGURE_DIR, "tb_pie_0.25.pdf"), family_color_map)

draw_pie_chart(Health_data,
               paste0("HC-enriched (n=", dplyr::n_distinct(Health_data$feature), ")"),
               file.path(FIGURE_DIR, "hc_pie_0.25.pdf"), family_color_map)

# ------------------------------------------------------------------------------
# 6. 病毒-宿主 (Viral-Host) 对应关联分析与导出
# ------------------------------------------------------------------------------
vir_host <- read.delim(VOTU_HOST_ANNOTATION, header = F, check.names = FALSE)

vir_sig_Patient_host <- vir_host[vir_host$V1 %in% Patient_data$feature, , drop = FALSE]
vir_sig_Health_host  <- vir_host[vir_host$V1 %in% Health_data$feature, , drop = FALSE]

colnames(vir_sig_Patient_host) <- c("V1", "host_num", "bac_genus")
colnames(vir_sig_Health_host)  <- c("V1", "host_num", "bac_genus")

vir_sig_Patient_host <- vir_sig_Patient_host %>%
  left_join(vir_tax %>% select(V1, V4) %>% rename(virus_family = V4), by = "V1")
vir_sig_Health_host  <- vir_sig_Health_host %>%
  left_join(vir_tax %>% select(V1, V4) %>% rename(virus_family = V4), by = "V1")

vir_sig_Patient_host$virus_family[is.na(vir_sig_Patient_host$virus_family)] <- "unclassified"
vir_sig_Health_host$virus_family[is.na(vir_sig_Health_host$virus_family)]   <- "unclassified"

write.table(vir_sig_Health_host, file.path(TABLE_DIR, "Health_sig_vir_host.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")
write.table(vir_sig_Patient_host, file.path(TABLE_DIR, "Patient_sig_vir_host.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")

#------宿主预测绘图
source("script/plot_virf_host.R")

my_manual_colors <- c(
  "Bacteroides" = "#DD9AC2", "Agathobaculum" = "#dda03c",
  "UBA9502" = "#2DD881", "Acetatifactor" = "#7FB7BE",
  "Agathobacter" = "#C874D9", "Angelakisella" = "#FFC8FB",
  "Clostridium_Q" = "#59C3C3", "Roseburia" = "#84E6F8",
  "Lachnospira" = "#52489C", "Other" = "#B6C9BB",
  "Gemmiger" = "#E2F1AF", "multiple" = "#7ADFBB",
  "Blautia_A" = "#CC3363", "Clostridium_M" = "#B9FFB7",
  "Parabacteroides" = "#EF946C"
)

plot_genus_stacked(
  data = vir_sig_Patient_host, min_occurrence = 6,
  manual_colors = my_manual_colors,
  plot_title = "TB-enriched vOTUs",
  output_file = file.path(FIGURE_DIR, "TB_vir_host.pdf"), width = 6, height = 5
)

plot_genus_stacked(
  data = vir_sig_Health_host, min_occurrence = 6,
  manual_colors = my_manual_colors,
  plot_title = "HC-enriched vOTUs",
  output_file = file.path(FIGURE_DIR, "Health_vir_host.pdf"), width = 6, height = 5
)
