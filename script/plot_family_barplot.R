# --------------------------------------------------------------
# 函数名称：plot_family_abundance_df
# 说明：从已有 data.frame 开始执行全部绘图流程
# --------------------------------------------------------------
plot_family_abundance <- function(abundance_df,   # 行 = Family，列 = Sample (相对丰度 %)
                                  meta_df,        # 至少包含列 Sample 与分组信息
                                  group_col = "Group",   # 用于 facet 的分组列名
                                  bar_colors = NULL,     # 前 top_n 种颜色 + Others（灰色），
                                  top_n = 10,
                                  y_range = NULL,
                                  cluster_samples = TRUE,# 是否按 Bray‑Curtis 聚类排序样本
                                  out_pdf = "family_abundance.pdf" ) {
  ## ----------------- 依赖包 -----------------
  if (!requireNamespace("tidyverse", quietly = TRUE)) {
    stop("Package 'tidyverse' is required but not installed.")
  }
  if (!requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required but not installed.")
  }
  library(tidyverse)
  library(vegan)
  
  ## ----------------- Step 3：宽表 → 长表 -----------------
  long_df <- abundance_df %>%
    as.data.frame() %>%
    rownames_to_column(var = "Family") %>%
    pivot_longer(cols = -Family,
                 names_to = "Sample",
                 values_to = "Abundance")
  
  ## ----------------- Step 4：挑选前 top_n 大科 -----------------
  top_families <- long_df %>%
    group_by(Family) %>%
    summarise(Total = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(Total)) %>%
    slice_head(n = top_n) %>%
    pull(Family)
  
  ## ----------------- Step 5：合并 “Others” -----------------
  df2 <- long_df %>%
    mutate(Family = ifelse(Family %in% top_families, Family, "Others")) %>%
    group_by(Sample, Family) %>%
    summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    mutate(Family = factor(Family,
                           levels = c(top_families, "Others")))
  
  ## ----------------- Step 6：样本聚类排序（可选） -----------------
  if (cluster_samples) {
    # 宽表（样本 × 科）用于计算 Bray‑Curtis 距离
    wide_mat <- df2 %>%
      pivot_wider(names_from = Family,
                  values_from = Abundance,
                  values_fill = 0) %>%
      column_to_rownames(var = "Sample")
    
    bc_dist <- vegdist(wide_mat, method = "bray")
    hc <- hclust(bc_dist, method = "average")
    sample_order <- rownames(wide_mat)[hc$order]
  } else {
    sample_order <- unique(df2$Sample)
  }
  
  ## ----------------- Step 7：合并分组信息 -----------------
  plot_df <- df2 %>%
    left_join(meta_df, by = c("Sample" = "Sample"))
  
  ## ----------------- Step 8：配色 -----------------
  if (is.null(bar_colors)) {
    # 默认配色：前 top_n 种配色 + Others 灰色
    default_pal <- scales::hue_pal()(top_n)
    bar_colors <- c(default_pal, "gray")
  }
  # 若颜色数量不足则循环复用
  if (length(bar_colors) < (top_n + 1)) {
    warning("Provided bar_colors fewer than needed; colors will be recycled.")
    bar_colors <- rep(bar_colors, length.out = top_n + 1)
  }
  # 按 Family 因子顺序命名颜色向量
  names(bar_colors) <- levels(plot_df$Family)
  
  ## ----------------- Step 9：绘图 -----------------
  p <- ggplot(plot_df,
              aes(x = factor(Sample, levels = sample_order),
                  y = Abundance,
                  fill = Family)) +
    geom_bar(stat = "identity", width = 0.8) +
    scale_fill_manual(values = bar_colors) +
    labs(x = "Sample",
         y = "Relative Abundance (%)",
         fill = "Family") +
    facet_grid(as.formula(paste0("~ ", group_col)),
               scales = "free_x", space = "free_x") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5),
          legend.position = "right")
  
  if (!is.null(y_range)){p <- p + coord_cartesian(ylim=y_range)}
  
  ## ----------------- Step 10：导出 PDF -----------------
  ggsave(filename = out_pdf,
         plot = p,
         device = "pdf",
         width = 10, height = 6, units = "in")
  
  ## 返回 ggplot 对象（便于继续编辑）
  invisible(p)
}
