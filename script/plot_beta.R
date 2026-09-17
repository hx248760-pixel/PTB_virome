library(vegan)
library(ggplot2)
library(patchwork)
library(dplyr)
library(ggpubr)

plot_pcoa_pro <- function(abd_dt, sample_dt,
                          colors = c("#56B4E9", "#D55E00", "#009E73", "#E69F00"),
                          comparisons = NULL,
                          xlim = NULL,
                          ylim = NULL,
                          adonis_out_name = NULL,
                          table_out_path = NULL,  # 导出 PERMANOVA 表格的文件路径
                          fig_out_path = NULL,    # 保存图片的完整文件路径
                          fig_width = 7,          # 保存图片的宽度
                          fig_height = 6          # 保存图片的高度
) {
  
  # 1. 预处理
  common_samples <- intersect(colnames(abd_dt), sample_dt$Sample)
  if (length(common_samples) == 0) {
    stop("错误: abd_dt 的列名与 sample_dt 的 Sample 列没有交集，请检查样本 ID。")
  }
  plot_df <- as.data.frame(t(abd_dt[, common_samples]))
  
  # 2. 动态检测分组逻辑
  if ("cohort" %in% colnames(sample_dt)) {
    message("检测到 'cohort' 列，将按 'cohort.Group' 进行复合分组。")
    sample_dt_processed <- sample_dt %>%
      mutate(Combined_Group = factor(paste(cohort, Group, sep = ".")))
  } else {
    message("未检测到 'cohort' 列，将直接按 'Group' 进行分组。")
    if (!"Group" %in% colnames(sample_dt)) {
      stop("错误: sample_dt 中缺少必需的 'Group' 列。")
    }
    sample_dt_processed <- sample_dt %>%
      mutate(Combined_Group = factor(Group))
  }
  
  # 3. 计算 PCoA
  dist_matrix <- vegdist(plot_df, method = "bray")
  pcoa_res <- cmdscale(dist_matrix, k = 2, eig = TRUE)
  var_exp <- round(pcoa_res$eig / sum(pcoa_res$eig) * 100, 2)
  
  # 4. 计算 PERMANOVA (提取 F 值, R2 和 p 值)
  meta_stat <- data.frame(Sample = rownames(plot_df)) %>%
    left_join(sample_dt_processed, by = "Sample")
  
  set.seed(123)
  perm_res <- adonis2(dist_matrix ~ Combined_Group, data = meta_stat, permutations = 999)
  
  # 提取统计指标
  f_val <- round(perm_res$F[1], 3)
  r2_val <- round(perm_res$R2[1], 3)
  p_val <- perm_res$`Pr(>F)`[1]
  
  p_text <- ifelse(p_val < 0.001, "P < 0.001", paste0("P = ", round(p_val, 3)))
  stat_label <- paste0("PERMANOVA\nF = ", f_val, "\nR² = ", r2_val, "\n", p_text)
  
  # 导出 adonis 统计表格到指定路径
  if (!is.null(table_out_path)) {
    perm_df <- as.data.frame(perm_res)
    perm_df <- cbind(Term = rownames(perm_df), perm_df)
    
    dir.create(dirname(table_out_path), recursive = TRUE, showWarnings = FALSE)
    write.table(perm_df, file = table_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("adonis2 结果表格已成功保存至: ", table_out_path)
  }
  
  # 赋值到环境变量（可选）
  if (!is.null(adonis_out_name)) {
    assign(adonis_out_name, perm_res, envir = parent.frame())
    message("adonis2 结果已赋值给变量 '", adonis_out_name, "'（位于调用环境）")
  }
  
  # 5. 整合绘图数据
  plot_data <- data.frame(
    Sample = rownames(pcoa_res$points),
    PCoA1 = pcoa_res$points[,1],
    PCoA2 = pcoa_res$points[,2]
  ) %>%
    left_join(sample_dt_processed, by = "Sample")
  
  # 6. 主图
  p_main <- ggplot(plot_data, aes(x = PCoA1, y = PCoA2, color = Combined_Group)) +
    geom_vline(xintercept = 0, color = "grey70", linetype = "dashed") +
    geom_hline(yintercept = 0, color = "grey70", linetype = "dashed") +
    geom_point(size = 3, alpha = 0.8) +
    stat_ellipse(aes(fill = Combined_Group), geom = "polygon", alpha = 0.1, show.legend = FALSE) +
    scale_color_manual(values = colors) +
    scale_fill_manual(values = colors) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    labs(
      x = paste0("PCoA1 (", var_exp[1], "%)"),
      y = paste0("PCoA2 (", var_exp[2], "%)"),
      tag = stat_label
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      legend.title = element_blank(),
      plot.tag.position = c(0.02, 0.98),
      plot.tag = element_text(size = 9, face = "plain", hjust = 0, vjust = 1, color = "black")
    )
  
  # 7. X轴边际箱线图
  p_x <- ggplot(plot_data, aes(x = Combined_Group, y = PCoA1, fill = Combined_Group)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons, label = "p.signif") +
    scale_fill_manual(values = colors) +
    coord_flip(ylim = xlim) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 8. Y轴边际箱线图
  p_y <- ggplot(plot_data, aes(x = Combined_Group, y = PCoA2, fill = Combined_Group)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(comparisons = comparisons, label = "p.signif") +
    scale_fill_manual(values = colors) +
    coord_cartesian(ylim = ylim) +
    theme_bw() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          legend.position = "none")
  
  # 9. 组合图像
  layout <- "
    A#
    BC
  "
  final <- p_x + p_main + p_y +
    plot_layout(design = layout, widths = c(4, 1), heights = c(1, 4))
  
  # 自动保存图片到 FIGURE_DIR
  if (!is.null(fig_out_path)) {
    dir.create(dirname(fig_out_path), recursive = TRUE, showWarnings = FALSE)
    ggsave(fig_out_path, plot = final, width = fig_width, height = fig_height)
    message("PCoA Beta多样性图已成功保存至: ", fig_out_path)
  }
  
  return(final)
}