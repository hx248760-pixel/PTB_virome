# ==============================================================================
# PTB Virome Project
# plot_virf_host.R
# Plot stacked barplots of predicted viral host genus grouped by virus family.
# ==============================================================================

library(dplyr)
library(ggplot2)

plot_genus_stacked <- function(data, 
                               min_occurrence = 6, 
                               manual_colors = NULL, 
                               plot_title = "Host Genus Distribution", 
                               output_file = NULL, 
                               width = 6, 
                               height = 5) {
  
  if (is.null(data) || nrow(data) == 0) {
    warning("输入数据为空，跳过绘图: ", plot_title)
    return(invisible(NULL))
  }
  
  # 1. 基础归类：host_num > 1 -> "multiple"
  processed_step1 <- data %>%
    mutate(
      virus_family = ifelse(is.na(virus_family) | virus_family == "", "unclassified", virus_family),
      host_num     = as.numeric(host_num),
      raw_genus    = case_when(
        host_num > 1 ~ "multiple",
        is.na(bac_genus) | bac_genus == "" ~ "Other",
        TRUE ~ bac_genus
      )
    )
  
  # 2. 统计单宿主 (host_num == 1) 的频次，筛选高频属
  single_host_counts <- processed_step1 %>%
    filter(raw_genus != "multiple" & raw_genus != "Other") %>%
    count(raw_genus, name = "genus_freq")
  
  top_single_genera <- single_host_counts %>%
    filter(genus_freq >= min_occurrence) %>%
    pull(raw_genus)
  
  # 3. 低频单宿主属归类为 "Other"
  final_df <- processed_step1 %>%
    mutate(
      plot_genus = case_when(
        raw_genus == "multiple" ~ "multiple",
        raw_genus %in% top_single_genera ~ raw_genus,
        TRUE ~ "Other"
      )
    )
  
  # 4. 按 病毒科 (virus_family) 与 宿主分类 (plot_genus) 汇总 vOTU 数量
  plot_summary <- final_df %>%
    count(virus_family, plot_genus, name = "vOTU_count")
  
  genus_rank <- plot_summary %>%
    group_by(plot_genus) %>%
    summarise(total = sum(vOTU_count), .groups = "drop") %>%
    arrange(desc(total))
  
  # 转换为 Factor，完全按总数量大小排列 Levels
  final_levels <- genus_rank$plot_genus
  plot_summary$plot_genus <- factor(plot_summary$plot_genus, levels = final_levels)
  
  # 5. 颜色配置与映射
  unique_genera <- levels(plot_summary$plot_genus)
  
  if (!is.null(manual_colors)) {
    use_colors <- manual_colors[names(manual_colors) %in% unique_genera]
    missing_genera <- setdiff(unique_genera, names(use_colors))
    if (length(missing_genera) > 0) {
      extra_colors <- rep("#B6C9BB", length(missing_genera))
      names(extra_colors) <- missing_genera
      use_colors <- c(use_colors, extra_colors)
    }
  } else {
    use_colors <- scales::hue_pal()(length(unique_genera))
    names(use_colors) <- unique_genera
  }
  
  # 如果没有在 manual_colors 中指定，赋予默认指定色
  if ("Other" %in% names(use_colors) && is.null(manual_colors["Other"])) {
    use_colors["Other"] <- "#B6C9BB"
  }
  if ("multiple" %in% names(use_colors) && is.null(manual_colors["multiple"])) {
    use_colors["multiple"] <- "#7ADFBB"
  }

  # ----------------------------------------------------------------------------
  # 6. 绘图（不使用描边线 color = NA）
  # ----------------------------------------------------------------------------
  p <- ggplot(plot_summary, aes(x = virus_family, y = vOTU_count, fill = plot_genus)) +
    geom_col(width = 0.6, color = NA) +
    scale_fill_manual(values = use_colors) +
    labs(
      x = NULL,
      y = "Number of vOTUs",
      title = plot_title,
      fill = "Genus"
    ) +
    theme_bw() +
    theme(
      plot.title         = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.text.x        = element_text(angle = 45, hjust = 1, vjust = 1, size = 10, face = "italic", color = "black"),
      axis.text.y        = element_text(size = 10, color = "black"),
      axis.title.y       = element_text(size = 11, face = "bold"),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.title       = element_text(size = 10, face = "bold"),
      legend.text        = element_text(size = 9, face = "italic"),
      legend.position    = "right"
    )
  
  # 7. 保存图片
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ggsave(output_file, plot = p, width = width, height = height)
    message("宿主属堆叠图已成功保存至: ", output_file)
  }
  
  return(p)
}
